"""Step 05 (Python parts): pySCENIC gene regulatory networks, HPV- epithelium.

Purpose: infer regulons (GRNBoost2 -> cisTarget motif pruning -> AUCell) in
  HPV- malignant and squamous normal epithelial cells; compare regulon
  activity between malignant and normal cells per patient; regulon
  specificity scores; stability over three runs (different seeds and
  malignant-cell subsamples); AUCell and tests use consensus regulons
  (regulon and target found in >= 2 of 3 runs, >= 10 targets).
Inputs:  data/processed/v2/step05_export/ (matrix.mtx genes x cells, raw counts,
         genes.csv, cells.csv, candidate_genes.csv; from scripts/05_TF_activity.R)
         cisTarget hg38 v10 resources and the hs_hgnc TF list (downloaded)
Outputs: data/processed/v2/step05_scenic/ (adjacencies, regulons, AUC matrix)
         results/step05_TF/tables/D_*.csv
Date: 2026-10-27
Run from the repo root (after the first Rscript run of 05_TF_activity.R):
  envs/py310/bin/python scripts/05_pyscenic.py
  STEP05_TIMING=50 envs/py310/bin/python scripts/05_pyscenic.py   (time 50 targets only)
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
import re
import sys
import time
import resource
import urllib.request
import multiprocessing as mp

import numpy as np
import pandas as pd
import scipy.io
from scipy import stats

N_WORKERS = 9
CAP_PER_PATIENT = 500
MIN_FRAC = 0.01
SEEDS = [1, 2, 3]
MIN_MAL, MIN_NORM = 20, 10
MIN_RUNS, MIN_TARGETS = 2, 10        # consensus regulons: found in >= 2 of 3 runs, >= 10 consensus targets
MAX_DOWNLOAD_GB = 4.0
CISTARGET = {
    "rank_500bp": "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/"
                  "gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
    "rank_10kb": "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/"
                 "gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
    "motifs": "https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
    "tfs": "https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/hs_hgnc_tfs.txt",
}


def r_paths(f="config/paths.R"):
    """Paths from config/paths.R (name = "value" lines)."""
    return dict(re.findall(r'^\s*(\w+)\s*=\s*"([^"]+)"', open(f).read(), flags=re.M))


P = r_paths()
EXPORT, SCENIC, CIS = P["step05_export"], P["step05_scenic"], P["cistarget_dir"]
TDIR = os.path.join(P["step05"], "tables")
for d in (SCENIC, CIS, TDIR):
    os.makedirs(d, exist_ok=True)
RUNTIME = []


def peak_mb():
    """Peak resident memory of this process and of finished children (MB; macOS reports bytes)."""
    s = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    c = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return round(s / 1e6), round(c / 1e6)


def log_time(stage, t0):
    s, c = peak_mb()
    RUNTIME.append({"stage": stage, "seconds": round(time.time() - t0, 1), "peak_main_MB": s, "peak_worker_MB": c})
    print(f"{stage}: {RUNTIME[-1]['seconds']} s (peak main {s} MB, worker {c} MB)", flush=True)


def bh(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    if ok.any():
        from statsmodels.stats.multitest import multipletests
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


# --------------------------------------------------------------------------
# D2. cisTarget resources (sizes checked before downloading)
# --------------------------------------------------------------------------
def fetch_resources():
    files = {k: os.path.join(CIS, os.path.basename(u)) for k, u in CISTARGET.items()}
    missing = {k: u for k, u in CISTARGET.items() if not os.path.exists(files[k])}
    sizes = {}
    for k, u in CISTARGET.items():
        req = urllib.request.Request(u, method="HEAD")
        with urllib.request.urlopen(req) as r:
            sizes[k] = int(r.headers.get("Content-Length", 0))
    total_gb = sum(sizes.values()) / 1e9
    pd.DataFrame({"resource": list(sizes), "url": [CISTARGET[k] for k in sizes], "bytes": list(sizes.values())}) \
        .to_csv(os.path.join(TDIR, "D_cistarget_resources.csv"), index=False)
    print(f"cisTarget resources: {total_gb:.2f} GB in total", flush=True)
    if total_gb > MAX_DOWNLOAD_GB:
        sys.exit(f"cisTarget resources {total_gb:.2f} GB exceed {MAX_DOWNLOAD_GB} GB; stopping")
    for k, u in missing.items():
        print(f"downloading {os.path.basename(u)}", flush=True)
        urllib.request.urlretrieve(u, files[k] + ".part")
        if os.path.getsize(files[k] + ".part") != sizes[k]:
            sys.exit(f"incomplete download: {u}")
        os.replace(files[k] + ".part", files[k])
    return files


# --------------------------------------------------------------------------
# GRNBoost2 with a process pool (as pySCENIC's arboreto_with_multiprocessing)
# --------------------------------------------------------------------------
_G = {}


def _init(ex, tf_mat, tf_names, genes, seed):
    _G.update(ex=ex, tf_mat=tf_mat, tf_names=tf_names, genes=genes, seed=seed)


def _one_target(j):
    from arboreto.core import infer_partial_network, SGBM_KWARGS, EARLY_STOP_WINDOW_LENGTH
    return infer_partial_network("GBM", SGBM_KWARGS, _G["tf_mat"], _G["tf_names"], _G["genes"][j],
                                 _G["ex"][:, j], False, EARLY_STOP_WINDOW_LENGTH, _G["seed"])


def grnboost2(ex, genes, tfs, seed, targets=None):
    tf_idx = [i for i, g in enumerate(genes) if g in tfs]
    tf_mat, tf_names = ex[:, tf_idx], [genes[i] for i in tf_idx]
    targets = range(len(genes)) if targets is None else targets
    _init(ex, tf_mat, tf_names, genes, seed)        # inherited by forked workers
    with mp.get_context("fork").Pool(N_WORKERS) as pool:
        res = list(pool.imap_unordered(_one_target, targets, chunksize=8))
    adj = pd.concat(res, ignore_index=True)
    return adj.sort_values("importance", ascending=False).reset_index(drop=True)


def main():
    t_all = time.time()
    # ---- inputs ------------------------------------------------------------
    t0 = time.time()
    genes_all = pd.read_csv(os.path.join(EXPORT, "genes.csv"))["gene"].tolist()
    cells = pd.read_csv(os.path.join(EXPORT, "cells.csv"))
    cand = pd.read_csv(os.path.join(EXPORT, "candidate_genes.csv"))["gene"].tolist()
    X = scipy.io.mmread(os.path.join(EXPORT, "matrix.mtx")).tocsr().T.tocsr()   # cells x genes
    assert X.shape == (len(cells), len(genes_all))
    log_time("load export", t0)

    files = fetch_resources()
    tfs = set(open(files["tfs"]).read().split())

    # ---- cell selection per run and gene set (defined on run 1) -------------
    def select(seed):
        rng = np.random.default_rng(seed)
        keep = []
        for (state, patient), idx in cells.groupby(["state", "patient"]).groups.items():
            idx = np.asarray(idx)
            if state == "malignant" and len(idx) > CAP_PER_PATIENT:
                idx = rng.choice(idx, CAP_PER_PATIENT, replace=False)
            keep.extend(idx.tolist())
        return np.sort(np.array(keep))

    sel = {s: select(s) for s in SEEDS}
    det = np.asarray((X[sel[SEEDS[0]]] > 0).mean(axis=0)).ravel()
    tot = np.asarray(X.sum(axis=0)).ravel()
    gmask = (det >= MIN_FRAC) | ((np.isin(genes_all, list(tfs)) | np.isin(genes_all, cand)) & (tot > 0))
    genes = [g for g, m in zip(genes_all, gmask) if m]
    gi = np.where(gmask)[0]
    lib = cells["n_counts"].to_numpy(float)

    def norm(rows):
        m = X[rows][:, gi].toarray().astype(np.float32)
        return np.log1p(m / lib[rows, None] * 1e4).astype(np.float32)

    sel_tab = pd.DataFrame([{"run": i + 1, "seed": s, "cells": len(sel[s]),
                             "malignant": int((cells.loc[sel[s], "state"] == "malignant").sum()),
                             "normal": int((cells.loc[sel[s], "state"] == "normal").sum()),
                             "patients": cells.loc[sel[s], "patient"].nunique()} for i, s in enumerate(SEEDS)])
    sel_tab["genes"] = len(genes)
    sel_tab["TFs_in_genes"] = len(set(genes) & tfs)
    sel_tab["candidates_in_genes"] = len(set(genes) & set(cand))
    sel_tab.to_csv(os.path.join(TDIR, "D_input_runs.csv"), index=False)
    print(sel_tab, flush=True)

    timing = int(os.environ.get("STEP05_TIMING", "0"))
    if timing:
        ex = norm(sel[SEEDS[0]])
        t0 = time.time()
        rng = np.random.default_rng(0)
        grnboost2(ex, genes, tfs, SEEDS[0], targets=rng.choice(len(genes), timing, replace=False).tolist())
        dt = time.time() - t0
        print(f"timing: {timing} targets in {dt:.1f} s -> full run ~{dt / timing * len(genes) / 60:.0f} min per seed", flush=True)
        return

    # ---- D3. GRNBoost2 -> ctx -> regulons, per run -------------------------
    from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
    from pyscenic.utils import modules_from_adjacencies
    from pyscenic.prune import prune2df, df2regulons
    dbs = [RankingDatabase(fname=files[k], name=k) for k in ("rank_500bp", "rank_10kb")]
    db_genes = set(dbs[0].genes)
    pd.DataFrame([{"genes_in_input": len(genes), "genes_in_ranking_db": len(set(genes) & db_genes)}]) \
        .to_csv(os.path.join(TDIR, "D_genes_in_cistarget.csv"), index=False)
    regulons_by_run = {}
    for run, seed in enumerate(SEEDS, start=1):
        reg_f = os.path.join(SCENIC, f"regulons_run{run}.csv")
        if os.path.exists(reg_f):
            regulons_by_run[run] = pd.read_csv(reg_f)
            continue
        ex = norm(sel[seed])
        adj_f = os.path.join(SCENIC, f"adjacencies_run{run}.tsv.gz")
        if os.path.exists(adj_f):
            adj = pd.read_csv(adj_f, sep="\t")
        else:
            t0 = time.time()
            adj = grnboost2(ex, genes, tfs, seed)
            adj.to_csv(adj_f, sep="\t", index=False)
            log_time(f"GRNBoost2 run {run} ({len(sel[seed])} cells, {len(genes)} genes)", t0)
        t0 = time.time()
        ex_df = pd.DataFrame(ex, columns=genes)
        modules = list(modules_from_adjacencies(adj, ex_df))
        df = prune2df(dbs, modules, files["motifs"], num_workers=N_WORKERS, client_or_address="custom_multiprocessing")
        regs = df2regulons(df)
        df.to_csv(os.path.join(SCENIC, f"ctx_motifs_run{run}.csv.gz"))
        rows = [{"regulon": r.name, "TF": r.transcription_factor, "target": g, "weight": w}
                for r in regs for g, w in r.gene2weight.items()]
        regulons_by_run[run] = pd.DataFrame(rows)
        regulons_by_run[run].to_csv(reg_f, index=False)
        log_time(f"ctx run {run} ({len(modules)} modules -> {len(regs)} regulons)", t0)
        del ex, ex_df, modules, df

    # ---- D4. consensus regulons over the three runs, AUCell on all exported cells
    # regulon kept if found in >= MIN_RUNS runs; target kept if found in >= MIN_RUNS runs;
    # regulon scored if it keeps >= MIN_TARGETS targets; weight = mean over runs
    from ctxcore.genesig import Regulon
    from pyscenic.aucell import aucell
    allr = pd.concat([df.assign(run=run) for run, df in regulons_by_run.items()], ignore_index=True)
    reg_runs = allr.groupby("regulon")["run"].nunique()
    cons = allr[allr["regulon"].isin(reg_runs[reg_runs >= MIN_RUNS].index)]
    cons = cons.groupby(["regulon", "TF", "target"]).agg(weight=("weight", "mean"), target_runs=("run", "nunique")).reset_index()
    cons = cons[cons["target_runs"] >= MIN_RUNS]
    size = cons.groupby("regulon")["target"].nunique()
    cons = cons[cons["regulon"].isin(size[size >= MIN_TARGETS].index)]
    cons.to_csv(os.path.join(SCENIC, "regulons_consensus.csv"), index=False)
    r1 = cons
    auc_f = os.path.join(SCENIC, "auc_consensus.csv.gz")
    if os.path.exists(auc_f):
        auc = pd.read_csv(auc_f, index_col=0)
    else:
        t0 = time.time()
        regs = [Regulon(name=n, gene2weight=dict(zip(d["target"], d["weight"])), transcription_factor=d["TF"].iloc[0],
                        gene2occurrence={}) for n, d in r1.groupby("regulon")]
        ex_all = pd.DataFrame(norm(np.arange(len(cells))), columns=genes, index=cells["cell"])
        auc = aucell(ex_all, regs, num_workers=N_WORKERS, seed=SEEDS[0])
        auc.to_csv(auc_f)
        log_time(f"AUCell ({auc.shape[0]} cells x {auc.shape[1]} regulons)", t0)
        del ex_all

    # ---- per-patient activity: paired (pairs) and all-patient tests ---------
    meta = cells.set_index("cell").loc[auc.index]
    long = auc.copy()
    long["patient"], long["state"] = meta["patient"].values, meta["state"].values
    n_cells = long.groupby(["patient", "state"]).size().rename("n_cells")
    pm = long.groupby(["patient", "state"]).mean()
    ok = n_cells[((n_cells.index.get_level_values("state") == "malignant") & (n_cells >= MIN_MAL)) |
                 ((n_cells.index.get_level_values("state") == "normal") & (n_cells >= MIN_NORM))].index
    pm = pm.loc[ok]
    pm_long = pm.reset_index().melt(id_vars=["patient", "state"], var_name="regulon", value_name="mean_auc")
    pm_long = pm_long.merge(n_cells.reset_index(), on=["patient", "state"])
    pm_long.to_csv(os.path.join(TDIR, "D_regulon_auc_per_patient.csv"), index=False)
    mal = pm.xs("malignant", level="state")
    nor = pm.xs("normal", level="state")
    pairs = sorted(set(mal.index) & set(nor.index))
    rows = []
    for rg in auc.columns:
        d = mal.loc[pairs, rg].values - nor.loc[pairs, rg].values
        p_pair = stats.wilcoxon(d).pvalue if np.any(d != 0) else np.nan
        p_all = stats.mannwhitneyu(mal[rg].values, nor[rg].values, alternative="two-sided").pvalue
        rows.append({"regulon": rg, "paired_diff": d.mean(), "paired_p": p_pair,
                     "all_diff": mal[rg].median() - nor[rg].median(), "all_p": p_all})
    ra = pd.DataFrame(rows)
    ra["paired_FDR"], ra["all_FDR"] = bh(ra["paired_p"]), bh(ra["all_p"])
    ra["TF"] = ra["regulon"].str.replace(r"\(.*\)$", "", regex=True)
    ra["n_targets"] = ra["regulon"].map(r1.groupby("regulon")["target"].nunique())

    from pyscenic.rss import regulon_specificity_scores
    rss = regulon_specificity_scores(auc, meta["state"])
    ra["RSS_malignant"] = ra["regulon"].map(rss.loc["malignant"])
    ra["RSS_normal"] = ra["regulon"].map(rss.loc["normal"])

    # motif support (all runs): number of enriched motifs, best NES, contexts
    ctx = pd.concat([pd.read_csv(os.path.join(SCENIC, f"ctx_motifs_run{run}.csv.gz"), header=[0, 1], index_col=[0, 1])
                     for run in regulons_by_run])
    ctx = ctx.sort_values(("Enrichment", "NES"), ascending=False)
    ctx_tf = ctx.index.get_level_values(0)
    sup = pd.DataFrame({"TF": ctx_tf, "NES": ctx[("Enrichment", "NES")].values,
                        "context": ctx[("Enrichment", "Context")].astype(str).values,
                        "motif": ctx.index.get_level_values(1)})
    sup = sup.groupby("TF").agg(n_motifs=("motif", "nunique"), best_NES=("NES", "max"),
                                best_motif=("motif", lambda m: m.iloc[0]),
                                databases=("context", lambda c: " ".join(sorted({k for k in ("rank_500bp", "rank_10kb") if any(k in x for x in c)}))))
    ra = ra.merge(sup, left_on="TF", right_index=True, how="left")

    # ---- D6. stability over the three runs ----------------------------------
    sets = {run: df.groupby("regulon")["target"].apply(set).to_dict() for run, df in regulons_by_run.items()}
    stab = []
    for rg in sorted(set().union(*[set(s) for s in sets.values()])):
        found = [run for run in sets if rg in sets[run]]
        jac = [len(sets[a][rg] & sets[b][rg]) / len(sets[a][rg] | sets[b][rg])
               for i, a in enumerate(found) for b in found[i + 1:]]
        stab.append({"regulon": rg, "runs_found": len(found), "runs": " ".join(map(str, found)),
                     "mean_target_jaccard": np.mean(jac) if jac else np.nan})
    stab = pd.DataFrame(stab)
    stab.to_csv(os.path.join(TDIR, "D_regulon_stability.csv"), index=False)
    ra = ra.merge(stab[["regulon", "runs_found", "mean_target_jaccard"]], on="regulon", how="left")
    ra["runs_found"] = ra["runs_found"].fillna(1).astype(int)

    cand_set = set(cand)
    ra["candidate_targets"] = ra["regulon"].map(
        r1.groupby("regulon")["target"].apply(lambda t: " ".join(sorted(set(t) & cand_set))))
    ra = ra.sort_values("all_p")
    ra.to_csv(os.path.join(TDIR, "D_regulon_activity.csv"), index=False)
    pd.DataFrame([{"run": run, "regulons": df["regulon"].nunique(), "targets_median": df.groupby("regulon")["target"].nunique().median()}
                  for run, df in regulons_by_run.items()]).to_csv(os.path.join(TDIR, "D_regulons_per_run.csv"), index=False)
    log_time("total", t_all)
    rt = pd.DataFrame(RUNTIME)
    rt_f = os.path.join(TDIR, "D_runtime.csv")
    if os.path.exists(rt_f):
        rt = pd.concat([pd.read_csv(rt_f), rt], ignore_index=True)
    rt.to_csv(rt_f, index=False)
    print(ra.head(25).to_string(), flush=True)


if __name__ == "__main__":
    main()
