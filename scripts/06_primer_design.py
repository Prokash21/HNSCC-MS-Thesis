"""Step 06: SYBR Green qRT-PCR primers for the HPV- panel genes.

Purpose: for each gene from scripts/06_qpcr_panel.R (pass 1), design qPCR pairs
  on the Ensembl canonical transcript (70-150 bp product, Tm ~60 C) with one
  primer spanning an exon-exon junction, preferring junctions shared by all
  main isoforms (MANE / GENCODE primary transcripts of the canonical biotype);
  check every pair against the whole human transcriptome (Ensembl cDNA + ncRNA):
  which transcripts it detects, product sizes, and products from other genes
  (paralogues, transcribed pseudogenes). Pairs should still be confirmed in
  NCBI Primer-BLAST before ordering.
Inputs:  results/step06_qPCR/tables/C_genes_for_primers.csv (gene, role)
         Ensembl release 116 REST (gene models, cDNA) and FASTA (cDNA, ncRNA)
Outputs: results/step06_qPCR/tables/primer_candidates.csv (every pair and its checks)
         results/step06_qPCR/tables/primer_best.csv       (best pair per gene)
         results/step06_qPCR/tables/primer_map.csv        (exon / primer coordinates)
Date: 2026-10-27
Run from the repo root (after the first Rscript run of 06_qpcr_panel.R):
  envs/py310/bin/python scripts/06_primer_design.py
"""
import bisect
import csv
import gzip
import hashlib
import json
import os
import re
import time
import urllib.request

import primer3


def r_paths(f="config/paths.R"):
    return dict(re.findall(r'^\s*(\w+)\s*=\s*"([^"]+)"', open(f).read(), flags=re.M))


P = r_paths()
TDIR = os.path.join(P["step06"], "tables")
REF_DIR = P["ensembl_dir"]
REST = "https://rest.ensembl.org"
RELEASE = 116
FASTA = {
    "cdna": f"https://ftp.ensembl.org/pub/release-{RELEASE}/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
    "ncrna": f"https://ftp.ensembl.org/pub/release-{RELEASE}/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz",
}
SEED = 12          # exact match required at the primer's 3' end
MAX_MISMATCH = 4   # over the whole primer, for a site to count as a binding site
MAX_PRODUCT = 2000 # off-target products longer than this are ignored

GLOBAL_ARGS = {
    "PRIMER_TASK": "generic",
    "PRIMER_PICK_LEFT_PRIMER": 1,
    "PRIMER_PICK_RIGHT_PRIMER": 1,
    "PRIMER_PICK_INTERNAL_OLIGO": 0,
    "PRIMER_OPT_SIZE": 20,
    "PRIMER_MIN_SIZE": 18,
    "PRIMER_MAX_SIZE": 24,
    "PRIMER_OPT_TM": 60.0,
    "PRIMER_MIN_TM": 58.0,
    "PRIMER_MAX_TM": 62.0,
    "PRIMER_PAIR_MAX_DIFF_TM": 2.0,
    "PRIMER_MIN_GC": 40.0,
    "PRIMER_MAX_GC": 60.0,
    "PRIMER_MAX_END_GC": 3,
    "PRIMER_MAX_POLY_X": 4,
    "PRIMER_MAX_SELF_ANY_TH": 45.0,
    "PRIMER_MAX_SELF_END_TH": 35.0,
    "PRIMER_PAIR_MAX_COMPL_ANY_TH": 45.0,
    "PRIMER_PAIR_MAX_COMPL_END_TH": 35.0,
    "PRIMER_MAX_HAIRPIN_TH": 24.0,
    "PRIMER_PRODUCT_SIZE_RANGE": [[70, 150]],
    "PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION": 4,
    "PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION": 7,
    "PRIMER_NUM_RETURN": 10,
}
DESIGN_ORDER = {"spans a junction shared by all main isoforms": 0, "spans an exon junction": 1,
                "primers in different exons": 2, "single exon: DNase treatment essential": 3}

COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def rest(path, payload=None, tries=6):
    """Ensembl REST call, cached on disk (the server returns intermittent HTTP 500s)."""
    cache = os.path.join(REF_DIR, f"rest_release{RELEASE}")
    os.makedirs(cache, exist_ok=True)
    f = os.path.join(cache, hashlib.md5((path + json.dumps(payload)).encode()).hexdigest() + ".json")
    if os.path.exists(f):
        return json.load(open(f))
    for i in range(tries):
        try:
            req = urllib.request.Request(
                REST + path,
                data=None if payload is None else json.dumps(payload).encode(),
                headers={"Content-Type": "application/json", "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=120) as r:
                out = json.loads(r.read())
            json.dump(out, open(f, "w"))
            return out
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(5 * (i + 1))


def transcripts_of(symbol):
    """All Ensembl transcripts of a gene, exons in transcript order, with cDNA."""
    gene = rest(f"/lookup/symbol/homo_sapiens/{symbol}?expand=1;mane=1")
    txs = gene["Transcript"]
    seqs = {}
    ids = [t["id"] for t in txs]
    for k in range(0, len(ids), 50):
        for rec in rest("/sequence/id", {"ids": ids[k:k + 50], "type": "cdna"}):
            seqs[rec["id"]] = rec["seq"].upper()
    out = []
    for t in txs:
        exons = sorted(t["Exon"], key=lambda e: e["start"], reverse=t["strand"] == -1)
        out.append({"id": t["id"], "canonical": bool(t.get("is_canonical")), "biotype": t.get("biotype", ""),
                    "mane": [m.get("type") for m in t.get("MANE", [])], "primary": bool(t.get("gencode_primary")),
                    "strand": t["strand"], "exons": [(e["start"], e["end"]) for e in exons],
                    "exon_lengths": [e["end"] - e["start"] + 1 for e in exons], "seq": seqs.get(t["id"], "")})
    return gene, out


def junctions_of(t):
    """Splice junctions as (donor, acceptor) genomic coordinates, in transcript order."""
    ex = t["exons"]
    if t["strand"] == 1:
        return [(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)]
    return [(ex[i][0], ex[i + 1][1]) for i in range(len(ex) - 1)]


def design(seq, exon_lengths, shared):
    """Primer3 pairs from all junction designs, best design first."""
    junctions, pos = [], -1
    for L in exon_lengths[:-1]:
        pos += L
        junctions.append(pos)            # index of the last base of each exon
    ## attempts in order; stop after the first junction-spanning set that gives pairs, but when only
    ## some junctions are shared, pool the shared-junction pairs with the any-junction pairs
    attempts = []
    if junctions:
        sj = [junctions[i] for i in shared]
        if len(sj) == len(junctions):
            attempts.append(("spans a junction shared by all main isoforms", {"SEQUENCE_OVERLAP_JUNCTION_LIST": sj}, True))
        else:
            if sj:
                attempts.append(("spans a junction shared by all main isoforms", {"SEQUENCE_OVERLAP_JUNCTION_LIST": sj}, False))
            attempts.append(("spans an exon junction", {"SEQUENCE_OVERLAP_JUNCTION_LIST": junctions}, True))
        attempts.append(("primers in different exons", {"SEQUENCE_TARGET": [[j, 2] for j in junctions]}, True))
    attempts.append(("single exon: DNase treatment essential", {}, True))
    pairs, seen = [], set()
    for label, extra, stop in attempts:
        seq_args = {"SEQUENCE_ID": "t", "SEQUENCE_TEMPLATE": seq, **extra}
        res = primer3.bindings.design_primers(seq_args, GLOBAL_ARGS)
        for i in range(res.get("PRIMER_PAIR_NUM_RETURNED", 0)):
            f, r = res[f"PRIMER_LEFT_{i}_SEQUENCE"], res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
            if (f, r) in seen:
                continue
            seen.add((f, r))
            lpos, llen = res[f"PRIMER_LEFT_{i}"]
            rpos, rlen = res[f"PRIMER_RIGHT_{i}"]       # rpos = 3'-most base
            pairs.append({
                "design": label, "forward": f, "reverse": r,
                "f_start": lpos, "f_len": llen, "r_start": rpos - rlen + 1, "r_len": rlen,
                "f_tm": round(res[f"PRIMER_LEFT_{i}_TM"], 1), "r_tm": round(res[f"PRIMER_RIGHT_{i}_TM"], 1),
                "f_gc": round(res[f"PRIMER_LEFT_{i}_GC_PERCENT"], 0), "r_gc": round(res[f"PRIMER_RIGHT_{i}_GC_PERCENT"], 0),
                "product_bp": res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"], "penalty": round(res[f"PRIMER_PAIR_{i}_PENALTY"], 2),
            })
        if pairs and stop:
            break
    return pairs


def load_transcriptome():
    os.makedirs(REF_DIR, exist_ok=True)
    names, chunks, starts, offset = [], [], [], 0
    for kind, url in FASTA.items():
        path = os.path.join(REF_DIR, os.path.basename(url))
        if not os.path.exists(path):
            print(f"Downloading Ensembl {kind} FASTA (release {RELEASE})...")
            urllib.request.urlretrieve(url, path)
        head, seq = None, []
        with gzip.open(path, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    if head:
                        names.append(head)
                        s = "".join(seq).upper()
                        starts.append(offset)
                        chunks.append(s)
                        offset += len(s) + 1
                    f = line[1:].split()
                    meta = dict(x.split(":", 1) for x in f[1:] if ":" in x and not x.startswith("chromosome"))
                    desc = line.split("description:", 1)[1].split(" [Source")[0].strip() if "description:" in line else ""
                    head = (f[0].split(".")[0], meta.get("gene_symbol", meta.get("gene", "?")),
                            meta.get("transcript_biotype", ""), desc)
                    seq = []
                else:
                    seq.append(line.strip())
            if head:
                names.append(head)
                s = "".join(seq).upper()
                starts.append(offset)
                chunks.append(s)
                offset += len(s) + 1
    return "#".join(chunks), starts, names


def sites(big, query, seed_at_end):
    """Positions in `big` where `query` binds: exact 3' seed, <= MAX_MISMATCH overall."""
    L = len(query)
    seed = query[-SEED:] if seed_at_end else query[:SEED]
    shift = L - SEED if seed_at_end else 0
    out, i = [], big.find(seed)
    while i != -1:
        s = i - shift
        if s >= 0:
            window = big[s:s + L]
            if "#" not in window:
                mm = sum(a != b for a, b in zip(window, query))
                if mm <= MAX_MISMATCH:
                    out.append((s, mm))
        i = big.find(seed, i + 1)
    return out


def products(big, starts, names, fwd, rev):
    """Transcripts where the pair gives a product (forward on the sense strand,
    reverse primer binding the reverse complement)."""
    f_sites = sites(big, fwd, seed_at_end=True)
    r_sites = sites(big, revcomp(rev), seed_at_end=False)   # rev's 3' end = rc's 5' end
    by_tx = {}
    for pos, mm in f_sites:
        k = bisect.bisect_right(starts, pos) - 1
        by_tx.setdefault(k, {"f": [], "r": []})["f"].append((pos, mm))
    for pos, mm in r_sites:
        k = bisect.bisect_right(starts, pos) - 1
        by_tx.setdefault(k, {"f": [], "r": []})["r"].append((pos, mm))
    hits = []
    for k, d in by_tx.items():
        for fp, fmm in d["f"]:
            for rp, rmm in d["r"]:
                size = rp + len(rev) - fp
                if 0 < size <= MAX_PRODUCT:
                    hits.append({"tx": names[k][0], "gene": names[k][1], "biotype": names[k][2], "size": size,
                                 "mm": fmm + rmm, "desc": names[k][3]})
    return hits


def main():
    genes = list(csv.DictReader(open(os.path.join(TDIR, "C_genes_for_primers.csv"))))
    print("Loading transcriptome...", flush=True)
    big, starts, names = load_transcriptome()
    print(f"  {len(names):,} transcripts", flush=True)
    rows, maps = [], []
    for g in genes:
        symbol, role = g["gene"], g["role"]
        gene, txs = transcripts_of(symbol)
        can = next(t for t in txs if t["canonical"])
        if sum(can["exon_lengths"]) != len(can["seq"]):
            print(f"  warning: {symbol} exon lengths do not add up to the cDNA length")
        ## main isoforms: canonical, MANE, and GENCODE primary transcripts of the canonical biotype
        main_tx = [t for t in txs if t["canonical"] or t["mane"] or (t["primary"] and t["biotype"] == can["biotype"])]
        main_ids = {t["id"] for t in main_tx}
        cj = junctions_of(can)
        shared = [i for i, j in enumerate(cj) if all(j in junctions_of(t) for t in main_tx)]
        pairs = design(can["seq"], can["exon_lengths"], shared)
        tx_ids = {t["id"] for t in txs if t["seq"]}
        maps.append({"gene": symbol, "kind": "transcript", "start": 0, "end": len(can["seq"]), "label": can["id"]})
        pos = 0
        for i, L in enumerate(can["exon_lengths"]):
            maps.append({"gene": symbol, "kind": "exon", "start": pos, "end": pos + L, "label": str(i + 1)})
            pos += L
        for rank, p in enumerate(pairs, 1):
            hits = products(big, starts, names, p["forward"], p["reverse"])
            on = [h for h in hits if h["gene"] == symbol and h["mm"] == 0]
            on_ids = {h["tx"] for h in on}
            ## readthrough genes that contain the target gene's own exons are the same locus, not off-target
            same_locus = [h for h in hits if h["gene"] != symbol and "readthrough" in h["desc"]
                          and symbol in h["desc"].replace(",", " ").split()]
            other = [h for h in hits if h["gene"] != symbol and h not in same_locus]
            off = sorted({f"{h['gene']} ({h['biotype']}, {h['size']} bp, {h['mm']} mm)" for h in other})
            readthrough = sorted({f"{h['gene']} ({h['desc']}, {h['size']} bp)" for h in same_locus})
            pseudo = sorted({h["gene"] for h in other if "pseudogene" in h["biotype"]})
            n_main = len(main_ids & on_ids)
            rows.append({"gene": symbol, "role": role, "gene_biotype": gene.get("biotype", ""),
                         "transcript": can["id"], "mane": "/".join(can["mane"]), "pair": rank, **p,
                         "main_isoforms": len(main_ids), "main_isoforms_detected": f"{n_main}/{len(main_ids)}",
                         "isoform_specific": n_main < len(main_ids),
                         "transcripts_detected": f"{len(on_ids & tx_ids)}/{len(tx_ids)}",
                         "on_target_sizes": " ".join(str(s) for s in sorted({h["size"] for h in on})),
                         "detects_canonical": can["id"] in on_ids,
                         "off_target_products": "; ".join(off) if off else "none",
                         "pseudogene_products": " ".join(pseudo) if pseudo else "none",
                         "readthrough_products": "; ".join(readthrough) if readthrough else "none",
                         "main_isoform_ids": " ".join(sorted(main_ids)),
                         "detected_transcript_ids": " ".join(sorted(on_ids & tx_ids)),
                         "shared_junctions": f"{len(shared)}/{len(cj)}"})
        print(f"  {symbol:10s} {len(pairs):2d} pairs; main isoforms {len(main_ids)}; "
              f"junctions shared {len(shared)}/{len(cj)}", flush=True)

    ## best pair per gene: no off-target product, detects the canonical transcript,
    ## most main isoforms, one product size, junction design, then Primer3 penalty
    def key(r):
        n_main = int(r["main_isoforms_detected"].split("/")[0])
        return (r["off_target_products"] != "none", not r["detects_canonical"], -n_main,
                len(r["on_target_sizes"].split()) > 1, DESIGN_ORDER[r["design"]], r["penalty"])
    best = {}
    for r in sorted(rows, key=key):
        best.setdefault(r["gene"], r)
    for r in rows:
        r["chosen"] = best.get(r["gene"]) is r

    cols = ["gene", "role", "gene_biotype", "transcript", "mane", "pair", "chosen", "design", "forward", "reverse",
            "f_tm", "r_tm", "f_gc", "r_gc", "product_bp", "penalty", "main_isoforms_detected", "isoform_specific",
            "transcripts_detected", "on_target_sizes", "detects_canonical", "off_target_products", "pseudogene_products",
            "readthrough_products", "main_isoform_ids", "detected_transcript_ids", "shared_junctions",
            "f_start", "f_len", "r_start", "r_len"]
    with open(os.path.join(TDIR, "primer_candidates.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(TDIR, "primer_best.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for g in genes:
            if g["gene"] in best:
                w.writerow(best[g["gene"]])
    for g in genes:
        r = best.get(g["gene"])
        if r:
            maps.append({"gene": g["gene"], "kind": "forward", "start": r["f_start"],
                         "end": r["f_start"] + r["f_len"], "label": r["forward"]})
            maps.append({"gene": g["gene"], "kind": "reverse", "start": r["r_start"],
                         "end": r["r_start"] + r["r_len"], "label": r["reverse"]})
    with open(os.path.join(TDIR, "primer_map.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["gene", "kind", "start", "end", "label"])
        w.writeheader()
        w.writerows(maps)

    print("\nBest pair per gene:")
    for g in genes:
        r = best.get(g["gene"])
        if r:
            print(f"  {g['gene']:10s} {r['product_bp']:4d} bp  {r['design']}; main {r['main_isoforms_detected']}; "
                  f"off-target: {r['off_target_products']}")
        else:
            print(f"  {g['gene']:10s} no pair found")


if __name__ == "__main__":
    main()
