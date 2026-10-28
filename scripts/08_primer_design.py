"""Step 08: final qPCR primers (reuse step 06 pairs, design the new genes).

Purpose: for every gene of the final list, (1) re-check the step 06 pair
  (same sequences, same canonical transcript and positions, Primer3
  re-evaluation with the step 06 settings, Ensembl 116 transcriptome check of
  on-target isoforms, off-target, pseudogene and readthrough products), or
  (2) design new pairs with exactly the step 06 code and settings (imported
  from scripts/06_primer_design.py) and pick the best pair with the step 06
  rule. Junction status is recomputed from primer positions on the canonical
  transcript.
Inputs:  results/step08_primers/tables/A_genes_for_primers.csv
         results/step06_qPCR/tables/primer_best.csv (step 06 pairs)
         Ensembl release 116 REST cache and cDNA/ncRNA FASTA (data/reference/ensembl/)
Outputs: results/step08_primers/tables/B_reuse_checks.csv, B_new_candidates.csv,
         B_primers_final.csv, B_primer_map.csv
Date: 2026-10-28
Run from the repo root:
  envs/py310/bin/python scripts/08_primer_design.py
"""
import csv
import importlib.util
import os

import primer3

spec = importlib.util.spec_from_file_location("step06", "scripts/06_primer_design.py")
s6 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s6)

P = s6.r_paths()
TDIR = os.path.join(P["step08"], "tables")
T06 = os.path.join(P["step06"], "tables")
MIN3, MIN5 = (s6.GLOBAL_ARGS["PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION"],
              s6.GLOBAL_ARGS["PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION"])
CHECK_FIELDS = ["main_isoforms_detected", "transcripts_detected", "on_target_sizes", "off_target_products",
                "pseudogene_products", "readthrough_products", "product_bp"]


def junction_ends(exon_lengths):
    out, pos = [], -1
    for L in exon_lengths[:-1]:
        pos += L
        out.append(pos)
    return out


def exon_of(pos, ends):
    return sum(pos > e for e in ends)


def junction_status(p, can, shared):
    """Which primer spans a canonical exon junction (step 06 overlap rules), else exon placement."""
    ends = junction_ends(can["exon_lengths"])
    f0, f1 = p["f_start"], p["f_start"] + p["f_len"] - 1
    r0, r1 = p["r_start"], p["r_start"] + p["r_len"] - 1
    spans = []
    for k, j in enumerate(ends):
        if j - f0 + 1 >= MIN5 and f1 - j >= MIN3:
            spans.append(("F", k))
        if j - r0 + 1 >= MIN3 and r1 - j >= MIN5:     # reverse primer: 3' end at r0
            spans.append(("R", k))
    if spans:
        sh = [s for s in spans if s[1] in shared]
        which = "+".join(sorted({s[0] for s in (sh or spans)}))
        return f"{which} spans an exon junction" + (" shared by all main isoforms" if sh else ""), True
    if exon_of(f0, ends) != exon_of(r1, ends):
        return "primers in different exons (intron-flanking)", False
    if len(can["exon_lengths"]) == 1:
        return "single-exon transcript: DNase and no-RT essential", False
    return "both primers in one exon: DNase and no-RT essential", False


def checks(symbol, p, txs, can, big, starts, names):
    """Transcriptome check of one pair, as in step 06 main()."""
    main_tx = [t for t in txs if t["canonical"] or t["mane"] or (t["primary"] and t["biotype"] == can["biotype"])]
    main_ids = {t["id"] for t in main_tx}
    tx_ids = {t["id"] for t in txs if t["seq"]}
    hits = s6.products(big, starts, names, p["forward"], p["reverse"])
    on = [h for h in hits if h["gene"] == symbol and h["mm"] == 0]
    on_ids = {h["tx"] for h in on}
    same_locus = [h for h in hits if h["gene"] != symbol and "readthrough" in h["desc"]
                  and symbol in h["desc"].replace(",", " ").split()]
    other = [h for h in hits if h["gene"] != symbol and h not in same_locus]
    off = sorted({f"{h['gene']} ({h['biotype']}, {h['size']} bp, {h['mm']} mm)" for h in other})
    rt = sorted({f"{h['gene']} ({h['desc']}, {h['size']} bp)" for h in same_locus})
    pseudo = sorted({h["gene"] for h in other if "pseudogene" in h["biotype"]})
    n_main = len(main_ids & on_ids)
    return {"main_isoforms_detected": f"{n_main}/{len(main_ids)}", "isoform_specific": n_main < len(main_ids),
            "transcripts_detected": f"{len(on_ids & tx_ids)}/{len(tx_ids)}",
            "on_target_sizes": " ".join(str(s) for s in sorted({h["size"] for h in on})),
            "detects_canonical": can["id"] in on_ids,
            "off_target_products": "; ".join(off) if off else "none",
            "pseudogene_products": " ".join(pseudo) if pseudo else "none",
            "readthrough_products": "; ".join(rt) if rt else "none",
            "main_isoform_ids": " ".join(sorted(main_ids)), "detected_transcript_ids": " ".join(sorted(on_ids & tx_ids))}


def primer3_check(seq, f, r):
    """Re-evaluate an existing pair with the step 06 Primer3 settings (task check_primers)."""
    args = dict(s6.GLOBAL_ARGS, PRIMER_TASK="check_primers", PRIMER_PICK_ANYWAY=0)
    args.pop("PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION", None)
    args.pop("PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION", None)
    res = primer3.bindings.design_primers({"SEQUENCE_ID": "t", "SEQUENCE_TEMPLATE": seq, "SEQUENCE_PRIMER": f,
                                           "SEQUENCE_PRIMER_REVCOMP": r}, args)
    ok = res.get("PRIMER_PAIR_NUM_RETURNED", 0) == 1
    return ok, ("pass" if ok else "; ".join(str(res.get(k, "")) for k in
                                            ("PRIMER_LEFT_EXPLAIN", "PRIMER_RIGHT_EXPLAIN", "PRIMER_PAIR_EXPLAIN")))


def oligo_stats(f, r):
    return {"f_tm": round(primer3.calc_tm(f), 1), "r_tm": round(primer3.calc_tm(r), 1),
            "f_gc": round(100 * (f.count("G") + f.count("C")) / len(f), 0),
            "r_gc": round(100 * (r.count("G") + r.count("C")) / len(r), 0),
            "f_hairpin_tm": round(primer3.calc_hairpin(f).tm, 1), "r_hairpin_tm": round(primer3.calc_hairpin(r).tm, 1),
            "f_homodimer_tm": round(primer3.calc_homodimer(f).tm, 1), "r_homodimer_tm": round(primer3.calc_homodimer(r).tm, 1),
            "heterodimer_tm": round(primer3.calc_heterodimer(f, r).tm, 1)}


def main():
    genes = list(csv.DictReader(open(os.path.join(TDIR, "A_genes_for_primers.csv"))))
    p06 = {r["gene"]: r for r in csv.DictReader(open(os.path.join(T06, "primer_best.csv")))}
    print("Loading transcriptome...", flush=True)
    big, starts, names = s6.load_transcriptome()
    print(f"  {len(names):,} transcripts", flush=True)
    final, reuse_rows, new_rows, maps = [], [], [], []
    for g in genes:
        symbol, group, source = g["gene"], g["role"], g["primer_source"]
        gene, txs = s6.transcripts_of(symbol)
        can = next(t for t in txs if t["canonical"])
        main_tx = [t for t in txs if t["canonical"] or t["mane"] or (t["primary"] and t["biotype"] == can["biotype"])]
        cj = s6.junctions_of(can)
        shared = [i for i, j in enumerate(cj) if all(j in s6.junctions_of(t) for t in main_tx)]
        maps.append({"gene": symbol, "kind": "transcript", "start": 0, "end": len(can["seq"]), "label": can["id"]})
        pos = 0
        for i, L in enumerate(can["exon_lengths"]):
            maps.append({"gene": symbol, "kind": "exon", "start": pos, "end": pos + L, "label": str(i + 1)})
            pos += L
        base = {"gene": symbol, "group": group, "primer_source": source, "gene_biotype": gene.get("biotype", ""),
                "transcript": can["id"], "mane": "/".join(can["mane"]), "shared_junctions": f"{len(shared)}/{len(cj)}"}
        if source.startswith("step 06"):
            old = p06[symbol]
            f, r = old["forward"], old["reverse"]
            fs, rs = can["seq"].find(f), can["seq"].find(s6.revcomp(r))
            p = {"design": old["design"], "forward": f, "reverse": r, "f_start": fs, "f_len": len(f),
                 "r_start": rs, "r_len": len(r), "product_bp": rs + len(r) - fs, "penalty": old["penalty"]}
            ck = checks(symbol, p, txs, can, big, starts, names)
            p3_ok, p3_msg = primer3_check(can["seq"], f, r)
            new_vals = {**ck, "product_bp": p["product_bp"]}
            diffs = [k for k in CHECK_FIELDS if str(new_vals[k]) != str(old[k])]
            st = oligo_stats(f, r)
            row = {**base, **p, **ck, **st, "primer3_check": p3_msg,
                   "transcript_unchanged": can["id"] == old["transcript"],
                   "positions_unchanged": fs == int(old["f_start"]) and rs == int(old["r_start"]),
                   "tm_unchanged": abs(st["f_tm"] - float(old["f_tm"])) <= 0.1 and abs(st["r_tm"] - float(old["r_tm"])) <= 0.1,
                   "transcriptome_check_unchanged": not diffs, "changed_fields": "; ".join(diffs) if diffs else "none"}
            row["reuse_status"] = "unchanged, passes" if (row["transcript_unchanged"] and row["positions_unchanged"] and p3_ok
                                                          and row["tm_unchanged"] and not diffs) else "CHECK"
            reuse_rows.append(row)
            chosen = row
        else:
            pairs = s6.design(can["seq"], can["exon_lengths"], shared)
            cands = []
            for rank, p in enumerate(pairs, 1):
                ck = checks(symbol, p, txs, can, big, starts, names)
                cands.append({**base, "pair": rank, **p, **ck, **oligo_stats(p["forward"], p["reverse"])})

            def key(r):   # step 06 selection rule
                n_main = int(r["main_isoforms_detected"].split("/")[0])
                return (r["off_target_products"] != "none", not r["detects_canonical"], -n_main,
                        len(r["on_target_sizes"].split()) > 1, s6.DESIGN_ORDER[r["design"]], r["penalty"])
            cands.sort(key=key)
            for k, c in enumerate(cands):
                c["chosen"] = k == 0
            new_rows += cands
            chosen = cands[0] if cands else None
            print(f"  {symbol:8s} {len(pairs):2d} pairs designed; junctions shared {len(shared)}/{len(cj)}", flush=True)
        if chosen:
            js, spans = junction_status(chosen, can, shared)
            chosen["junction_status"] = js
            chosen["spans_junction"] = spans
            final.append(chosen)
            maps.append({"gene": symbol, "kind": "forward", "start": chosen["f_start"],
                         "end": chosen["f_start"] + chosen["f_len"], "label": chosen["forward"]})
            maps.append({"gene": symbol, "kind": "reverse", "start": chosen["r_start"],
                         "end": chosen["r_start"] + chosen["r_len"], "label": chosen["reverse"]})

    def write(name, rows, cols):
        with open(os.path.join(TDIR, name), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
            w.writeheader()
            w.writerows(rows)

    core = ["gene", "group", "primer_source", "gene_biotype", "transcript", "mane", "design", "junction_status", "spans_junction",
            "forward", "reverse", "f_tm", "r_tm", "f_gc", "r_gc", "f_hairpin_tm", "r_hairpin_tm", "f_homodimer_tm",
            "r_homodimer_tm", "heterodimer_tm", "product_bp", "penalty", "main_isoforms_detected", "isoform_specific",
            "transcripts_detected", "on_target_sizes", "detects_canonical", "off_target_products", "pseudogene_products",
            "readthrough_products", "main_isoform_ids", "detected_transcript_ids", "shared_junctions",
            "f_start", "f_len", "r_start", "r_len"]
    write("B_reuse_checks.csv", reuse_rows, core + ["primer3_check", "transcript_unchanged", "positions_unchanged",
                                                    "tm_unchanged", "transcriptome_check_unchanged", "changed_fields",
                                                    "reuse_status"])
    write("B_new_candidates.csv", new_rows, ["pair", "chosen"] + core)
    write("B_primers_final.csv", final, core + ["reuse_status"])
    write("B_primer_map.csv", maps, ["gene", "kind", "start", "end", "label"])
    print("\nFinal pairs:")
    for r in final:
        print(f"  {r['gene']:8s} {r['product_bp']:4d} bp  {r['junction_status']}; main {r['main_isoforms_detected']}; "
              f"off-target {r['off_target_products']}; {r.get('reuse_status', 'new')}")


if __name__ == "__main__":
    main()
