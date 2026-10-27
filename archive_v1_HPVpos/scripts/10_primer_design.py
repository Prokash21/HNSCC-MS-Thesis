"""Design SYBR Green qRT-PCR primers for the panel from scripts/08.

For each gene: take the Ensembl canonical transcript and ask Primer3 for
qPCR pairs (70-150 bp product, Tm ~60 C) where one primer spans an exon-exon
junction, so genomic DNA is not amplified. Every pair is then checked against
the whole human transcriptome (Ensembl cDNA + ncRNA): which transcripts of the
gene it detects, and whether it could also amplify another gene. Confirm the
chosen pairs in NCBI Primer-BLAST before ordering.

Outputs: HNSCC-HPV/results/05_primers/tables/
           primer_candidates.csv  every designed pair with its checks
           primer_order_sheet.csv the best pair per gene
           primer_map.csv         exon / primer coordinates for the map figure
         data/reference/ensembl/  cached transcriptome FASTA (not in git)

Run from the repo root:
  uv run --with primer3-py python scripts/10_primer_design.py
"""
import bisect
import csv
import gzip
import json
import os
import time
import urllib.request

import primer3

GENES = [
    ("SYCP2", "Panel: positive control"),
    ("ZFR2", "Panel: novelty"),
    ("MIR9-3HG", "Panel: novelty (lncRNA)"),
    ("KLHL35", "Panel: novelty"),
    ("GZMK", "TME: immune marker"),
    ("COL1A2", "TME: stromal marker"),
    ("SFRP2", "TME: stromal marker"),
    ("RPS18", "Reference gene"),
    ("SDHA", "Reference gene"),
    ("PSMC5", "Reference gene"),
    ("RPLP0", "Reference gene (backup)"),
    ("FOXD3-AS1", "Backup: lncRNA"),
    ("RIBC2", "Backup: most novel, low expression"),
]
OUT_DIR = "HNSCC-HPV/results/05_primers/tables"
REF_DIR = "data/reference/ensembl"
REST = "https://rest.ensembl.org"
FASTA = {
    "cdna": "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
    "ncrna": "https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz",
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
    "PRIMER_NUM_RETURN": 8,
}

COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def rest(path, payload=None, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(
                REST + path,
                data=None if payload is None else json.dumps(payload).encode(),
                headers={"Content-Type": "application/json", "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read())
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(2 * (i + 1))


def transcripts_of(symbol):
    """All Ensembl transcripts of a gene, exons in transcript order, with cDNA."""
    gene = rest(f"/lookup/symbol/homo_sapiens/{symbol}?expand=1")
    txs = gene["Transcript"]
    seqs = {}
    ids = [t["id"] for t in txs]
    for k in range(0, len(ids), 50):
        for rec in rest("/sequence/id", {"ids": ids[k:k + 50], "type": "cdna"}):
            seqs[rec["id"]] = rec["seq"].upper()
    out = []
    for t in txs:
        exons = sorted(t["Exon"], key=lambda e: e["start"], reverse=t["strand"] == -1)
        out.append({"id": t["id"], "canonical": bool(t.get("is_canonical")),
                    "biotype": t.get("biotype", ""), "exon_lengths": [e["end"] - e["start"] + 1 for e in exons],
                    "seq": seqs.get(t["id"], "")})
    return gene, out


def design(seq, exon_lengths):
    """Primer3 pairs; prefer a primer spanning a junction, then intron-flanking."""
    junctions, pos = [], -1
    for L in exon_lengths[:-1]:
        pos += L
        junctions.append(pos)            # index of the last base of each exon
    attempts = []
    if junctions:
        attempts.append(("primer spans exon junction", {"SEQUENCE_OVERLAP_JUNCTION_LIST": junctions}))
        attempts.append(("primers in different exons", {"SEQUENCE_TARGET": [[j, 2] for j in junctions]}))
    attempts.append(("single exon: DNase treatment essential", {}))
    for label, extra in attempts:
        seq_args = {"SEQUENCE_ID": "t", "SEQUENCE_TEMPLATE": seq, **extra}
        res = primer3.bindings.design_primers(seq_args, GLOBAL_ARGS)
        n = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
        if n:
            pairs = []
            for i in range(n):
                lpos, llen = res[f"PRIMER_LEFT_{i}"]
                rpos, rlen = res[f"PRIMER_RIGHT_{i}"]       # rpos = 3'-most base
                pairs.append({
                    "design": label,
                    "forward": res[f"PRIMER_LEFT_{i}_SEQUENCE"],
                    "reverse": res[f"PRIMER_RIGHT_{i}_SEQUENCE"],
                    "f_start": lpos, "f_len": llen, "r_start": rpos - rlen + 1, "r_len": rlen,
                    "f_tm": round(res[f"PRIMER_LEFT_{i}_TM"], 1),
                    "r_tm": round(res[f"PRIMER_RIGHT_{i}_TM"], 1),
                    "f_gc": round(res[f"PRIMER_LEFT_{i}_GC_PERCENT"], 0),
                    "r_gc": round(res[f"PRIMER_RIGHT_{i}_GC_PERCENT"], 0),
                    "product_bp": res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
                    "penalty": round(res[f"PRIMER_PAIR_{i}_PENALTY"], 2),
                })
            return pairs, junctions
    return [], junctions


def load_transcriptome():
    os.makedirs(REF_DIR, exist_ok=True)
    names, chunks, starts, offset = [], [], [], 0
    for kind, url in FASTA.items():
        path = os.path.join(REF_DIR, os.path.basename(url))
        if not os.path.exists(path):
            print(f"Downloading Ensembl {kind} FASTA...")
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
                    head = (f[0].split(".")[0], meta.get("gene_symbol", meta.get("gene", "?")))
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
                    hits.append({"tx": names[k][0], "gene": names[k][1], "size": size, "mm": fmm + rmm})
    return hits


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print("Loading transcriptome...")
    big, starts, names = load_transcriptome()
    print(f"  {len(names):,} transcripts")
    rows, maps = [], []
    for symbol, role in GENES:
        gene, txs = transcripts_of(symbol)
        can = next(t for t in txs if t["canonical"])
        if sum(can["exon_lengths"]) != len(can["seq"]):
            print(f"  warning: {symbol} exon lengths do not add up to the cDNA length")
        pairs, junctions = design(can["seq"], can["exon_lengths"])
        tx_ids = {t["id"] for t in txs if t["seq"]}
        maps.append({"gene": symbol, "kind": "transcript", "start": 0, "end": len(can["seq"]),
                     "label": can["id"]})
        pos = 0
        for i, L in enumerate(can["exon_lengths"]):
            maps.append({"gene": symbol, "kind": "exon", "start": pos, "end": pos + L, "label": str(i + 1)})
            pos += L
        for rank, p in enumerate(pairs, 1):
            hits = products(big, starts, names, p["forward"], p["reverse"])
            on = sorted({h["tx"] for h in hits if h["gene"] == symbol and h["mm"] == 0})
            off = sorted({f"{h['gene']} ({h['size']} bp, {h['mm']} mm)" for h in hits if h["gene"] != symbol})
            rows.append({"gene": symbol, "role": role, "gene_biotype": gene.get("biotype", ""),
                         "transcript": can["id"], "pair": rank, **p,
                         "transcripts_detected": f"{len([t for t in on if t in tx_ids])}/{len(tx_ids)}",
                         "detects_canonical": can["id"] in on,
                         "off_target_products": "; ".join(off) if off else "none"})
        print(f"  {symbol:10s} {len(pairs)} pairs ({pairs[0]['design'] if pairs else 'none'})")

    ## best pair per gene: no off-target product, detects the canonical
    ## transcript, most transcripts detected, then Primer3 penalty
    def key(r):
        n_det = int(r["transcripts_detected"].split("/")[0])
        return (r["off_target_products"] != "none", not r["detects_canonical"], -n_det, r["penalty"])
    best = {}
    for r in sorted(rows, key=key):
        best.setdefault(r["gene"], r)
    for r in rows:
        r["chosen"] = best.get(r["gene"]) is r

    cols = ["gene", "role", "gene_biotype", "transcript", "pair", "chosen", "design", "forward", "reverse",
            "f_tm", "r_tm", "f_gc", "r_gc", "product_bp", "penalty", "transcripts_detected",
            "detects_canonical", "off_target_products", "f_start", "f_len", "r_start", "r_len"]
    with open(os.path.join(OUT_DIR, "primer_candidates.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    order_cols = ["gene", "role", "forward", "reverse", "f_tm", "r_tm", "f_gc", "r_gc", "product_bp",
                  "design", "transcript", "transcripts_detected", "off_target_products"]
    with open(os.path.join(OUT_DIR, "primer_order_sheet.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=order_cols, extrasaction="ignore")
        w.writeheader()
        for symbol, _ in GENES:
            if symbol in best:
                w.writerow(best[symbol])
    for symbol, _ in GENES:
        r = best.get(symbol)
        if r:
            maps.append({"gene": symbol, "kind": "forward", "start": r["f_start"],
                         "end": r["f_start"] + r["f_len"], "label": r["forward"]})
            maps.append({"gene": symbol, "kind": "reverse", "start": r["r_start"],
                         "end": r["r_start"] + r["r_len"], "label": r["reverse"]})
    with open(os.path.join(OUT_DIR, "primer_map.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["gene", "kind", "start", "end", "label"])
        w.writeheader()
        w.writerows(maps)

    print("\nBest pair per gene:")
    for symbol, _ in GENES:
        r = best.get(symbol)
        if r:
            print(f"  {symbol:10s} F {r['forward']:24s} R {r['reverse']:24s} {r['product_bp']:4d} bp  "
                  f"Tm {r['f_tm']}/{r['r_tm']}  {r['design']}; detects {r['transcripts_detected']}; "
                  f"off-target: {r['off_target_products']}")
        else:
            print(f"  {symbol:10s} no pair found")


if __name__ == "__main__":
    main()
