# Script for selecting better representative genes from the mmseqs2 clustering (no partial genes)

#!/usr/bin/env python3
from collections import defaultdict

# --- Snakemake IO ---
tsv_file = snakemake.input['tsv']
genes_file = snakemake.input['genes']
lengths_file = snakemake.input['lengths']
nr_genes_file = snakemake.output['nr_genes']
corrected_tsv_file = snakemake.output['corrected_tsv']

# --- 1. Load gene lengths ---
lengths = {}
with open(lengths_file) as f:
    for line in f:
        if not line.strip():
            continue
        gene, length = line.strip().split("\t")[:2]
        lengths[gene] = int(length)

# --- 2. Parse cluster table ---
clusters = defaultdict(list)
with open(tsv_file) as f:
    for line in f:
        if not line.strip():
            continue
        parts = line.strip().split("\t")
        if len(parts) < 2:
            continue
        cluster, gene = parts[0], parts[1]
        clusters[cluster].append(gene)

# --- 3. Refine representatives (exclude short/partial genes) ---
corrected = []
for cluster, genes in clusters.items():
    valid_genes = [g for g in genes if g in lengths]
    valid_genes.sort(key=lambda g: lengths[g])
    if not valid_genes:
        continue
    cutoff = int(len(valid_genes) * 0.75)  # keep genes ≥75% length rank in cluster
    for i, g in enumerate(valid_genes, 1):
        if i >= cutoff:
            corrected.append((cluster, g))

# --- 4. Write corrected cluster table ---
with open(corrected_tsv_file, "w") as out:
    for cluster, gene in corrected:
        out.write(f"{cluster}\t{gene}\n")

# --- 5. Extract representative sequences ---
wanted = {gene for _, gene in corrected}
seqs = {}
with open(genes_file) as f:
    seq_id = None
    seq_buf = []
    for line in f:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if seq_id and seq_id in wanted:
                seqs[seq_id] = "".join(seq_buf)
            seq_id = line[1:].split()[0]
            seq_buf = []
        else:
            seq_buf.append(line)
    if seq_id and seq_id in wanted:
        seqs[seq_id] = "".join(seq_buf)

# --- 6. Write representative FASTA ---
with open(nr_genes_file, "w") as out:
    for gene in sorted(seqs.keys()):
        seq = seqs[gene]
        out.write(f">{gene}\n")
        for i in range(0, len(seq), 60):
            out.write(seq[i:i+60] + "\n")

