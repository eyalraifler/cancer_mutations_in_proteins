KRAS_DOMAINS = {
    (1, 166): "G-domain",
    (10, 14): "P-loop",
    (30, 40): "Switch I",
    (58, 72): "Switch II",
    (167, 188): "Hypervariable(HVR)",
    (185, 188): "CAAX motif"
}

EGFR_DOMAINS = {
    (1, 165): "Domain I",
    (713, 979): "Kinase Domain",
    # ... שאר הדומיינים
}

# files paths
DATA_PATHS = {
    "EGFR_FASTA": "data/EGFR.fasta.txt",
    "KRAS_FASTA": "data/KRAS.fasta.txt",
    "EGFR_MUTATIONS": "data/EGFR_mutations.tsv",
    "KRAS_MUTATIONS": "data/KRAS_mutations.tsv"
}