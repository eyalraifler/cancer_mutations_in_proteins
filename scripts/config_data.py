KRAS_DOMAINS = {
    (10, 14): "P-loop",
    (30, 40): "Switch I",
    (58, 72): "Switch II",
    (185, 188): "CAAX motif"
}

EGFR_DOMAINS = {
    (1, 24): "Signal Peptide",
    (25, 189): "Domain I (L1)",
    (190, 333): "Domain II (CR1)",
    (334, 504): "Domain III (L2)",
    (505, 645): "Domain IV (CR2)",
    (646, 668): "Transmembrane (TM)",
    (669, 712): "Juxtamembrane (JM)",
    (713, 973): "Tyrosine Kinase Domain",
    (741, 750): "P-loop (Nucleotide Binding)",
    (850, 875): "Activation Loop (A-loop)",
    (974, 1210): "C-terminal Tail"
}

# files paths
DATA_PATHS = {
    "EGFR_FASTA": "data/EGFR.fasta.txt",
    "KRAS_FASTA": "data/KRAS.fasta.txt",
    "EGFR_MUTATIONS": "data/EGFR_mutations.tsv",
    "KRAS_MUTATIONS": "data/KRAS_mutations.tsv"
}