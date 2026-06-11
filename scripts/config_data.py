#for running the regression with different domain definitions, 
# we can switch between these dictionaries in main.py

#זווית 1
KRAS_DOMAINS_BY_DEFENITION = {
    (1, 166): "G Domain",
}

#זווית 2
KRAS_ALL_REGIONS_WITH_FUNCTION = {
    (10, 14): "P-loop",
    (30, 40): "Switch I",
    (58, 72): "Switch II",
    (185, 188): "CAAX motif"
}

#זווית 1
EGFR_DOMAINS_BY_DEFENITION = {
    (25, 189):  "Domain I (Ligand binding)",
    (190, 333): "Domain II (Dimerization interface)",
    (334, 504): "Domain III (Ligand binding)",
    (505, 645): "Domain IV (CR2)",
    (713, 973): "Tyrosine Kinase Domain",
}


#זווית 2
EGFR_ALL_REGIONS_WITH_FUNCTION = {
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


#רק האזורים המשמעותיים בכל דומיין
EGFR_DOMAINS_INPORTANT_ONES = {
    (572, 574): "Extracellular Domain",
    (279, 280): "Extracellular Domain",
    (712, 724): "Kinase Domain",
    (746, 753): "Kinase Domain",
    (858, 858): "Kinase Domain"
}

# files paths
DATA_PATHS = {
    "EGFR_FASTA": "data/EGFR.fasta.txt",
    "KRAS_FASTA": "data/KRAS.fasta.txt",
    "EGFR_MUTATIONS": "data/EGFR_mutations.tsv",
    "KRAS_MUTATIONS": "data/KRAS_mutations.tsv"
}