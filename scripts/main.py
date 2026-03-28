import re

EGFR_mutations_dict = {}
KRAS_mutations_dict = {}
'''
key: index of mutation
value: number of times appeared
'''


def get_list_of_all_protein_changes(mutations_file):
    ''''
    accepst: a mutation file
    return: a list of all protein changes in the mutation file 
    '''
    protein_changes = []
    for line in mutations_file:
        if line[0:4] == "Gene":
            continue
        line = line.strip()
        line = line.split("\t")
        change = line[5]
        protein_changes.append(change)
    return protein_changes


def index_of_mutation(protein_change):
    '''
    accepts: a line from the mutation file, which contains the protein change 
    return: the index of the mutation, which is the number in the middle of the protein change
    '''
    
    




if __name__ == "__main__":
    EGFR_protein = open("data/EGFR.fasta.txt", "r")
    Kras_protein = open("data/KRAS.fasta.txt", "r")
    EGFR_mutations = open("data/EGFR_mutations.tsv", "r")
    KRAS_mutations = open("data/KRAS_mutations.tsv", "r")
    
    EGFR_mutations_list = get_list_of_all_protein_changes(EGFR_mutations)
    print(EGFR_mutations_list)
