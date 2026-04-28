import re
import pandas as pd




def get_protein_sequence(protein_file):
    """
    קלט: קובץ חלבון פתוח (פורמט FASTA).
    פלט: רצף החלבון כמחרוזת.
    """
    sequence = ""
    for line in protein_file:
        if line.startswith('>'):
            continue
        sequence += line.strip()
    return sequence

def get_protein_sequence_length(protein_file):
    """
    קלט: קובץ חלבון פתוח (פורמט FASTA).
    פלט: אורך רצף החלבון (מספר שלם).
    """
    sequence = get_protein_sequence(protein_file)
    return len(sequence)


def get_list_of_all_protein_changes(mutations_file):
    """
    קלט: קובץ מוטציות פתוח (פורמט TSV).
    פלט: רשימה של מחרוזות המייצגות את שינויי החלבון (למשל: 'L858R').
    """
    protein_changes = []
    for line in mutations_file:
        if line[0:4] == "Gene":
            continue
        line = line.strip()
        line = line.split("\t")
        if len(line) > 5:  # לוודא שהעמודה קיימת
            change = line[5]
            protein_changes.append(change)
    return protein_changes


def get_mutation_indices(mutation_list):
    """
    קלט: רשימת מחרוזות של מוטציות.
    פלט: רשימה של רשימות, כאשר כל תת-רשימה מכילה את האינדקסים המספריים שנמצאו במוטציה.
    """
    all_indices = []
    for mut in mutation_list:
        numbers = re.findall(r'\d+', mut)
        # הופך את כל המספרים שמצאנו במוטציה ל-integers
        indices = [int(n) for n in numbers]
        if indices:
            all_indices.append(indices)
    return all_indices


def update_domain_counts(indices_list, domain_dict):
    """
    קלט: רשימת אינדקסים (פלט של get_mutation_indices) ומילון דומיינים לעדכון.
    פלט: אין (מעדכן את המילון הקיים In-place).
    """
    for indices in indices_list:
        # עוברים על כל טווח שקיים במילון (למשל (713, 979))
        for (start, end) in domain_dict.keys():
            # בודקים אם לפחות אחד מהאינדקסים של המוטציה נמצא בתוך הטווח
            if any(start <= idx <= end for idx in indices):
                domain_dict[(start, end)] += 1
