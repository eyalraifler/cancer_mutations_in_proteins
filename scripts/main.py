import re

# מילון עבור חלבון KRAS
# המפתח: טווח חומצות אמינו (התחלה, סוף)
# הערך: מספר המוטציות שנמצאו בטווח זה
KRAS_domains = {
    (1, 165): 0,
    (10, 17): 0,
    (30, 76): 0,
    (166, 185): 0,
    (185, 188): 0
}

# מילון עבור חלבון EGFR
EGFR_domains = {
    (1, 165): 0,
    (166, 309): 0,
    (310, 480): 0,
    (481, 644): 0,
    (645, 668): 0,
    (669, 712): 0,
    (713, 979): 0,
    (980, 1210): 0
}


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
        if len(line) > 5:  # בדיקת בטיחות שהעמודה קיימת
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


if __name__ == "__main__":
    # פתיחת קבצי הנתונים
    EGFR_protein = open("data/EGFR.fasta.txt", "r")
    Kras_protein = open("data/KRAS.fasta.txt", "r")
    EGFR_mutations = open("data/EGFR_mutations.tsv", "r")
    KRAS_mutations = open("data/KRAS_mutations.tsv", "r")
    
    # 1. חילוץ רשימת המוטציות
    EGFR_mutations_list = get_list_of_all_protein_changes(EGFR_mutations)
    KRAS_mutations_list = get_list_of_all_protein_changes(KRAS_mutations)

    # 2. הפיכת המוטציות למספרים
    egfr_indices = get_mutation_indices(EGFR_mutations_list)
    kras_indices = get_mutation_indices(KRAS_mutations_list)

    # 3. עדכון המילונים שהגדרת למעלה
    update_domain_counts(egfr_indices, EGFR_domains)
    update_domain_counts(kras_indices, KRAS_domains)

    # הדפסת התוצאות לראות שזה עבד
    print("EGFR Domain Counts:", EGFR_domains)
    print("KRAS Domain Counts:", KRAS_domains)

    # סגירת קבצים
    EGFR_mutations.close()
    KRAS_mutations.close()
    EGFR_protein.close()
    Kras_protein.close()