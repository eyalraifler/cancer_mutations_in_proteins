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


def build_regression_dataframe(protein_name, seq_len, domains_dict, mutation_indices):
    """
    יוצר טבלה שבה כל שורה היא עמדה בחלבון (1 עד אורך החלבון).
    """
    # הופך רשימת רשימות [[12], [858]] לרשימה אחת שטוחה [12, 858]
    flat_mutations = [idx for sublist in mutation_indices for idx in sublist]
    
    rows = []
    for pos in range(1, seq_len + 1):
        # האם העמדה בתוך דומיין כלשהו? (0 או 1)
        is_in_domain = 0
        for (start, end) in domains_dict.keys():
            if start <= pos <= end:
                is_in_domain = 1
                break
        
        # האם יש מוטציה בעמדה הזו? (0 או 1)
        has_mutation = 1 if pos in flat_mutations else 0
        
        # כמה מוטציות יש בעמדה הזו?
        mut_count = flat_mutations.count(pos)
        
        rows.append({
            'protein': protein_name,
            'pos': pos,
            'is_in_domain': is_in_domain,
            'has_mutation': has_mutation,
            'mutation_count': mut_count
        })
    
    return pd.DataFrame(rows)


def count_positions_outside_domains(df, protein_name):
    """
    Count how many amino acid positions are outside all defined domains.
    """
    outside = df[(df['protein'] == protein_name) & (df['is_in_domain'] == 0)]
    return len(outside)


def count_mutations_outside_domains(df, protein_name):
    """
    Count how many mutations occur outside all defined domains.
    """
    outside = df[(df['protein'] == protein_name) & (df['is_in_domain'] == 0)]
    return int(outside['mutation_count'].sum())


def create_domain_mutations_summary(df, protein_name, domains_dict):
    """
    Create a summary of mutation counts per domain for a given protein.
    Returns a list of dictionaries with domain and mutation count information.
    """
    summary = []
    
    for (start, end), domain_name in domains_dict.items():
        # Filter dataframe for this domain
        domain_data = df[(df['pos'] >= start) & (df['pos'] <= end) & (df['protein'] == protein_name)]
        
        # Count mutations in this domain
        mutation_count = domain_data['mutation_count'].sum()
        
        summary.append({
            'protein': protein_name,
            'domain': domain_name,
            'range': f"{start}-{end}",
            'mutations': int(mutation_count)
        })
    
    return summary


def save_domain_mutations_report(df_egfr, df_kras, egfr_domains, kras_domains):
    """
    Save a report of mutations per domain for each protein.
    """
    filename = "results/domain_mutations.txt"
    
    egfr_summary = create_domain_mutations_summary(df_egfr, "EGFR", egfr_domains)
    kras_summary = create_domain_mutations_summary(df_kras, "KRAS", kras_domains)
    
    with open(filename, 'w') as f:
        f.write("MUTATION COUNT BY DOMAIN\n")
        f.write("=" * 60 + "\n\n")
        
        # EGFR summary
        f.write("EGFR PROTEIN\n")
        f.write("-" * 60 + "\n")
        f.write(f"{'Domain Name':<35} {'Range':<15} {'Mutations':<10}\n")
        f.write("-" * 60 + "\n")
        for entry in egfr_summary:
            f.write(f"{entry['domain']:<35} {entry['range']:<15} {entry['mutations']:<10}\n")

        egfr_outside_aa = count_positions_outside_domains(df_egfr, "EGFR")
        egfr_outside_mutations = count_mutations_outside_domains(df_egfr, "EGFR")
        f.write("\n")
        f.write(f"amount_of_amino_acids_not_in_domain = {egfr_outside_aa}\n")
        f.write(f"mutations_not_in_domain = {egfr_outside_mutations}\n")
        
        f.write("\n" + "=" * 60 + "\n\n")
        
        # KRAS summary
        f.write("KRAS PROTEIN\n")
        f.write("-" * 60 + "\n")
        f.write(f"{'Domain Name':<35} {'Range':<15} {'Mutations':<10}\n")
        f.write("-" * 60 + "\n")
        for entry in kras_summary:
            f.write(f"{entry['domain']:<35} {entry['range']:<15} {entry['mutations']:<10}\n")

        kras_outside_aa = count_positions_outside_domains(df_kras, "KRAS")
        kras_outside_mutations = count_mutations_outside_domains(df_kras, "KRAS")
        f.write("\n")
        f.write(f"amount_of_amino_acids_not_in_domain = {kras_outside_aa}\n")
        f.write(f"mutations_not_in_domain = {kras_outside_mutations}\n")
        
        f.write("\n" + "=" * 60 + "\n")
    
    print(f"Domain mutations report saved to: {filename}")