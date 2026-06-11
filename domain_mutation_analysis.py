import pandas as pd

df = pd.read_csv("results/protein_mutation_data.csv")

for protein, group in df.groupby("protein"):
    in_domain = group[group["is_in_domain"] == 1]
    not_in_domain = group[group["is_in_domain"] == 0]

    in_domain_count = len(in_domain)
    in_domain_mutated = in_domain["has_mutation"].sum()

    not_in_domain_count = len(not_in_domain)
    not_in_domain_mutated = not_in_domain["has_mutation"].sum()

    print(f"=== {protein} ===")
    print(f"Amino acids IN domain:       {in_domain_count}")
    print(f"  - With mutation:           {in_domain_mutated}")
    print(f"  - Without mutation:        {in_domain_count - in_domain_mutated}")
    print()
    print(f"Amino acids NOT in domain:   {not_in_domain_count}")
    print(f"  - With mutation:           {not_in_domain_mutated}")
    print(f"  - Without mutation:        {not_in_domain_count - not_in_domain_mutated}")
    print()
