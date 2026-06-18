import pandas as pd
from config_data import KRAS_DOMAINS_BY_DEFENITION, KRAS_ALL_REGIONS_WITH_FUNCTION, EGFR_DOMAINS_BY_DEFENITION, EGFR_ALL_REGIONS_WITH_FUNCTION, DATA_PATHS
from logistic_regression import run_logistic_regression
from processor import (
    get_protein_sequence_length, 
    get_list_of_all_protein_changes, 
    get_mutation_indices,
    build_regression_dataframe,
    save_domain_mutations_report
)

'''
For KRAS choose:
KRAS_DOMAINS_BY_DEFENITION OR KRAS_ALL_REGIONS_WITH_FUNCTION

For EGFR choose:
EGFR_DOMAINS_BY_DEFENITION OR EGFR_ALL_REGIONS_WITH_FUNCTION
'''

CHOSEN_KRAS_DOMAINS_DICTIONARY = KRAS_ALL_REGIONS_WITH_FUNCTION
CHOSEN_EGFR_DOMAINS_DICTIONARY = EGFR_ALL_REGIONS_WITH_FUNCTION


if __name__ == "__main__":
    # יצירת DataFrame עבור EGFR
    with open(DATA_PATHS["EGFR_FASTA"], "r") as f_seq, open(DATA_PATHS["EGFR_MUTATIONS"], "r") as f_mut:
        egfr_len = get_protein_sequence_length(f_seq)
        egfr_mut_list = get_list_of_all_protein_changes(f_mut)
        egfr_indices = get_mutation_indices(egfr_mut_list)
        
        df_egfr = build_regression_dataframe("EGFR", egfr_len, CHOSEN_EGFR_DOMAINS_DICTIONARY, egfr_indices)

    # יצירת DataFrame עבור KRAS
    with open(DATA_PATHS["KRAS_FASTA"], "r") as f_seq, open(DATA_PATHS["KRAS_MUTATIONS"], "r") as f_mut:
        kras_len = get_protein_sequence_length(f_seq)
        kras_mut_list = get_list_of_all_protein_changes(f_mut)
        kras_indices = get_mutation_indices(kras_mut_list)
        
        df_kras = build_regression_dataframe("KRAS", kras_len, CHOSEN_KRAS_DOMAINS_DICTIONARY, kras_indices)

    # איחד הDtaFrames של EGFR ו-KRAS לDataFrame אחד
    final_df = pd.concat([df_egfr, df_kras])
    #המרה לקובץ CSV
    final_df.to_csv("results/protein_mutation_data.csv", index=False)
    
    print("Success! Data table created.")
    
    # יוצרת סיכום עבור כל חלבון של כמות המוטציות בכל דומיין עם טווחי הדומיין
    #  וכמות העמדות וסכום המוטציות שלא נמצאים בדומיין כלשהוא 
    save_domain_mutations_report(df_egfr, df_kras, CHOSEN_EGFR_DOMAINS_DICTIONARY, CHOSEN_KRAS_DOMAINS_DICTIONARY)
    
    # איפוס קובץ התוצאות של הרגרסיה הלוגיסטית  
    open('results/logistic_regression_results.txt', 'w').close()
    
    # הרצת רגרסיה לוגיסטית על כל הנתונים המשולבים ובנוסף על כל חלבון בנפרד
    run_logistic_regression(final_df, analysis_name="Combined (EGFR & KRAS)")
    run_logistic_regression(df_egfr, analysis_name="EGFR Only")
    run_logistic_regression(df_kras, analysis_name="KRAS Only")
    
