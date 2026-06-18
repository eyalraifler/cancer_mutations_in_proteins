import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report

def run_logistic_regression(df, analysis_name="Combined"):
    print(f"\n" + "="*40)
    print(f"Running Logistic Regression for: {analysis_name}")
    print("="*40)

    X = df[['is_in_domain']]
    y = df['has_mutation']

    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)

    model = LogisticRegression(class_weight='balanced')
    model.fit(X_train, y_train)

    intercept = model.intercept_[0]
    coef = model.coef_[0][0]
    odds_ratio = np.exp(coef)

    # p-value via statsmodels Logit (MLE, no regularization)
    X_sm = sm.add_constant(X)
    sm_result = sm.Logit(y, X_sm).fit(disp=0)
    p_value = sm_result.pvalues['is_in_domain']

    print(f"Intercept: {intercept:.4f}")
    print(f"Coefficient (is_in_domain): {coef:.4f}")
    print(f"Odds Ratio: {odds_ratio:.4f}")
    print(f"P-value (is_in_domain): {p_value:.4e}")

    y_pred = model.predict(X_test)
    report = classification_report(y_test, y_pred)
    print(report)

    with open('results/logistic_regression_results.txt', 'a') as f:
        f.write(f"========================================\n")
        f.write(f"Analysis: {analysis_name}\n")
        f.write(f"========================================\n")
        f.write(f"Intercept: {intercept:.4f}\n")
        f.write(f"Coefficient (is_in_domain): {coef:.4f}\n")
        f.write(f"Odds Ratio: {odds_ratio:.4f}\n")
        f.write(f"P-value (is_in_domain): {p_value:.4e}\n\n")
        f.write("Classification Report:\n")
        f.write(report)
        f.write("\n\n")
        
