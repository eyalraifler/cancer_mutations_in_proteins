import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report

def run_logistic_regression():
    # 1. טעינת הנתונים
    df = pd.read_csv('results/protein_mutation_data.csv')

    # 2. הגדרת המשתנים
    # אנחנו מתמקדים ב-is_in_domain כמשתנה המנבא העיקרי
    X = df[['is_in_domain']] 
    y = df['has_mutation']

    # 3. חלוקה לסט אימון ובדיקה עם שמירה על יחס המוטציות (stratify)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)

    # 4. בניית המודל עם משקולות מאוזנות
    # זה קריטי כי יש מעט מאוד מוטציות ביחס לאתרים רגילים
    model = LogisticRegression(class_weight='balanced')
    model.fit(X_train, y_train)

    # 5. הצגת התוצאות
    intercept = model.intercept_[0]
    coef = model.coef_[0][0]

    print(f"Intercept: {intercept:.4f}")
    print(f"Coefficient (is_in_domain): {coef:.4f}")

    # חישוב ה-Odds Ratio
    import numpy as np
    odds_ratio = np.exp(coef)
    print(f"Odds Ratio: {odds_ratio:.4f}")

    y_pred = model.predict(X_test)
    report = classification_report(y_test, y_pred)
    print(report)

    with open('results/logistic_regression.txt', 'w') as f:
        f.write(f"Intercept: {intercept:.4f}\n")
        f.write(f"Coefficient (is_in_domain): {coef:.4f}\n")
        f.write(f"Odds Ratio: {odds_ratio:.4f}\n\n")
        f.write("Classification Report:\n")
        f.write(report)