import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report

def run_logistic_regression():
    df = pd.read_csv('results/protein_mutation_data.csv')

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

    with open('results/logistic_regression.txt', 'w') as f:
        f.write(f"Intercept: {intercept:.4f}\n")
        f.write(f"Coefficient (is_in_domain): {coef:.4f}\n")
        f.write(f"Odds Ratio: {odds_ratio:.4f}\n")
        f.write(f"P-value (is_in_domain): {p_value:.4e}\n\n")
        f.write("Classification Report:\n")
        f.write(report)
        
        
'''
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report, confusion_matrix

# 1. טעינת הנתונים מהקובץ שלך
# (הנחת עבודה: הקובץ נקרא 'protein_data.csv')
df = pd.read_csv("protein_data.csv")

# 2. הנדסת פיצ'רים: יצירת משתנה רציף הלוקח בחשבון את כמות המוטציות
# אנחנו משתמשים ב-log1p כדי להתמודד עם ערכי 0 ועם ה"זנב הארוך" של מאגר COSMIC
df["log_mutation_count"] = np.log1p(df["mutation_count"])

# 3. הגדרת המשתנים למודל
# המשתנים הבלתי תלויים (X) - כרגע מתבססים על כמות המוטציות (וניתן להוסיף פיצ'רים מבניים בעתיד)
X = df[["log_mutation_count"]]

# המשתנה התלוי (y) - המטרה הבינארית שלנו: האם המיקום בדומיין או לא
y = df["is_in_domain"]

# 4. חלוקה לסט אימון וסט בדיקה (כדי שנוכל להעריך את ביצועי המודל)
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y
)

# 5. נרמול הנתונים (סקיילינג) - שלב קריטי ברגרסיה לוגיסטית כאשר עובדים עם משתנים רציפים
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# 6. הגדרת המודל ואימונו
model = LogisticRegression()
model.fit(X_train_scaled, y_train)

# 7. חיזוי והערכת ביצועים על סט הבדיקה
y_pred = model.predict(X_test_scaled)

print("--- דוח ביצועי המודל (Classification Report) ---")
print(classification_report(y_test, y_pred))

# 8. הדפסת המקדם (Coefficient) של המודל כדי להבין את הקשר הביולוגי
intercept = model.intercept_[0]
coef_mutation = model.coef_[0][0]

print("\n--- משוואת המודל והמקדמים ---")
print(f"הקבוע (Intercept): {intercept:.4f}")
print(f"המקדם של log_mutation_count: {coef_mutation:.4f}")

if coef_mutation > 0:
    print(
        "פירוש ביולוגי: ככל שמיקום צובר יותר מוטציות ב-COSMIC, כך עולה הסיכוי שהוא בתוך הדומיין הפונקציונלי."
    )
else:
    print(
        "פירוש ביולוגי: לא נמצא קשר חיובי ישיר בין כמות המוטציות להימצאות בדומיין במודל הנוכחי."
    )
'''
