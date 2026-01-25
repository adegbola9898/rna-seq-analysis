
import pandas as pd
import warnings
warnings.filterwarnings('ignore')
warnings.filterwarnings('ignore', category=DeprecationWarning)

filtered_count_data = pd.read_csv('/content/filtered_count_data.csv')
filtered_count_data

col_data = pd.read_csv('/content/col_data.csv')
col_data

# Transpose the filtered_count_data dataframe
filtered_count_data_transposed = filtered_count_data.T
filtered_count_data_transposed.reset_index(inplace=True)
filtered_count_data_transposed.columns = filtered_count_data_transposed.iloc[0] # Set the first row as column header
filtered_count_data_transposed = filtered_count_data_transposed.drop(filtered_count_data_transposed.index[0]) # Drop the first row


conditions = col_data[['Conditions']] # Extract 'Conditions' column to append it later
merged_df = pd.merge(col_data.drop(columns=['Conditions']), filtered_count_data_transposed, left_on='Sample', right_on='gene_id', how='right')
merged_df_final = pd.concat([merged_df, conditions], axis=1) # Append 'Conditions' column at the end

merged_df_final

data = merged_df_final.iloc[:, 1:] #drops the Sample column of the merged dataframe

data.describe()

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

# Exclude non-numerical columns
numerical_data = data.drop(columns=['gene_id', 'Conditions'])

# Standardize the data
scaler = StandardScaler()
scaled_data = scaler.fit_transform(numerical_data)

# Perform PCA
pca = PCA(n_components=0.95)  # Keep 95% of variance
principal_components = pca.fit_transform(scaled_data)

# Number of components chosen
num_components = pca.n_components_

# Variance explained by each of the selected components
explained_variance = pca.explained_variance_ratio_

num_components, explained_variance

# Scatter Plot of the First Two Principal Components
plt.figure(figsize=(10, 5))

plt.subplot(1, 2, 1)
plt.scatter(principal_components[:, 0], principal_components[:, 1])
plt.title('PCA - First two components')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')

# Scree Plot
plt.subplot(1, 2, 2)
plt.plot(range(1, num_components + 1), explained_variance, marker='o', linestyle='--')
plt.title('Scree Plot')
plt.xlabel('Number of Components')
plt.ylabel('Variance Explained')

plt.tight_layout()
plt.show()

from sklearn.svm import SVC
from sklearn.model_selection import cross_val_score

# Prepare the features (PCA components) and labels
X = principal_components  # PCA-reduced data
y = data['Conditions'].values

# Initialize the SVM model with a linear kernel
svm_model = SVC(kernel='linear', C=1)

# Train the SVM model on the entire dataset
svm_model.fit(X, y)

# Perform cross-validation to estimate the model's performance
cross_val_scores = cross_val_score(svm_model, X, y, cv=5)

# Output the mean cross-validation score
print("Mean cross-validation accuracy:", cross_val_scores.mean())

from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score

# Prepare the features and labels
X = principal_components  # PCA-reduced data
y = data['Conditions'].values

# Initialize classifiers
classifiers = {
    "Random Forest": RandomForestClassifier(n_estimators=100),
    "Gradient Boosting": GradientBoostingClassifier(n_estimators=100),
    "K-Nearest Neighbors": KNeighborsClassifier(n_neighbors=3),
    "Logistic Regression": LogisticRegression(max_iter=1000)
}

# Perform cross-validation and print out the mean accuracy for each classifier
for name, clf in classifiers.items():
    scores = cross_val_score(clf, X, y, cv=5)
    print(f"{name}: Mean cross-validation accuracy = {scores.mean():.2f}")



from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score

# Prepare the features and labels
X = principal_components  # PCA-reduced data
y = data['Conditions'].values

# Initialize classifiers
classifiers = {
    "Random Forest": RandomForestClassifier(n_estimators=100),
    "Gradient Boosting": GradientBoostingClassifier(n_estimators=100),
    "K-Nearest Neighbors": KNeighborsClassifier(n_neighbors=3),
    "Logistic Regression": LogisticRegression(max_iter=1000)
}

# Perform cross-validation and print out the mean accuracy for each classifier
for name, clf in classifiers.items():
    scores = cross_val_score(clf, X, y, cv=5)
    print(f"{name}: Mean cross-validation accuracy = {scores.mean():.2f}")
