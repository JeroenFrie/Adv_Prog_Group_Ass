# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 20:27:09 2023

@author: 20201969
"""
# Import shizzle

import pandas as pd
from sklearn.ensemble import VotingClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
import sklearn.preprocessing as sp
from CSV_Load import CSV_Loader

#%% Preparing the data
# load dataframe for training
known_molecules = CSV_Loader("Descriptors_Vals_2D_3D.csv")

# Split its ass
X = known_molecules.iloc[:, 2:]
y = known_molecules.iloc[:, 1]

# Scaling
scaler_type = sp.StandardScaler()
scaler_type.fit(X)

# Determine the input dimension based on the number of variables
input_dim = X.shape[1]

# Loading the dataset that we want to know the inhibition of
unknown = CSV_Loader("Unknown_Mol_Desc.csv")

# Splitting of the name column to get the descriptor values only
unknown_desc = unknown.iloc[:, 1:]

#%% Train the whole data set with Logistic regression
# Create individual classifiers
classifier1 = RandomForestClassifier(n_estimators=400)
classifier2 = LogisticRegression(max_iter=2000)

# Create the ensemble classifier using majority voting with probability estimation
ensemble_classifier = VotingClassifier(
    estimators=[('rfc', classifier1), ('sq', classifier2)],
    voting='soft',  # Use 'soft' for weighted voting with predicted probabilities
    flatten_transform=True,  # Enable probability estimation
)

# Train the ensemble classifier
ensemble_classifier.fit(X, y)

#%% Predicting the values for the unknown data
# Make probability predictions
y_pred_prob = ensemble_classifier.predict_proba(unknown_desc)

# Get the continuous output probabilities
continuous_output = y_pred_prob[:, 1] 

output_linked = pd.DataFrame({'Molecule': unknown.iloc[:, 0], 'Continuous Output': continuous_output})
# Sort the DataFrame based on the 'Continuous Output' column in descending order
output = output_linked.sort_values(by='Continuous Output', ascending=False)

# Select the top 100 rows
top_100_df = output.head(100)

# Save the top 100 rows to an Excel file
top_100_df.to_excel('top_100_molecules.xlsx', index=False)