# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 00:22:17 2023

@author: 20201969
"""

"""
In this file the possible neural networks are listed we can use for 
            determining whether a inhibitor actually inhibits ALDH1
"""
#%% MLP

import pandas as pd
import numpy as np
from keras.models import Sequential
from keras.layers import Dense
from sklearn.model_selection import train_test_split
from sklearn.metrics import balanced_accuracy_score
#%%
# Load in the data frame
data = pd.read_csv('Descriptors_Vals_2D_3D.csv')
#%% Only run this part once for splitting the data and such
# Assuming you have variables (X) and corresponding binary labels (y) as your training data
X = data.iloc[:, 2:]    # All the columns except the molecule names and whether they inhibit or not
y = data.iloc[:, 1]  # The column which contains whether the molecule inhibits

# Determine the input dimension based on the number of variables
input_dim = X.shape[1]

#Splitting data in training and validation set
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

#%%Run this part for the actual model
# Create a Sequential model
model = Sequential()

# Add the input layer and the first hidden layer
model.add(Dense(64, activation='relu', input_shape=(input_dim,)))

# Add more hidden layers
model.add(Dense(64, activation='relu'))
model.add(Dense(64, activation='relu'))

# The output layer for a binary output,
#   if a continuous output is desired: change 'sigmoid' to 'relu'
model.add(Dense(1, activation='sigmoid'))

# Compile the model
model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

# Training the model
model.fit(X_train, y_train, epochs=10, batch_size=32, validation_split=0.2)

# Print the model summary
model.summary()

# Evaluate the model on the validation data
loss, accuracy = model.evaluate(X_test, y_test)
print("Validation Loss:", loss)
print("Validation Accuracy:", accuracy)

# Make predictions on the validation data
y_pred = model.predict(X_test)

# Convert probabilities to binary predictions
y_pred_binary = (y_pred > 0.5).astype(int)

# Calculate balanced accuracy
balanced_acc = balanced_accuracy_score(y_test, y_pred_binary)
print("Balanced Accuracy:", balanced_acc)


#TODO comparative method for true value and predicted value (RMES)
#TODO retrain on all the data
#TODO use sensitivity also

#%% Random Forest
from sklearn.ensemble import RandomForestClassifier

rf_model = RandomForestClassifier(n_estimators=400)
rf_model.fit(X_train, y_train)
rf_predictions = rf_model.predict(X_test)
balanced_acc = balanced_accuracy_score(y_test, rf_predictions)
print("Balanced Accuracy:", balanced_acc)

