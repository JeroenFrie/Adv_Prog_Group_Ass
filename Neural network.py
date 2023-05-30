"""
In this file the possible neural networks are listed we can use for 
            determining whether a inhibitor actually inhibits ALDH1
"""
#%% MLP
import numpy as np
from keras.models import Sequential
from keras.layers import Dense

import numpy as np
from keras.models import Sequential
from keras.layers import Dense

# Assuming you have variables (X) and corresponding binary labels (y) as your training data
X = np.array([[var1, var2, var3, ...], [var1, var2, var3, ...], ...])  # Replace with your variables
y = np.array([label1, label2, label3, ...])  # Replace with your binary labels

# Determine the input dimension based on the number of variables
input_dim = X.shape[1]

# Create a Sequential model
model = Sequential()

# Add the input layer and the first hidden layer
model.add(Dense(64, activation='relu', input_shape=(input_dim,)))

# Add more hidden layers
model.add(Dense(32, activation='relu'))
model.add(Dense(16, activation='relu'))

# Add the output layer
model.add(Dense(1, activation='sigmoid'))

# Compile the model
model.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy'])

# Print the model summary
model.summary()

# Train the model using the training data
model.fit(X, y, epochs=10, batch_size=32, validation_split=0.2)


#%%

