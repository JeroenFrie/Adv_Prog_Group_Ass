from sklearn.model_selection import KFold

def cross_validation_XY(X, Y, k=5, scoring='accuracy'):
    """
    Input:
    X = The input features
    Y = The target variable
    k = Number of folds for cross-validation
    scoring = The scoring metric for evaluation (default: 'accuracy')

    Output:
    train_test_sets_XY = List of tuples containing training and test sets for each fold
    """

    train_test_sets_XY = []

    for train, test in KFold(n_splits=k).split(X):
        Y_train, Y_test = Y.iloc[train], Y.iloc[test]
        X_train, X_test = X.iloc[train], X.iloc[test]

        train_test_sets_XY.append((X_train, X_test, Y_train, Y_test))

    return train_test_sets_XY

