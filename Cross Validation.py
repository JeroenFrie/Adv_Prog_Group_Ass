from sklearn.model_selection import KFold

def cross_validation_Y(Y, k=5, scoring = 'accuracy'):
    """
    Input:
    Y = The y-coordinate values, the values of a descriptor (for instance molecular weight)
    k = Number of fold of the crossvalidation

    Output:
    train_test_sets = This is a list of tuples, each tuple contains training and test sets for each k-fold
    """

    train_test_sets_Y = []

    for train, test in KFold(n_splits=k).split(Y):
        Y_train, Y_test = Y[train], Y[test]


        train_test_sets_Y.append((None, None, Y_train, Y_test))

    return train_test_sets_Y

def cross_validation_XY(X, Y, k=5, scoring = 'accuracy'):
    """
    Input:
    Y = The y-coordinate values, the values of a descriptor (for instance molecular weight)
    k = Number of fold of the crossvalidation

    Output:
    train_test_sets = This is a list of tuples, each tuple contains training and test sets for each k-fold
    """

    train_test_sets_XY = []

    for train, test in KFold(n_splits=k).split(X):
        Y_train, Y_test = Y[train], Y[test]
        X_train, X_test = X[train], X[test]

        train_test_sets_XY.append((X_train, X_test, Y_train, Y_test))

    return train_test_sets_XY

