def has_same_known_substructure(A, B):
    if A.shape != B.shape:
        return False

    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            if A[i, j] == 0 and B[i, j] != 0:
                return False
            if B[i, j] == 0 and A[i, j] != 0:
                return False
            if A[i, j] == 1 and B[i, j] != 1:
                return False
            if B[i, j] == 1 and A[i, j] != 1:
                return False

    return True
