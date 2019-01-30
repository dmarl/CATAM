from math import *
from numpy import matrix as m

def mateq(A, B):
    if len(A) != len(B):
        return 0
    for i in range(len(A)):
        for j in range(len(B)):
            if A[i, j] != B[i, j]:
                return 0
    return 1

def ord(A, p):
    id = A* (A.I)
    B = A
    c = 1
    while mateq(B, id) == 0:
        B = B*A%p
        c = c + 1
    return c