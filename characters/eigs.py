from math import *
from numpy import matrix as m, reshape, identity as Id, linalg as la, conjugate as con, random as r

def innprod(cl, char1, char2):
    order = 0
    for i in cl:
        order += len(cl[i])
    innprod = 0
    for i in range(len(char1)):
        innprod += len(cl[i])*char1[i, 0]*con(char2[i, 0])
    innprod /= order
    return innprod

def l2d(l):
    D = {}
    for i in l:
        D[i] = i
    return D

def mat(sc, i):
    d = len(sc)
    M = Id(d)
    for j in range(d):
        for k in range(d):
            M[j, k] = sc[i][j][k]
    return M

def eigvecs(sc):
    d = len(sc)
    M = 0*Id(d)
    while len(l2d(la.eig(M)[0])) != d:
        M = 0*Id(d)
        weights = r.rand(d)
        for i in range(d):
            M += mat(sc, i)*weights[i]
            
    return la.eig(M)[1]

def char(vecs, cl, i):
    d = len(vecs)
    v = []
    for j in range(d):
        v.append(vecs[j, i]/len(cl[j]))
    vec = m(v).reshape(d, 1)
    vec /= vec[d-1, 0]
    return vec/sqrt(innprod(cl, vec, vec))

def compilechars(sc, cl):
    d = len(sc)
    ct = Id(d+1)*0
    vecs = eigvecs(sc)
    for i in range(d):
        ct[0, i+1] = len(cl[i])
        for j in range(d):
            v = char(vecs, cl, j)
            for k in range(d):
                ct[j+1, k+1] = v[k, 0]
    table = Id(d+1)*0
    l = [[ct[i+1, d], i+1] for i in range(d)]
    l.sort()
    for j in range(d):
        table[j+1] = ct[l[j][1]]
    for k in range(d):
        table[k+1, 0] = k+1
        table[0, k+1] = len(cl[k])
    return table
