from math import *
from numpy import matrix as m, reshape, identity as Id

def hmat(m, p):
    h = 0
    d = len(m)
    for i in range(d):
        for j in range(d):
            h += m[i, j]*p**(i*d+j)
    return h

def baseprod(group, bi, bj, p):
    prod = {}
    for g in group:
        prod[hmat(g, p)] = [0, g]
    for i in bi:
        for j in bj:
            prod[hmat(i*j%p, p)][0] += 1
    return prod

def sc(group, classdict, p):
    #returns the array of structure constants \nu_ijk of G
    r = len(classdict)
    classes = []
    sc = [[[0 for i in range(r)] for j in range(r)] for k in range(r)]
    for key in classdict.keys():
        l = []
        for elmnt in classdict[key]:
            l.append(elmnt)
        classes.append(l)
    for i in range(r):
        for j in range(r):
            for k in range(r):
                bibj = baseprod(group, classes[i], classes[j], p)
                sc[i][j][k] = bibj[hmat(classes[k][0], p)][0]
    return sc
