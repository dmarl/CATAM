from math import *
from numpy import matrix as m, reshape, identity as Id
from tree import schreier_tree as st, genset as gs

def hvec(v, p):
    h = 0
    d = len(v)
    for i in range(d):
        h = h + v[i, 0]*p**(i)
    return h

def hmat(m, p):
    h = 0
    d = len(m)
    for i in range(d):
        for j in range(d):
            h += m[i, j]*p**(i*d+j)
    return h

def compilereps(S, a, p):
    reps = {}
    T = st(S, a, p)
    g = Id(len(a))
    h = hvec(a, p)
    for key in T.nodes:
        q = T.nodes[key].data
        while hvec(q, p) != h:
            edge = find(T, q, p).label
            g = g*edge%p
            q = inv(edge, p)*q%p
        reps[hmat(g, p)] = g
    return reps


def trans(base, sgs, p):
    stabs = []
    trans = []
    for j in range(len(base)):
        T = gs(sgs, base, p, j)
        stabs.append(T)
    for k in range(len(stabs)):
        stabset = stabs[k]
        trans.append(compilereps(stabset, base[k], p))
    return trans
