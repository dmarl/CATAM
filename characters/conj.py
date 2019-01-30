from math import *
from numpy import matrix as m, reshape, identity as Id
from group import hmat

def tr(m, p):
    tr = 0
    for i in range(len(m)):
        tr = tr + m[i, i]
    return tr%p

def isconj(g, h, group, p):
    if tr(g, p) != tr(h, p):
        return 0
    for k in group:
        if hmat(k*g%p, p) == hmat(h*k%p, p):
            return 1
    return 0     
        
def conj(group, p):
    classes = {} 
    for g in group:
        c = 0
        for key in classes.keys():
            if isconj(g, classes[key][0], group, p):
                classes[key].append(g)
                c = 1
        if c == 0:
            classes[len(classes)+1] = [g]
    return classes

def ordercl(cl):
    D = dict(cl)
    classes = []
    while len(D) > 0:
        c = [0, 0]
        for i in D:
            if len(D[i]) >= c[0]:
                c = [len(D[i]), i]
        classes.append(D[c[1]])
        del D[c[1]]
    classdict = {}
    for i in range(len(classes)):
        classdict[i] = classes[i]
    return classdict

        