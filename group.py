from math import *
from numpy import matrix as m, reshape, identity as Id
from itertools import product

def hmat(m, p):
    h = 0
    d = len(m)
    for i in range(d):
        for j in range(d):
            h += m[i, j]*p**(i*d+j)
    return h

def hvec(v, p):
    h = 0
    d = len(v)
    for i in range(d):
        h = h + v[i, 0]*p**(i)
    return h

class Node(object):
    def __init__(self, data, label):
        self.data = data
        self.label = label
        self.children = {}
        self.parent = []
    def add_child(self, obj, p):
        self.children[hvec(obj.data, p)] = obj
    def add_parent(self, obj):
        self.parent.append(obj)

class Tree(object):
    def __init__(self, root):
        self.nodes = {-1: root}
    def add_node(self, obj, p):
        self.nodes[hvec(obj.data, p)] = obj

def vecin(v, S, p):
    h = hvec(v, p)
    for key in S.keys():
        if h == key:
            return 1
    return 0

def matin(m, S, p):
    h = hmat(m, p)
    for key in S.keys():
        if h == key:
            return 1
    return 0

   
def find(tree, node, p):
    h = hvec(node, p)
    for key in tree.nodes.keys():
        if hvec(tree.nodes[key].data, p) == h:
            return tree.nodes[key]

def schreier_tree(S, a, p):
    root = Node(a, None)
    tree = Tree(root)
    points = {hvec(a, p): a}
    gen = {hvec(a, p): a}
    while gen != {}:
        children = {}
        for key1 in gen.keys():
            for key2 in S.keys():
                q = S[key2]*gen[key1]%p
                child = Node(q, S[key2])
                if vecin(q, points, p) == 0:
                    points[hvec(q, p)] = q
                    child.add_parent(find(tree, gen[key1], p))
                    find(tree, gen[key1], p).add_child(child, p)
                    tree.add_node(child, p)
                    children[hvec(q, p)] = q
        gen = children
    return tree

def inv(M, p):
    N = M
    while hmat(M*N%p, p)!= hmat(Id(len(M)), p):
        N = M*N % p
    return N

def cosetrep(S, a, q, p):
    T = schreier_tree(S, a, p)
    g = Id(len(a))
    h = hvec(a, p)
    while hvec(q, p) != h:
        edge = find(T, q, p).label
        g = g*edge%p
        q = inv(edge, p)*q%p
    return g

def vec(n, d, p):
    if n == 0:
        return m([0]*d).reshape(d, 1)
    v = []
    while n:
        v.append(n%p)
        n /= p
    if len(v) < d:
        for i in range(d-len(v)):
            v.append(0)
    return m(v[::-1]).reshape(d, 1)

def newbp(M):
    d = len(M)
    for i in range(d):
        for j in range(d):
            if i != j:
                if M[i, j] != 0:
                    v = [0]*d
                    v[j] = 1
                    return m(v).reshape(d, 1)
    for i in range(d):
        for j in range(d):
            if i != j:
                if M[i, i] != M[j, j]:
                    u = [0]*d
                    w = [0]*d
                    u[i] = 1
                    w[j] = 1
                    v = [u[i]+w[i] for i in range(d)]
                    return m(v).reshape(d, 1)
    v = [0]*d
    v[0] = 1
    return m(v).reshape(d, 1)

def doesstab(m, B, p):
    if len(B) == 0:
        return 1
    for i in range(len(B)):
        if hvec(m*B[i]%p, p) != hvec(B[i], p):
            return 0
    return 1

def initBSGS(S, B, p):
    base = B
    sgs = {}
    for key in S.keys():
        if doesstab(S[key], base, p):
            point = newbp(S[key])
            base.append(point)
        sgs[hmat(S[key], p)] = S[key]
        sinv = inv(S[key], p)
        sgs[hmat(sinv, p)] = sinv
    return (sgs, base)

def genset(S, B, p, i):
    if i == 0:
        return S
    T = {}
    for key in S.keys():
        c = 0
        for j in range(i):
            if hvec(S[key]*B[j]%p, p) != hvec(B[j], p):
                c = 1
        if c == 0:
            T[hmat(S[key], p)] = S[key]
    return T

def reducegen(S, B, p, m):
    for i in range(m):
        T = genset(S, B, p, i)
        for key1 in T.keys():
            T1 = T[key1]
            for key2 in T.keys():
                T2 = T[key2]
                if hmat(T1, p) != hmat(T2, p):
                    if hvec(T1*B[i]%p, p) == hvec(T2*B[i]%p, p):
                        S[hmat(T1*inv(T2, p)%p, p)] = T1*inv(T2, p)%p
                        del S[key2]
    return S

def schreiergen(T, root, point, s, p):
    t1 = cosetrep(T, root, point, p)
    t2 = cosetrep(T, root, s*point%p, p)
    return inv(t2, p)*s*t1%p
    
def extendBSGS(B, S, i, p):
    d = len(B[0])
    Base = {}
    for vec in B:
        Base[hvec(vec, p)] = vec
    h = hmat(Id(d), p)
    T = genset(S, B, p, i)
    tree = schreier_tree(T, B[i], p)
    for key1 in tree.nodes.keys():
        for key2 in T.keys():
            gen = schreiergen(T, B[i], tree.nodes[key1].data, T[key2], p)
            if hmat(gen, p) != h and matin(gen, S, p) == 0:
                S[hmat(gen, p)] = gen
                if hmat(gen*gen, p) != h and matin(inv(gen, p), S, p) == 0:
                    S[hmat(inv(gen, p), p)] = inv(gen, p)
                if doesstab(gen, B, p):
                    newp = newbp(gen)
                    if vecin(newp, Base, p) == 0:
                        B.append(newp)
                        Base[hvec(newp, p)] = newp
    return(B, S)

def BSGS(S, p):
    sgs, base = initBSGS(S, [], p)
    l = len(base)
    for i in range(l):
        base, sgs = extendBSGS(base, sgs, i, p)
        sgs = reducegen(sgs, base, p, i)
    return (base, sgs)

def compilereps(S, a, p):
    reps = {}
    T = schreier_tree(S, a, p)
    g = Id(len(a))
    h = hvec(a, p)
    for key in T.nodes.keys():
        g = cosetrep(S, a, T.nodes[key].data, p)
        reps[hmat(g, p)] = g
    reps[hmat(Id(len(a)), p)] = Id(len(a))
    repslist = []
    for key in reps.keys():
        repslist.append(reps[key])
    return repslist

def trans(base, sgs, p):
    stabs = []
    trans = []
    for j in range(len(base)):
        T = genset(sgs, base, p, j)
        stabs.append(T)
    for k in range(len(stabs)):
        trans.append(compilereps(stabs[k], base[k], p))
    return trans

def order(trans):
    translen = []
    for i in trans:
        translen.append(len(i))
    order = 1
    for i in translen:
        order *= i
    return order

def elmnts(trans, p, d):
    fact = list(product(*trans))
    elmnts = []
    for i in fact:
        g = Id(d)
        for j in i:
            g = g*j%p
        elmnts.append(g)
    return elmnts

m1 = m('1 0 2; 0 0 1; 1 2 0')
m2 = m('0 1 2; 2 2 2; 0 0 1')
m3 = m('1 2 4; 2 4 2; 1 3 0')
m4 = m('2 1 2; 0 2 1; 2 0 3')
m5 = m('1 0 0; 0 1 5; 4 4 0')
m6 = m('5 4 2; 0 1 0; 1 1 2')
m7 = m('0 0 1 1; 1 1 1 0; 0 1 1 1; 1 1 0 0')
m8 = m('1 0 1 0; 1 1 1 0; 1 1 1 1; 0 1 1 1')
m9 = m('1 1; 0 1')
m10 = m('1 0; 1 1')

set1 = {hmat(m1, 3):m1, hmat(m2, 3):m2}
set2 = {hmat(m3, 5):m3, hmat(m4, 5):m4}
set3 = {hmat(m5, 7):m5, hmat(m6, 7):m6}
set4 = {hmat(m7, 2):m7, hmat(m8, 2):m8}
set5 = {hmat(m9, 2):m9, hmat(m10, 2):m10}