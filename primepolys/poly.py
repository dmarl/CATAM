from math import *
from numpy import random as r

def mult(a, b, p):
    ab = [0 for i in range(len(a)+len(b)-1)]
    for i in range(len(a)):
        for j in range(len(b)):
            ab[i+j] += a[i]*b[j]
            ab[i+j] %= p
    return ab

def diff(a, b, p):
    l = max(len(a), len(b))
    A = [0 for i in range(l)]
    B = [0 for i in range(l)]
    for i in range(len(a)):
        A[i] = a[i]
    for i in range(len(b)):
        B[i] = b[i]
    C = [(A[i]-B[i])%p for i in range(l)]
    return C

def add(a, b, p):
    l = max(len(a), len(b))
    A = [0 for i in range(l)]
    B = [0 for i in range(l)]
    for i in range(len(a)):
        A[i] = a[i]
    for i in range(len(b)):
        B[i] = b[i]
    C = [(A[i]+B[i])%p for i in range(l)]
    return C

def modexp(a, b, p):
    #computes a**b mod p with binary modular exponentiation
    modexp = 1
    while b != 0:
        if b%2 == 1:
            modexp = modexp*a%p
        a = a*a%p
        b /= 2
    return modexp

def inv(a, p):
    #inverse of a mod p
    return modexp (a, p-2, p)

def skim(a, p):
    while a[0]%p == 0 and len(a) > 1:
        del a[0]
    return a
    
def divpoly(a, b, p):
    #returns the quotient q and remainder
    #r of polynomial division of a by b modulo p
    A = a
    a = skim(a, p)
    b =  skim(b, p)
    if len(a) < len(b):
        return 0, a
    q = [0 for i in range(len(a)-len(b)+1)]
    for i in range(len(q)):
        k = a[0]*inv(b[0], p)%p
        q[i] = k
        a = diff(a, mult(b, [k]+[0 for i in range(len(q)-i-1)], p), p)
        del a[0]
    q = skim(q, p)
    return q, skim(diff(A, mult(b, q, p), p), p)

def polygcd(a, b, p):
    c = b
    if len(a) < len(b):
        b = a
        a = c
    while divpoly(a, b, p)[1] != [0]:
        c = b
        b = divpoly(a, b, p)[1]
        a = c
    if len(b) == 1:
        return [1]
    return b

def pol(p):
    return [1]+[0 for i in range(p-2)]+[-1]+[0]

def polyexp(q, n, p):
    exp = [1]
    while n != 0:
        if n%2 == 1:
            exp = mult(exp, q, p)
        q = mult(q, q, p)
        n /= 2
    return exp

def squroot(a, p):
    f = [1, 0, -a%p]
    h = f
    c = 0
    while len(h) == len(f) or len(h) == 1:
        b = r.randint(0, p)
        g = diff(polyexp([1, b], (p-1)/2, p), [0 for i in range((p-1)/2)]+[1], p)
        h = polygcd(f, g, p)
        c += 1
    H = divpoly(f, h, p)[0]
    return h[1]*inv(h[0], p)%p, H[1]*inv(H[0], p)%p, c

def roots(f, p):
    f = polygcd([1]+[0 for i in range(p-2)]+[-1]+[0], f, p)
    if len(f) == 1:
        print "poly has no solution mod", p
        return
    h = f
    if len(f) == 2:
        return -(f[1]*inv(f[0], p))%p
    while len(h) != 2:
        b = r.randint(1, p)
        g = diff(polyexp([1, b], (p-1)/2, p), [0 for i in range((p-1)/2)]+[1], p)
        h = polygcd(f, g, p)
    H = divpoly(f, h, p)[0]
    return roots(h, p), roots(H, p)

def compileroots(l):
    roots = []
    for i in l:
        roots.append(i)
    return roots
        