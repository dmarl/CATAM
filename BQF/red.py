from math import *
from itertools import product
from numpy import random as r
 
def gcd(a, b):
    a = abs(a)
    b = abs(b)
    if b > a:
        c = b
        b = a
        a = c
    if b == 0:
        return a
    while a%b != 0:
        c = b
        b = a%b
        a = c
    return b

def GCD(a, b, c):
    a = abs(a)
    b = abs(b)
    c = abs(c)
    d = gcd(a, b)
    return gcd(c, d)

def bezout(a, b):
    flip = 0
    if a < b:
        flip = 1
        c = b
        b = a
        a = c
    if b == 0:
        return [a, 1, 0]
    s = [1, 0]
    t = [0, 1]
    while a%b != 0:
        q = a/b
        c = b
        b = a%b
        a = c
        S = s[0] - q*s[1]
        T = t[0] - q*t[1]
        s = [s[1], S]
        t = [t[1], T]
    if flip == 0:
        return [b, s[1], t[1]]
    else:
        return [b, t[1], s[1]]

def modexp(a, b, p):
    #computes a**b mod p with binary modular exponentiation
    modexp = 1
    while b != 0:
        if b%2 == 1:
            modexp = modexp*a%p
        a = a*a%p
        b /= 2
    return modexp

def frms(d):
    frms = []
    lim = int(floor(sqrt(-d/3.0)))
    for a in range(1, lim+1):
        for b in range(-a, a+1):
            if (b**2-d)%(4*a) == 0:
                c = (b**2-d)/(4*a)
                prim = 0
                if GCD(a, b, c) == 1:
                    prim = 1
                frms.append([a, b, c, prim])
    return frms

def tabulate(n):
    for i in range(1, n+1):
        if i%4 == 0 or i%4 == 3:
            c = 0
            for j in frms(-i):
                if j[-1] == 1:
                    c += 1
            print "h(-", i, ") = ", c, ", # forms = ", len(frms(-i))
    return

def elmntreduce(f, ops):
    while ((f[2] < f[0]) or (f[0] == f[2] and f[1] < 0)) or (f[1] <= -f[0]) or (f[1] > f[0]):
        if (f[2] < f[0]) or (f[0] == f[2] and f[1] < 0):
            ops.append("S")
            f = [f[2], -f[1], f[0]]
        if f[1] <= -f[0]:
            ops.append("T")
            f = [f[0], f[1]+2*f[0], f[0]+f[1]+f[2]]
        if f[1] > f[0]:
            ops.append("T^-1")
            f = [f[0], f[1]-2*f[0], f[0]-f[1]+f[2]]
    return f, ops

def totient(n):
    c = 0
    for i in range(1, n+1):
        if gcd(n, i) == 1:
            c += 1
    return c

def compose(f, g):
    a1 = f[0]
    a2 = g[0]
    b1 = f[1]
    b2 = g[1]
    c1 = f[2]
    d = b1**2-4*a1*c1
    m = bezout(a1, (b1+b2)/2)[0]
    x = bezout(a1, (b1+b2)/2)[1]
    y = bezout(a1, (b1+b2)/2)[2]
    n = gcd(m, a2)
    z = modexp(m/n, totient(a2/n)-1, a2/n)*(x*(b2-b1)/2-c1*y)
    z %= (a2/n)
    A = a1*a2/(n**2)
    B = b1 + 2*a1*z/n
    C = (B**2-d)/(4*A)
    return [A, B, C]

def elmntorder(f):
    order = 0
    g = f
    c = 0
    while g != f or c != 1:
        order += 1
        g = compose(g, f)
        g = elmntreduce(g, [])[0]
        c = 1
    return order
        
def factor(n, l):
    for i in range(2, int(sqrt(n))+1):
        if n%i == 0:
            l.append(i)
            return factor(n/i, l)
    l.append(n)
    return l

def partition(n):
    if n == 0:
        yield []
        return
    #partitions of n-1 used to generate those of n
    for p in partition(n-1):
        yield [1]+p
        #if p is a singleton add [p[0]+1] to list; otherwise if the first two elements of p are strictly increasing, add 1 to first element
        if p and (len(p) < 2 or p[1] > p[0]):
            yield [p[0]+1]+p[1:]

def partorders(part, p):
    orderdict = {}
    gens = []
    for i in part:
        gens.append([j+1 for j in range(p**i)])
    group = [k for k in product(*gens)]
    orders = [len(i) for i in gens]
    for g in group:
        orderg = max(len(gens[i])/(gcd(g[i], len(gens[i]))) for i in range(len(gens)))
        if orderg not in orderdict:
            orderdict[orderg] = 1
        else:
            orderdict[orderg] += 1
    return orderdict

def isprimepower(a, p):
    if a == 1:
        return 1
    while a%p == 0:
        a /= p
    if a == 1:
        return 1
    return 0
    
def class_group(d):
    Forms = frms(d)
    forms = []
    for form in Forms:
        if form[-1] == 1:
            forms.append([form[0], form[1], form[2]])
    order = len(forms)
    factlist = factor(order, [])
    factdict = {}
    for i in factlist:
        if i not in factdict:
            factdict[i] = 1
        else:
            factdict[i] += 1
    orderdict = {}
    if len(forms) == 1:
        return [1, [1]]
    for i in forms:
        if elmntorder(i) not in orderdict:
            orderdict[elmntorder(i)] = 1
        else:
            orderdict[(elmntorder(i))] += 1
    presentation = []
    for prime in factdict:
        primeorderdict = {}
        for i in orderdict:
            if isprimepower(i, prime):
                primeorderdict[i] = orderdict[i]
        for part in partition(factdict[prime]):
            if partorders(part, prime) == primeorderdict:
                presentation.append([prime, part])
    return presentation

def cgp(d):
    pres = class_group(d)
    print "class group for", d, " = "
    for j in range(len(pres)):
        for i in range(len(pres[j][1])):
            if len(pres[j][1][i:]) == 1 and len(pres[j:]) == 1:
                print "C", pres[j][0]**pres[j][1][i]
            else:
                print "C", pres[j][0]**pres[j][1][i], "x"
    return

def amb(d):
    ambs = []
    forms = [[i[0], i[1], i[2]] for i in frms(d)]
    for i in forms:
        if elmntreduce([i[0], -i[1], i[2]], [])[0] == i:
            ambs.append(i)
    return ambs

def isamb(f):
    if elmntreduce([f[0], -f[1], f[2]], [])[0] == f:
        return 1
    return 0

def cgexp(f, exp):
    d = f[1]**2-4*f[0]*f[2]
    if d%4 == 0:
        g = [1, 0, -d/4]
    else:
        g = [1, 1, (1-d)/4]
    while exp != 0:
        if exp%2 == 1: 
            g = compose(g, f)
            g = elmntreduce(g, [])[0]
        f = compose(f, f)
        f = elmntreduce(f, [])[0]
        exp /= 2
    return g

def isprime(p):
    for i in range(2, int(sqrt(p))+1):
        if p%i == 0:
            return 0
    return 1

def plist(B):
    l = []
    for i in range(2, B):
        if isprime(i):
            l.append(i)
    return l

def genform(d):
    while 1:
        a = r.randint(2, int(sqrt(-d)))
        for b in range(int(sqrt(a))):
            if (b**2 - d)%(4*a) == 0:
                c = (b**2-d)/(4*a)
                return elmntreduce([a, b, c], [])[0]

def bqrfact(N, K, B):
    l = plist(B)
    for k in range(2, K+1):
        n = -k*N
        if n%4 == 0 or n%4 == 1:
            for i in range(1, 10):
                f = genform(n)
                print f
                for p in l:
                    exp = p
                    while exp < B:
                        f = cgexp(f, exp)
                        print x
                        if isamb(f):
                            if (f[1] == 0):
                                if f[0] == 1:
                                    fac = gcd(f[2], N)
                                    if fac != N:
                                        return fac
                                        fac = gcd(f[0], N)
                                        if fac != N and fac != 1:
                                            return fac
                            if (f[0] != 1) and (f[0] == f[1]) or (f[0] == f[2]):
                                fac = gcd(f[0], N)
                                if fac != N and fac != 1:
                                    return fac
                        exp *= p
                        for j in range(1, 10):
                            if isamb(f):
                               if (f[1] == 0):
                                   if f[0] == 1:
                                       fac = gcd(f[2], N)
                                       if fac != N:
                                           return fac
                                           fac = gcd(f[0], N)
                                           if fac != N and fac != 1:
                                               return fac
                               if (f[0] != 1) and (f[0] == f[1]) or (f[0] == f[2]):
                                   fac = gcd(f[0], N)
                                   if fac != N and fac != 1:
                                       return fac
                            f = cgexp(f, 2)
    return 0