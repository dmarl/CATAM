from math import *

def gcd(a, b):
    if b > a:
        c = b
        b = a
        a = c
    while a%b != 0:
        c = b
        b = a%b
        a = c
    return b

def modexp(a, b, p):
    #computes a**b mod p with binary modular exponentiation
    modexp = 1
    while b != 0:
        if b%2 == 1:
            modexp = modexp*a%p
        a = a*a%p
        b /= 2
    return modexp

def fact(a, b):
    #highest power of b dividing a
    exp = 0
    while a%b == 0:
        exp += 1
        a /= b
    return exp

def nonqr(p):
    #returns a (minimal) nonqr mod p
    for i in range(2, p):
        if legsym(i, p) == -1:
            return i

def inv(a, p):
    #inverse of a mod p
    return modexp(a, p-2, p)

def legsym(a, p):
    #returns (a/p) using euler's criterion
    if a%p == 0:
        return 0
    legsym = modexp(a, (p-1)/2, p)
    if legsym == p-1:
        return -1
    return 1

def jacsym(a, p):
    #recursive procedure for (a/p)
    if a%p == 0:
        return 0
    a = a%p
    exp = fact(a, 2)
    b = a*2**-exp
    c = 2**exp
    if b == 1:
        return (-1)**(exp*(p**2-1)/8)
    return (-1)**(((p-1)*(b-1)/4) + exp*(p**2-1)/8)*jacsym(p, b)
    
def sqroot(a, p):
    if legsym(a, p) == -1:
        print a, "is not a qr mod", p
        return
    if p%4 == 3:
        x = modexp(a, (p+1)/4, p)
        return x, p-x
    if p%8 == 5:
        x = modexp(a, (p+3)/8, p)
        if x*x%p == a:
            return x, p-x
        else:
            x = x*modexp(2, (p-1)/4, p)%p
            return x, p-x
    exp = fact(p-1, 2)
    s = (p-1)/(2**exp)
    n = nonqr(p)
    b = modexp(n, s, p)
    if modexp(a, (p-1)/4, p) == 1:
        digits = [0]
    else:
        digits = [1]
    for i in range(3, exp+1):
        r = 0
        for j in range(len(digits)):
            r += digits[j]*2**j
        m = modexp(n, r*(p-1)/(2**(i-1)), p)
        m = inv(m, p)
        t = modexp(a, (p-1)/(2**i), p)*m%p
        if t == 1:
            digits.append(0)
        else:
            digits.append(1)
    r = 0
    for j in range(len(digits)):
        r += digits[j]*2**j
    print r
    y = modexp(b, r, p)
    x = modexp(a, (s+1)/2, p)*inv(y, p)%p
    return x, p-x
