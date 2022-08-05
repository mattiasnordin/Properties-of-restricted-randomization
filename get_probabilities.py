#This code has been tested with Python version 3.10.4
from sympy import symbols, binomial, simplify #Tested with version 1.10.1
from itertools import permutations
import numpy as np #Tested with version 1.22.3


def gen_pr_name(inp, nr):
    out = 'p_'
    for i in inp:
        if i == -1:
            out += '0'
        else:
            out += '1'
    if nr == 2:
        out += "'"
    elif nr == 3:
        out += "''"
    elif nr == 4:
        out += "'''"
    return out


def combs_eq(ab, nr):
    '''Generate different versions of equation A27 (nr=2),
    equation A28 (nr=3) or equation A29 (nr=4)'''
    c1 = 0
    c2 = 0
    c3 = 0
    for i in range(nr):
        c1 += ab[0][i] == ab[1][i] == 1
        c2 += ab[0][i] > ab[1][i]
        c3 += ab[0][i] < ab[1][i]
        
    return (binomial(N-nr, N2-c1-c2) * 
            binomial(N2-c1-c2, N2-u-c1) * 
            binomial(N-nr-N2+c1+c2, u-c3) /
            (binomial(N, N2) * binomial(N2, N2-u) * binomial(N2, u)))


def get_simplify2(ai, bi, aj, bj):
    nr = 2
    out = list(set(permutations(
        np.repeat([1, -1], [nr, nr], axis=0), nr)))
    out2 = list(permutations(out, 2))
    for i in out:
        out2.append((i, i))

    v = 0

    for i in out2:
        q = combs_eq(i, nr)
        if i[0][0] == ai and i[0][1] == aj and i[1][0] == bi and i[1][1] == bj:
            v += q
    
    return  simplify(v)


def get_simplify3(ai, aj, bi, bk):
    nr = 3
    out = list(set(permutations(
        np.repeat([1, -1], [nr, nr], axis=0), nr)))
    out2 = list(permutations(out, 2))
    for i in out:
        out2.append((i, i))

    v = 0

    for i in out2:
        q = combs_eq(i, nr)
        if i[0][0] == ai and i[0][1] == aj and i[1][0] == bi and i[1][2] == bk:
            v += q
    
    return  simplify(v)


def get_simplify4(ai, aj, bk, bl):
    nr = 4
    out = list(set(permutations(
        np.repeat([1, -1], [nr, nr], axis=0), nr)))
    out2 = list(permutations(out, 2))
    for i in out:
        out2.append((i, i))

    v = 0

    for i in out2:
        q = combs_eq(i, nr)
        if i[0][0] == ai and i[0][1] == aj and i[1][2] == bk and i[1][3] == bl:
            v += q
    
    return  simplify(v)


N, u= symbols('N u')
N2 = N/2

inp = list(set(permutations(np.repeat([1, -1], [4, 4], axis=0), 4)))
out = {}
for i in inp:
    out[gen_pr_name(i, 2)] = get_simplify2(i[0], i[1], i[2], i[3])
    out[gen_pr_name(i, 3)] = get_simplify3(i[0], i[1], i[2], i[3])
    out[gen_pr_name(i, 4)] = get_simplify4(i[0], i[1], i[2], i[3])

#Print all probabilities
for i, j in out.items():
    print(i, j)