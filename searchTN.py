# -*- coding: utf-8 -*-
"""
Created on Sun Dec  8 10:02:17 2024

@author: erkanuslu
"""
from __future__ import division
from __future__ import print_function
from scipy.stats import binom

import random
import functools
import pandas as pd
import numpy as np
import math
import itertools


def _eval_at(poly, x, prime):
    """Evaluates polynomial (coefficient tuple) at x, used to generate a
    shamir pool in make_random_shares below.
    """
    accum = 0
    for coeff in reversed(poly):
        accum *= x
        accum += coeff
        accum %= prime
    return accum

def _eval_at_unreversed(poly, x, prime):
    """Evaluates polynomial (coefficient tuple) at x, used to generate a
    shamir pool in make_random_shares below.
    """
    accum = 0
    for coeff in poly:
        accum *= x
        accum += coeff
        accum %= prime
    return accum

def _extended_gcd(a, b):
    """
    Division in integers modulus p means finding the inverse of the
    denominator modulo p and then multiplying the numerator by this
    inverse (Note: inverse of A is B such that A*B % p == 1). This can
    be computed via the extended Euclidean algorithm
    http://en.wikipedia.org/wiki/Modular_multiplicative_inverse#Computation
    """
    x = 0
    last_x = 1
    y = 1
    last_y = 0
    while b != 0:
        quot = a // b
        a, b = b, a % b
        x, last_x = last_x - quot * x, x
        y, last_y = last_y - quot * y, y
    return last_x

def _divmod(num, den, p):
    """Compute num / den modulo prime p

    To explain this, the result will be such that:
    den * _divmod(num, den, p) % p == num
    """
    
    inv= _extended_gcd(den, p)
    if (inv * den) % p != 1:
        print ("hata", den, "   ", inv )
    
    return int(int(num%p) * int(inv%p) % p)

def _lagrange_interpolate(x, x_s, y_s, p):
    """
    Find the y-value for the given x, given n (x, y) points;
    k points will define a polynomial of up to kth order.
    """
    k = len(x_s)
    assert k == len(set(x_s)), "points must be distinct"
    def PI(vals):  # upper-case PI -- product of inputs
        accum = 1
        for v in vals:
            accum *= v
        return accum
    nums = []  # avoid inexact division
    dens = []
    for i in range(k):
        others = list(x_s)
        cur = others.pop(i)  
        nums.append(PI(x - o for o in others))
        dens.append(PI(cur - o for o in others))
    den = PI(dens)
    num = sum([_divmod(nums[i] * den * y_s[i] % p, dens[i], p)
               for i in range(k)])
    return (_divmod(num, den, p) + p) % p

def _lagrange_interpolate_shares(x, shares, p):
    """
    Find the y-value for the given x, given n (x, y) points;
    k points will define a polynomial of up to kth order.
    """
    x_s, y_s = zip(*shares)
    k = len(x_s)
    assert k == len(set(x_s)), "points must be distinct"
    def PI(vals):  # upper-case PI -- product of inputs
        accum = 1
        for v in vals:
            accum *= v
        return accum
    nums = []  # avoid inexact division
    dens = []
    for i in range(k):
        others = list(x_s)
        cur = others.pop(i)  
        nums.append(PI(x - o for o in others))
        dens.append(PI(cur - o for o in others))
    den = PI(dens)
    num = sum([_divmod(nums[i] * den * y_s[i] % p, dens[i], p)
               for i in range(k)])
    return (_divmod(num, den, p) + p) % p

def _lagrange_interpolate_polynomial(x_s, y_s, p):
    """
    Find the y-value for the given x, given n (x, y) points;
    k points will define a polynomial of up to kth order.
    """
    k = len(x_s)
    assert k == len(set(x_s)), "points must be distinct"
    def PI(vals):  # upper-case PI -- product of inputs
        accum = 1
        for v in vals:
            accum *= v
        return accum
    nums = []  # avoid inexact division
    dens = []
    for i in range(len(x_s)):
        others = list(x_s)
        cur = others.pop(i)  
        temp_num = []
        for j in others:
            temp_num.append([1,-j])
        nums.append(poly_accum(temp_num,p))
        dens.append(PI(cur - o for o in others))
        
    #den = PI(dens) 
    #dens = np.array(dens) 
    #dens = dens % p
    #print (dens)
    for i in range(len(x_s)):
        nums[i] = (nums[i]*y_s[i]) % p
    result = _divmodpoly(nums,dens,p)
    result_poly = [0 for i in range(len(result[0]))]
    for i in range(len(result)):
        result_poly += result[i] % p
    return result_poly % p

def recover_secret(shares, prime):
    """
    Recover the secret from share points
    (points (x,y) on the polynomial).
    """
    if len(shares) < 3:
        raise ValueError("need at least three shares")
    x_s, y_s = zip(*shares)
    return _lagrange_interpolate(0, x_s, y_s, prime)  
def recover_polynomial(shares, prime):
    """
    Recover the secret from share points
    (points (x,y) on the polynomial).
    """
    if len(shares) < 3:
        raise ValueError("need at least three shares")
    x_s, y_s = zip(*shares)
    return _lagrange_interpolate_polynomial(x_s, y_s, prime)  
 
    def PI(vals):  # upper-case PI -- product of inputs
        accum = 1
        for v in vals:
            accum *= v 
        return accum
def poly_accum(poly,p):  # upper-case PI -- product of inputs
    poly = np.array(poly)
    result = poly[0]
    for v in range(1,len(poly)):
        result = np.convolve(result,poly[v])
    return result
def _divmodpoly(numPoly, denPoly, p):
    resultPoly = []
    for i in range(len(numPoly)):
        result = []
        for j in range(len(numPoly[0])):
            print(i, "    ", j, "   ", denPoly[i])
            result.append( (_divmod(numPoly[i][j], denPoly[i], p) + p) % p)
        resultPoly.append(result)    
    return np.array(resultPoly,dtype = np.int64)


nums = []  # avoid inexact division
dens = []
for i in range(3):
    others = list(x_s)
    cur = others.pop(i)  
    temp_num = []
    for j in others:
        temp_num.append([1,-j])
    nums.append(poly_accum(temp_num,p))
    dens.append(PI(cur - o for o in others))
den = PI(dens) 
for i in range(len(x_s)):
    nums[i] = nums[i]*y_s[i]
result = _divmodpoly(nums,dens,p)
result_poly = [0 for i in range(len(result[0]))]
for i in range(len(result)):
    result_poly += result[i] % p
result_poly = result_poly % p

p = 8380417
tempCoefs = {0,1,2,p-2, p-1}
result11q0 = []
result11q1 = []
result11q2 = []
result11q3 = []
for i1 in tempCoefs:
    for i2 in tempCoefs:
        for i3 in tempCoefs:
            for i4 in tempCoefs:
                for i5 in tempCoefs:
                    for i6 in tempCoefs:
                        for i7 in tempCoefs:
                            for i8 in tempCoefs:
                                for i9 in tempCoefs:
                                     for i10 in tempCoefs:
                                         for i11 in tempCoefs:
                                             temp = [[1,i1], [2,i2], [3,i3], [4,i4], [5,i5], [6,i6], [7,i7], [8,i8], [9,i9], [10,i10],[11,i11]]
                                             recovered = recover_secret(temp,p)
                                             if recovered in tempCoefs:
                                                 result11q0.append([recovered, [i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11]]) 
                                                 point1 = _lagrange_interpolate_shares(12,temp,p)
                                                 if point1 in tempCoefs:
                                                     result11q1.append([recovered, [i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11]]) 
                                                     point2 = _lagrange_interpolate_shares(13,temp,p)
                                                     if point2 in tempCoefs:
                                                         result11q2.append([recovered, [i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11]]) 
                                                         point3 = _lagrange_interpolate_shares(14,temp,p)
                                                         if point3 in tempCoefs:
                                                             result11q3.append([recovered, [i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11]]) 
                                             
                                                
                                                
                                                point2 = _lagrange_interpolate_shares(12,temp,p)
                                             point3 = _lagrange_interpolate_shares(13,temp,p)
                                             if point1 in tempCoefs and point2 in tempCoefs and point3 in tempCoefs:
                                                result10q3.append([recovered, [i1,i2,i3,i4,i5,i6,i7,i8,i9,i10]]) 
                               
                                    result8q.append([recover_polynomial(temp,p), recovered, [i1,i2,i3,i4,i5,i6,i7,i8]]) 
                                       if recovered in tempCoefs:
                                           if _eval_at_unreversed(recover_polynomial(temp,p),8,p) in tempCoefs: 
                                               result8q.append([recover_polynomial(temp,p), recovered, [i1,i2,i3,i4,i5,i6,i7,i8]])   
num = sum([_divmod(nums[i] * y_s[i] % p, dens[i], p)
            for i in range(k)])
 