#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 18 20:42:42 2025

@author: erkanuslu
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 23:12:53 2024

@author: erkanuslu
"""

from __future__ import division
from __future__ import print_function

import random
import functools
import numpy as np
from sympy import Matrix
from scipy.linalg import lu
import hashlib



# 12th Mersenne Prime
_PRIME = 7

_RINT = functools.partial(random.SystemRandom().randint, 0)

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

def make_random_shares(secret, minimum, shares, prime=_PRIME):
    """
    Generates a random shamir pool for a given secret, returns share points.
    """
    if minimum > shares:
        raise ValueError("Pool secret would be irrecoverable.")
    poly = [secret] + [_RINT(prime - 1) for i in range(minimum - 1)]
    points = [(i, _eval_at(poly, i, prime))
              for i in range(1, shares + 1)]
    return points

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
    return last_x, last_y

def _divmod(num, den, p):
    """Compute num / den modulo prime p

    To explain this, the result will be such that:
    den * _divmod(num, den, p) % p == num
    """
    inv, _ = _extended_gcd(den, p)
    return num * inv

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

def recover_secret(shares, prime=_PRIME):
    """
    Recover the secret from share points
    (points (x,y) on the polynomial).
    """
    if len(shares) < 3:
        raise ValueError("need at least three shares")
    x_s, y_s = zip(*shares)
    return _lagrange_interpolate(0, x_s, y_s, prime)


def convol(a,b,prime):
    res = []
    for i in range(4):
        for j in range(4):
            res.append(a[i]*b[j])
    res[1] += res[4]
    res[2] += res[8]
    res[3] += res[12]
    res[6] += res[9]
    res[7] += res[13]
    res[11] += res[14]
    res = np.delete(res,[4,8,12,9,13,14])
    return (res %  prime)  
    """ x1^2 + x1x2 + x1x3 + x1x4 + x2^2 + x2x3 + x2x4 + x3^2 + x3x4 + x4^2 """



def T_convol(y):
    return [convol(y[0], y[2],7),convol(y[0], y[3],7),convol(y[1], y[2],7),convol(y[1], y[3],7)]
    """ x1x3 + x1x4 + x2x3 + x2x4"""

def make_random_polynomial_shares(poly, minimum, shares, prime):
    result = list()
    for i in range(len(poly)):
        result.append(make_random_shares(poly[i], minimum, shares, prime))
    result_polys = list()
    for i in range(len(result[0])):
        temp = list()
        for j in range(len(result)):
            temp.append([i+1, result[j][i][1]%prime])
        result_polys.append(temp)    
    return result_polys
def recover_polynomial_secret(shares_poly, threshold, prime):
    result_poly = list()
    for i in range(len(shares_poly[0])):
        temp = list()
        for j in range(threshold):
            temp.append(shares_poly[j][i])
        result_poly.append(recover_secret(temp, prime))   
    return result_poly  



def linear_transformation(arr,convol,prime):
    result = [0 for element in range(len(convol[0]))]
    for i in range(len(convol)):
        result += arr[i]*convol[i]
    return result % prime    

def linear_transformation_shares(polyShares,y,prime):
    result_poly = list()
    for i in range(len(polyShares)):
        temp = [0 for element in range(len(y[0]))]
        for j in range(len(polyShares[0])):
            temp += polyShares[i][j][1] * y[j]
        temp2 = list()
        for j in range(len(temp)):
            temp2.append([i+1, temp[j] % prime])
        
        result_poly.append(temp2)
            
    return result_poly       
         
def eval_at(f,x3,x4,prime):
    result_poly = []
    for j in range(2):
        result_poly.append((f[2*j] * x3 + f[(2*j)+1] * x4) % prime)
    return result_poly   
def eval_at_shares(f_shares,x3x4,prime):
    result_poly = list()
    for i in range(len(f_shares)):
        temp = []
        for j in range(2):
            temp.append((f_shares[i][2*j][1] * x3x4[0] + f_shares[i][(2*j)+1][1] * x3x4[1]) % prime)
        result_poly.append(temp)
    
    return result_poly   

def eval_at_shares_x3x4(f_shares,x3x4_shares,prime):
    result_poly = list()
    for i in range(len(f_shares)):
        temp = []
        for j in range(2):
            temp.append([f_shares[i][0][0],(f_shares[i][2*j][1] * x3x4_shares[i][0][1] + f_shares[i][(2*j)+1][1] * x3x4_shares[i][1][1]) % prime])
        result_poly.append(temp)
    
    return result_poly   

def eval_shared_messages(f1_share,f2_share,sign_share,prime):
    result = []
    for i in range(len(f1_share)):
        temp = [(f1_share[i][0][1]*sign_share[i][0][1]*sign_share[i][2][1] + f1_share[i][1][1]*sign_share[i][0][1]*sign_share[i][3][1] + f1_share[i][2][1]*sign_share[i][1][1]*sign_share[i][2][1] + f1_share[i][3][1]*sign_share[i][1][1]*sign_share[i][3][1])%prime
                , (f2_share[i][0][1]*sign_share[i][0][1]*sign_share[i][2][1] + f2_share[i][1][1]*sign_share[i][0][1]*sign_share[i][3][1] + f2_share[i][2][1]*sign_share[i][1][1]*sign_share[i][2][1] + f2_share[i][3][1]*sign_share[i][1][1]*sign_share[i][3][1])%prime]
        result.append(temp)
    return result   
            

def modular_inverse(a, mod):
    # Genişletilmiş Öklid algoritması ile modüler ters bulma
    t, new_t = 0, 1
    r, new_r = mod, a
    while new_r != 0:
        quotient = r // new_r
        t, new_t = new_t, t - quotient * new_t
        r, new_r = new_r, r - quotient * new_r
    if r > 1:
        raise ValueError(f"{a} has no modular inverse under mod {mod}")
    if t < 0:
        t += mod
    return t

def solve_modular_system(coeff_matrix, result_vector, mod):
    coeff_matrix = np.array(coeff_matrix, dtype=int)
    result_vector = np.array(result_vector, dtype=int)
    n = len(result_vector)

    # Gauss eliminasyon uygulayarak çözüm bulma
    for i in range(n):
        # i. sütunda sıfır olmayan en büyük elemanı bul
        if coeff_matrix[i, i] == 0:
            for j in range(i + 1, n):
                if coeff_matrix[j, i] != 0:
                    # Satır değişimi yaparak sıfır olmayan bir eleman getir
                    coeff_matrix[[i, j]] = coeff_matrix[[j, i]]
                    result_vector[[i, j]] = result_vector[[j, i]]
                    break
            else:
                raise ValueError("Sistemde çözüm yok veya çözüm tekil.")
        
        # Modüler tersini al ve kendini bir yap
        inv = modular_inverse(coeff_matrix[i, i], mod)
        coeff_matrix[i] = (coeff_matrix[i] * inv) % mod
        result_vector[i] = (result_vector[i] * inv) % mod

        # Diğer satırlardan bu elemanı çıkar
        for j in range(n):
            if i != j:
                factor = coeff_matrix[j, i]
                coeff_matrix[j] = (coeff_matrix[j] - factor * coeff_matrix[i]) % mod
                result_vector[j] = (result_vector[j] - factor * result_vector[i]) % mod

    return result_vector




def matrix_shares_mult(matrix,vector_shares,prime):
    result = []
    for i in range(len(vector_shares)):
        temp = []
        for j in range(len(vector_shares[0])):
            temp.append(vector_shares[i][j][1])
        temp2 = np.dot(matrix,temp)
        temp3 = []
        for z in range(len(temp2)):
            temp3.append([i+1,temp2[z]%prime])
        result.append(temp3)    
    return result   

def concat(a,b):
    return  int(str(a) + str(b))


def matrix_vector_mod_product(matrix, vector, mod):
    """
    4x4 matris ve 4x1 vektörün mod işlemine göre çarpımını yapar.
    
    Args:
        matrix (list of list of int): 4x4 boyutunda matris.
        vector (list of int): 4x1 boyutunda vektör.
        mod (int): Mod değeri.
    
    Returns:
        list of int: Mod işlemine göre sonuçlanan 4x1 vektör.
    """
    # Matris ve vektörün NumPy dizisine dönüştürülmesi
    matrix = np.array(matrix)
    vector = np.array(vector)
    
    # Matris ve vektör çarpımı
    result = np.dot(matrix, vector)
    
    # Mod işlemi
    result_mod = result % mod
    
    # Listeye dönüştürülmesi
    return result_mod


def hash_to_two_integers(input_data, mod_value):
    
    
    if isinstance(input_data, str):
        bytes_data = input_data.encode()  # Stringi bytes'a dönüştür
    elif isinstance(input_data, (int, float)):
        bytes_data = str(input_data).encode()  # Sayıyı string yap ve encode et
    elif isinstance(input_data, (list, dict, set)):
        bytes_data = json.dumps(input_data).encode()  # Yapıları JSON stringine çevir ve encode et
    elif isinstance(input_data, bytes):
        bytes_data = input_data  # Eğer zaten bytes ise doğrudan kullan
    else:
        raise TypeError("Desteklenmeyen veri tipi!")
    
    # SHA-256 hash oluştur


    
    # Verilen veriyi SHA-256 ile hashle
    hash_output = hashlib.sha256(bytes_data).hexdigest()
    
    # Hash çıktısını ikiye böl
    midpoint = len(hash_output) // 2
    part1 = hash_output[:midpoint]
    part2 = hash_output[midpoint:]
    
    # Her iki parçayı da integer'a çevir
    int1 = int(part1, 16) % mod_value
    int2 = int(part2, 16) % mod_value
    
    return [int1, int2]

def eval_verification(sign,pk,q):
    return ((pk[0]*sign[0]*sign[0])+(pk[1]*sign[0]*sign[1])+(pk[2]*sign[0]*sign[2])+(pk[3]*sign[0]*sign[3]) +(pk[4]*sign[1]*sign[1])+(pk[5]*sign[1]*sign[2])+(pk[6]*sign[1]*sign[3])+(pk[7]*sign[2]*sign[2])+(pk[8]*sign[2]*sign[3])+(pk[9]*sign[3]*sign[3]))%q 
def eval_verification_Shares(sign,pk_Shares,q):
    x_s, pk = zip(*pk_Shares)
    return ((pk[0]*sign[0]*sign[0])+(pk[1]*sign[0]*sign[1])+(pk[2]*sign[0]*sign[2])+(pk[3]*sign[0]*sign[3]) +(pk[4]*sign[1]*sign[1])+(pk[5]*sign[1]*sign[2])+(pk[6]*sign[1]*sign[3])+(pk[7]*sign[2]*sign[2])+(pk[8]*sign[2]*sign[3])+(pk[9]*sign[3]*sign[3]))%q 
def vt_uov_Keygen(f1,f2,T,q):
    
    
    T = T_convol(T)
    f1_T = linear_transformation(f1,T,q)
    f2_T = linear_transformation(f2,T,q)
    f1_Shares = make_random_polynomial_shares(f1,3,6,q)
    f2_Shares = make_random_polynomial_shares(f2,3,6,q)
    f1_T_Shares = linear_transformation_shares(f1_Shares,T,q)
    f2_T_Shares = linear_transformation_shares(f2_Shares,T,q)
    sk_shares=[f1_Shares,f2_Shares]
    pk_shares = [f1_T_Shares,f2_T_Shares]
    
    return sk_shares,pk_shares

def vt_uov_Sign(f1,f2,sk_shares,T,M,salt,q):
    x3x4 = [2,3]
    T_inv = []
    A = Matrix(T) # keyMatrix is your basic matrix ndrarray format
    T_inv.append(np.array(A.inv_mod(q)))
    t = hash_to_two_integers(concat(M,salt),q)
    f1_x3x4 = eval_at(f1,x3x4[0],x3x4[1],q)
    f2_x3x4 = eval_at(f2,x3x4[0],x3x4[1],q)
    s = solve_modular_system([f1_x3x4,f2_x3x4], t, q)
    s = np.append(s, x3x4)
    s_inv = matrix_vector_mod_product(T_inv[0], s, 7)
    sign = [s_inv,salt]
    
    salt_Shares = make_random_shares(salt, 3, 6, prime=8380417)
    t_Shares = []
    for i in range(len(salt_Shares)):
        t_Shares.append(hash_to_two_integers(concat(M,salt_Shares[i][1]),q))
    f1_shares_x3x4 = eval_at_shares(sk_shares[0], x3x4, q)
    f2_shares_x3x4 = eval_at_shares(sk_shares[1], x3x4, q)
    sign_Shares = []
    for i in range(len(f1_shares_x3x4)):
        sTemp = solve_modular_system([f1_shares_x3x4[i],f2_shares_x3x4[i]], t_Shares[i], q)
        sTemp = np.append(sTemp,x3x4)
        sTemp = matrix_vector_mod_product(T_inv, sTemp, 7)
        sign_Shares.append([sTemp,salt_Shares[i][1]])
    return sign , sign_Shares, t, t_Shares

def vt_uov_Verify(sign_Shares, M, pk_Shares,q):
    t_verification_Shares = []
    for i in range(len(sign_Shares)):
        t_verification_Shares.append(hash_to_two_integers(concat(M,sign_Shares[i][1]),q))
    t_prime_Shares = []
    for i in range(len(sign_Shares)):
        temp = [eval_verification_Shares(sign_Shares[i][0][0],pk_Shares[0][i],q),eval_verification_Shares(sign_Shares[i][0][0],pk_Shares[1][i],q)]
        t_prime_Shares.append(temp)
    return   t_verification_Shares == t_prime_Shares   
def main():
    """Main function"""
    
    q = 7
    f1 = [4,6,2,1]
    f2 = [2,6,0,3] 
    T = [[1,2,4,1],[2,3,0,6],[4,4,2,1],[6,3,1,1]]
    T_inv = []
    A = Matrix(T) # keyMatrix is your basic matrix ndrarray format
    T_inv.append(np.array(A.inv_mod(q)))
   
    M = 1423
    salt = 123
   
    sk_Shares, pk_Shares = vt_uov_Keygen(f1, f2, T, 7)
    sign , sign_Shares, t ,t_Shares = vt_uov_Sign(f1,f2,sk_Shares,T,M,salt,q)
    vt_uov_Verify(sign_Shares, M, pk_Shares,q)
   
    
    

if __name__ == '__main__':
    main()
