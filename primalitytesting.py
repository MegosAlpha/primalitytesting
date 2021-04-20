# Primality Testing
# Library of Primality Testing Algorithms, written for a Number Theory paper.

import random
from math import sqrt
import functools
from sympy.ntheory.residue_ntheory import jacobi_symbol
import sys

# GCD adjustment
sys.setrecursionlimit(1500)

class SieveOfEratosthenes:
    """Implementation of the Lazy / Iterative Version of the Sieve of Eratosthenes"""
    start = 3
    primes = [2]

    # Uses Fibonacci's Observation (sqrt bound)
    def next(self):
        while True:
            for p in self.primes:
                if p > sqrt(self.start):
                    self.primes.append(self.start)
                    self.start += 2
                    return self.primes[-1]
                if self.start % p == 0:
                    self.start += 2
                    break

# Undependable random algorithm
def fermat_test_random(n,k):
    if n == 2:
        return True
    for i in range(k):
        base = random.randint(2,n-1)
        if pow(base,n-1,n) != 1:
                return False
    return True

# Notably, some of the Carmichael Numbers don't fail here.
# This is because it will fail on its factors. In other words,
# the Fermat test becomes a horribly suboptimal means of factorization.
def fermat_test_all_bases(n):
    if n == 2:
        return True
    for i in range(2,n):
        if pow(i,n-1,n) != 1:
            return False
    return True

# Simple recursive GCD
def gcd(a,b):
    if b > a:
        return gcd(b,a)
    if b == 0:
        return a
    return gcd(b, a%b)

# Miller-Rabin
# Can be massively optimized by testing with the first few primes,
# since those bases are sufficient for some quite large composite numbers.
def miller_rabin(n,k):
    if n == 2:
        return True
    bases = set()
    if k >= n-1:
        bases = set(range(1,n))
    else:
        while len(bases) < k-1:
            bases.add(random.randint(1,n-1))
    l = 0
    m = n-1
    while m % 2 == 0:
        m = m//2
        l += 1
    for b in bases:
        fme_results = [pow(b,m*(2**i),n) for i in range(l+1)]
        if fme_results[-1] != 1:
            return False
        for res in fme_results:
            gcdres = gcd(res, n)
            if not (gcdres == 1 or gcdres == n):
                return False
    return True

# Deterministic version of the Miller-Rabin Test
# Not to be confused with the ERH-dependent Miller Test
# This is very slow as numbers get larger
def miller_rabin_deterministic(n):
    if n == 2 or n == 3 or n == 5:
        return True
    bases = range(1,n//2+3)
    l = 0
    m = n-1
    while m % 2 == 0:
        m = m//2
        l += 1
    for b in bases:
        fme_results = [pow(b,m*(2**i),n) for i in range(l+1)]
        if fme_results[-1] != 1:
            return False
        for res in fme_results:
            gcdres = gcd(res, n)
            if not (gcdres == 1 or gcdres == n):
                return False
    return True

# Everything below here is the Baillie-PSW
def miller_rabin_base_2(n):
    l = 0
    m = n-1
    while m % 2 == 0:
        m = m//2
        l += 1
    fme_results = [pow(2,m*(2**i),n) for i in range(l+1)]
    if fme_results[-1] != 1:
        return False
    for res in fme_results:
        if gcd(res, n) != 1 and gcd(res, n) != p:
            return False
    return True

# Iterative form of the Lucas Sequence
# Compute a single term, with caching.
lsU_cache = {}
def lsU(k, P, Q):
    seq = [0,1]
    if Q in lsU_cache.keys():
        seq = lsU_cache[Q]
    while len(seq) < k + 1:
        seq.append(P*seq[-1] - Q*seq[-2])
    lsU_cache[Q] = seq
    return seq[-1]

# Simple implementation of the Lucas-Selfridge Test without noting
# the Strong Lucas Psuedoprime.
def lucas_selfridge(n):
    # Select an element D from the Selfridge sequence 5, -7, 9, -11, 13...
    D = 5
    Q = (1 - D) // 4
    # Second condition expresses n | Q.
    while jacobi_symbol(D,n) != -1 or (Q//n) * n == Q:
        if D < 0:
            D = -D + 2
        else:
            D = -(D + 2)
        Q = (1 - D) // 4
    P = 1
    Q = (1-D)//4
    l = 0
    # There is a stronger variety of this test, but it is slower, and
    # more difficult to implement. As far as I know, there are no
    # known psuedoprimes either way.
    return lsU(n+1, P, Q) % n == 0

# Baillie-PSW
def baillie_psw(n):
    # Special Cases
    if n == 1:
        # 1 is not prime.
        return False
    if n == 2:
        return True
    # Other even numbers
    if n & 1 == 0:
        return False
    # The test directly
    return miller_rabin_base_2(n) and lucas_selfridge(n)
