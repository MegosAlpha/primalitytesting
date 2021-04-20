# Implementation of the Agrawal-Kayal-Saxena Test
# Note: Only works with SageMath. Kept as .py for Syntax Highlighting, but change to
# .sage or copy to a SageMath notebook for actual use
from sympy import perfect_power
def order(r,n):
    k = 1
    v = pow(n,k,r)
    while v != 1:
        k += 1
        v = pow(n,k,r)
    return k


def aks(n):
    if n % 2 == 0:
        return n==2
    if perfect_power(n) != False:
        return False
    r = 2
    l2n = log(n,2).n()
    ordx = order(r,n)
    while ordx <= l2n**2:
        r += 1
        # SageMath has built-in GCD
        if gcd(r,n) != 1:
            if r < n:
                return False
            else:
                continue
        ordx = order(r,n)
    if n <= r:
        return True
    # This is where SageMath is important
    zn = PolynomialRing(Integers(n),'v')
    v = zn.gen()
    znbtr = zn.quotient(v^r - 1, 'x')
    x = znbtr.gen()
    for a in range(1, floor(sqrt(euler_phi(r))*l2n)):
        if (x+a)^n != (x^n + a):
            return False
    return True
