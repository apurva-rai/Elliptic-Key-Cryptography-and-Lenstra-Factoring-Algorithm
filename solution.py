import math
import random

def xgcd(a, b):
    """
    return (g, x, y) such that a*x + b*y = g = gcd(a, b)
    credit https://en.wikibooks.org/wiki/Algorithm_Implementation/Mathematics/Extended_Euclidean_algorithm
    """
    x0, x1, y0, y1 = 0, 1, 1, 0
    while a != 0:
        (q, a), b = divmod(b, a), a
        y0, y1 = y1, y0 - q * y1
        x0, x1 = x1, x0 - q * x1
    return b, x0, y0

def modinv(a, b):
    """
    return x such that (x * a) % b == 1
    credit https://en.wikibooks.org/wiki/Algorithm_Implementation/Mathematics/Extended_Euclidean_algorithm
    """
    g, x, _ = xgcd(a, b)
    return x % b, g

def elliptic_add(P, Q, a, b, n):
    """
    INPUT
     - Points 'P' = (x1, y1) and 'Q' = (x2, y2)
     - Coefficients 'a' and 'b' defining the elliptic curve y^2 = x^3 + ax + b
     - Modulus of congruence 'n'
    OUTPUT
     - Point R = P + Q
    """
    #First, check to see if either point is the identity; return the other if so
    if (P == 0):
        return Q
    elif (Q == 0):
        return P

    #Neither point is the identity; get the coordinate pair of each
    try:
        x1, y1 = P
        x2, y2 = Q
    except TypeError:
        return None

    #We need to make sure these points are actually solutions to the given elliptic curve;
    #we implement this for you. The function will throw an error if this does not hold
    assert((y1**2 - x1**3 - a*x1 - b) % n == 0),"P is not a point on the curve"
    assert((y2**2 - x2**3 - a*x2 - b) % n == 0),"Q is not a point on the curve"
    num = den = 0

    if (P == Q):
        if (2*y1%n == 0):
            return (0,0)
        else:
            num = (3*(x1**2) + a)%n
            den = (2*y1)%n
    else:
        if ((x2-x1)%n == 0):
            return 0
        else:
            num = (y2-y1)%n
            den = (x2- x1)%n

    if math.gcd(den,n) != 1:
        return None
    else:
        den = modinv(den,n)[0]

    slope = (num * den) % n

    xR = (pow(slope,2,n) - x1 - x2) % n
    yR = (slope*(x1- (pow(slope,2,n) - x1 - x2)) - y1) % n

    return (xR,yR)

def elliptic_mul(d, P, a, n):
    """
    INPUT
     - Positive integer 'd' with which to multiply 'P'
     - Point 'P'
     - Coefficient 'a' which when taken with 'P' and 'n' uniquely determines the coefficient 'b' defining the elliptic curve y^2 = x^3 + ax + b
     - Modulus of congruence 'n'
    OUTPUT
     - Point d*P
    """
    #Fixing coefficient 'a' and point 'P', we compute constant term 'b' for the curve
    if (P != 0):
        b = (P[1]**2 - P[0]**3 - a*P[0]) % n

    R = P
    while d > 1:
        R = elliptic_add(P, R, a, b, n)
        d-= 1

    return R

"""
All versions of the Lenstra algorithm computes successive multiples of a point P. The algorithms
differ in which successive multiples they use. In particular, on iteration i, we compute f(i)P, where
f(i) is the function giving the multiple to use on iteration i.
"""
def lenstra(P, a, n):
    """
    INPUT
     - Point 'P'
     - Coefficient 'a' which when taken with 'P' and 'n' uniquely determines the coefficient 'b' defining the elliptic curve y^2 = x^3 + ax + b
     - Modulus of congruence 'n'
    OUTPUT
     - Computed factor of 'n'

    Hint: One possible choice for f(i) is f(i) = i, that is, at step i we calculate i!*P
    """
    i = 0
    while (i < n):
        tempR = elliptic_mul(math.factorial(i), P, a, n)
        i += 1
        if(tempR == None):
            return n
        try:
            if(math.gcd((2*tempR[1])%n,n) != 1):
                print(i,math.factorial(i))
                return (math.gcd(2*tempR[1],n)%n)
        except TypeError:
            return n


def rand_elliptic(n):
    """
    Generates a random elliptic curve y^2 = x^3 + ax + b
    and a point P0 on that curve
    """
    state = False
    while(state != True):
        P0 = random.randint(0, n - 1), random.randint(0, n - 1)
        a = random.randint(0, n - 1)
        if (P0 != 0):
            b = (P0[1]**2 - P0[0]**3 - a*P0[0]) % n
        state = ((P0[1]**2)%n == (P0[0]**3)%n + ((a*P0[0])%n) + b)
    return P0, a


def lenstra_random(n):
    """
    After you define the `lenstra` function above, this function will run until a proper factor is found by generating random curves
    """
    #Set g to n to guarantee entry into the loop
    g = n

    #Repeat until g is not n(i.e., a "proper factor" of n)
    while (g == n):
        P0, a = rand_elliptic(n)
        g = lenstra(P0, a, n)
        b = (P0[1]**2 - P0[0]**3 - a*P0[0]) % n

        print(P0,a,b)

    #Temporary fix to prevent error in the test cases
    #return n if g is None else g
    #You will want to comment out the above line and uncomment the below line
    return g

print(lenstra_random(978361))
