import json
from fractions import Fraction
from math import factorial
P=json.load(open('pieces.json'))
K,T,S=P['K'],P['T'],P['S']
def prod(f,n):
    r=1
    for k in range(1,n+1): r*=f(k)
    return r
def S_closed(n):
    return Fraction(9*n+16,16)*prod(lambda k: Fraction(6*(3*k-1)*(3*k+1),(k+2)*(2*k+3)), n)
def T_closed(n):
    return Fraction(7*n+5,5)*prod(lambda k: Fraction(6*(2*k-1)*(6*k-1)*(6*k+1),(k+1)*(4*k+3)*(4*k+5)), n)
print("S product form (9n+16)/16 * ... verified:", all(S_closed(n)==S[n] for n in range(47)))
print("T product form (7n+5)/5  * ... verified:", all(T_closed(n)==T[n] for n in range(47)))

def S_fact(n):
    return Fraction((9*n+16)*2**(2*n+2)*3*factorial(3*n+1)*factorial(n+1),
                    16*factorial(n)*factorial(n+2)*factorial(2*n+3))
print("S factorial form  (9n+16)*3*4^(n+1)*(3n+1)!*(n+1)! / (16*n!*(n+2)!*(2n+3)!):",
      all(S_fact(n)==S[n] for n in range(47)))
