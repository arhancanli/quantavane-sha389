# Pure Python point count (no PARI): supersingular primes of 389a1 (y^2 + y = x^3 + x^2 - 2x) below a bound.
# a_p = -sum_x legendre(4x^3 + 4x^2 - 8x + 1, p) for odd p != 389. Supersingular for p >= 5 iff a_p = 0.
import sys
B = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
def primes(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]
def ap(p):
    chi = [-1] * p; chi[0] = 0
    for y in range(1, (p + 1) // 2): chi[y * y % p] = 1
    return -sum(chi[(4*x*x*x + 4*x*x - 8*x + 1) % p] for x in range(p))
extra = [int(a) for a in sys.argv[2:]]
ss = [p for p in primes(B) if p >= 5 and p != 389 and ap(p) == 0]
print("ss primes <", B, ":", ss)
for p in extra: print("check", p, "a_p =", ap(p))
