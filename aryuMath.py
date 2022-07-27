import math


def binomial_coefficients(a, b):
    """Get binomial coefficients aCb"""
    n, m = a, b
    if n < 2 * m:
        m = n - m
    if n == m or m == 0:
        return 1
    tmp = [n - i for i in range(m)]
    for i in range(1, m + 1):
        for j in range(m):
            if not tmp[j] % i:
                tmp[j] //= i
                break
    result = 1
    for i in tmp:
        result *= i
    return result


def farey(N):
    """Get order N Farey series"""
    res = [[[0, 1], [1, 1]]]
    for idx in range(N - 1):
        tmp = []
        for f in range(len(res[idx]) - 1):
            bunshi = res[idx][f][0] + res[idx][f + 1][0]
            bunbo = res[idx][f][1] + res[idx][f + 1][1]
            GCD_Num = math.gcd(bunshi, bunbo)
            bunshi = bunshi // GCD_Num
            bunbo = bunbo // GCD_Num
            tmp.append(res[idx][f])
            if bunbo > idx + 2:
                continue
            tmp.append([bunshi, bunbo])
        tmp.append([1, 1])
        res.append(tmp)
    return res


def bernoulli_q(num):
    """calc Bernoulli Number In Q on 0 ~ N"""
    res = [[1, 1]]
    for i in range(1, num + 1):
        numer, denom = res[0][0] * binomial_coefficients(i + 1, 0), res[0][1]
        for j in range(1, i):
            numer = numer * res[j][1] + denom * res[j][0] * binomial_coefficients(
                i + 1, j
            )
            denom = denom * res[j][1]
            if numer:
                g = math.gcd(numer, denom)
                numer //= g
                denom //= g
        numer *= -1
        denom *= i + 1
        if (numer < 0 and denom < 0) or (numer > 0 and denom < 0):
            numer = -numer
            denom = -denom
        if not numer:
            denom = 1
        else:
            g = math.gcd(abs(numer), abs(denom))
            numer //= g
            denom //= g
        res.append([numer, denom])
    return res


def sqrt_continued_fraction(num):
    """Get Continued fraction expansion of sqrt(num)"""
    res = []
    rootnum = int(math.sqrt(num))
    ak = rootnum
    res.append(ak)
    if ak**2 == num:
        return res
    Pk = 0
    Qk = 1
    first = True
    while True:
        Pk1 = ak * Qk - Pk
        Qk1 = (num - Pk1**2) // Qk
        if first:
            TP = Pk1
            TQ = Qk1
            first = False
        elif TP == Pk1 and TQ == Qk1:
            break
        ak = (Pk1 + rootnum) // Qk1
        Pk = Pk1
        Qk = Qk1
        res.append(ak)
    return res
