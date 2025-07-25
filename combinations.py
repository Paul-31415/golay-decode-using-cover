def enumerate_combs(n,ones):
    v = ((1<<ones)-1)<<(n-ones)
    while v >> ones:
        yield v
        lsb = ((v^(v-1))+1)>>1
        if lsb > 1:
            v = v-lsb+(lsb>>1)
        else:
            lcb = ((v^(v+1))+1)>>1
            m = lcb-1
            v &= ~m
            lsb = ((v^(v-1))+1)>>1
            v = v-lsb+(lsb>>1)
            v += m*(lsb//lcb)>>1
    yield v
import math
def nth_comb(n,length,ones,comb=None):
    if length == 1:
        return (ones>0)*1
    comb = math.comb(length,ones) if comb is None else comb
    z = (comb * ones)//length
    return (nth_comb(n-z*(n>=z),length-1,ones-(n<z),z if n<z else comb-z)<<1) | (n<z)

    # len l, set s
    # if we put a one, there are l-1 choose s-1 rest
    #            else:           l-1 choose s
    #            a ncr b =  a!/(b!(a-b)!)      dec a: * (a-b)/a
    #                                        dec a&b: * b/a
    #
    #                     (l-1)!/((s-1)!(l-1-s+1)!)
    #                  =  (l-1)!/(s!(l-1-s)!)  * s / (l-s)
    #
    #nth_comb results:
    #  000111
    #  001011
    #  010011
    #  100011
    #  001101
    #  010101
    #  100101
    #  011001
    #  101001
    #  110001
    #  001110
    #  010110
    #  100110
    #  011010
    #  101010
    #  110010
    #  011100
    #  101100
    #  110100
    #  111000








def calc_comb_low_overhead(i,j):
    # 6 choose 2 = 6 * 5 / (1 * 2)
    #  a = 2
    #  b = 4
    a = min(j,i-j)
    b = min(j,i-j)
    r = ov = 1
    pow2 = 0
    pow3 = 0
    for n in range(1,a+1):
        num,denom = i+1-n,n
        #g = math.gcd(num,denom)
        #num //= g
        #denom //= g
        while num & 1 == 0:
            num >>= 1
            pow2 += 1
        while denom & 1 == 0:
            denom >>= 1
            pow2 -= 1
        while num % 3 == 0:
            num //= 3
            pow3 += 1
        while denom % 3 == 0:
            denom //= 3
            pow3 -= 1                        
        r *= num
        ov = max(ov,r)
        assert r%denom == 0, f"division issue with {r} and {denom}"
        r //= denom
    return (r<<pow2)*(3**pow3), ov
    

#check for overflow conditions for c impl
def nth_comb_c(n, length, setbits, length_choose_setbits):
    r = 0
    while length > 1 and setbits > 0:
        g = math.gcd(length,setbits)
        if length%g: print("length/g",length,g)
        if length_choose_setbits%(length//g): print("length_choose_setbits/(length/g)",length_choose_setbits,length//g)
        z = (length_choose_setbits//(length//g)) * (setbits//g)
        if z >>64 : print("z ov",z,length,setbits)
        length -= 1
        if n >= z:
            n -= z
            length_choose_setbits -= z
        else:
            length_choose_setbits = z
            setbits -= 1
            r |= 1<<length
    if length == 1 and setbits>0: r |= 1
    return r
        
