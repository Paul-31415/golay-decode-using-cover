golay_mat = [0b100000000000100111110001,
             0b010000000000010011111010,
             0b001000000000001001111101,
             0b000100000000100100111110,
             0b000010000000110010011101,
             0b000001000000111001001110,
             0b000000100000111100100101,
             0b000000010000111110010010,
             0b000000001000011111001001,
             0b000000000100001111100110,
             0b000000000010010101010111,
             0b000000000001101010101011]
#⠑⢄⠀⠀⠀⠀⡑⢌⠻⣿⣦⡪⠀
#⠀⠀⠑⢄⠀⠀⣿⣦⡑⢌⠻⡪
#⠀⠀⠀⠀⠑⢄⡨⡻⡻⡢⡱⣮

golay_mat = [
    0b100000000000101000111011,
    0b010000000000110100011101,
    0b001000000000011010001111,
    0b000100000000101101000111,
    0b000010000000110110100011,
    0b000001000000111011010001,
    0b000000100000011101101001,
    0b000000010000001110110101,
    0b000000001000000111011011,
    0b000000000100100011101101,
    0b000000000010010001110111,
    0b000000000001111111111110,]
def transpose_bitmat(m,w=23):
    r = [0]*w
    for i in range(w):
        for j in range(len(m)):
            r[i] |= ((m[j]>>i)&1)<<j
    return r
    
golay_pcheck_mat =  transpose_bitmat([
    0b100111000111,
    0b101011011001,
    0b101101101010,
    0b101110110100,
    0b110011101100,
    0b110101110001,
    0b110110011010,
    0b111001010110,
    0b111010100011,
    0b111100001101,
    0b011111111111],12)
golay_pcheck_rmat = transpose_bitmat([1<<i for i in range(10,-1,-1)]+golay_pcheck_mat,11)
def gf2vm(v,m):
    r = 0
    for i in range(len(m)):
        r ^= ((v>>i)&1)*m[len(m)-1-i]
    return r

def errs(b=23,me=3):
    from combinations import enumerate_combs
    for i in range(me+1):
        yield from enumerate_combs(b,i)

def correct_via_pcheck(w):
    def parity(w):
        return gf2vm(w,golay_pcheck_mat)^(w>>12)
    
    for c in errs():
        if parity(w^c) == 0:
            return w^c
import functools
@functools.cache
def can_use_as_parity(w):
    mat = list(golay_pcheck_rmat)
    #do gauss on
    #choosing rows (bitposes) ~w to be identity
    #eliminating on cols
    r = 0
    fail = 0
    for b in range(23):
        if (w>>b)&1:
            for i in range(r,11):
                if (mat[i]>>b)&1:
                    mat[r],mat[i] = mat[i],mat[r]
                    break
            else:
                fail |= (~mat[r]>>b)&1
            for i in range(11):
                if i != r and (mat[i]>>b)&1:
                    mat[i] ^= mat[r]
            r += 1
    return w,mat,fail
    
    
    
        

# find >=4 sets of 11 bits to use as parity so that
# any combination of 3 bit errors must all occur together in at least 1 group
# . . . . . . . . . . . . . . . . . . . . . . .
# 1 1 1 1 1 1 1 1 1 1 1
#                         2 2 2 2 2 2 2 2 2 2 2
#   3   3   3   3   3   3   3   3   3   3   3
#                    
#                     E E E
#this means that for any two errors, all remaining error locations must be covered
#which means from 23 choose 3 (= 1771)
#   group them together in to sets which or to 11 set bits
#    and cover all of them using the fewest sets.

def find_golay_sets():
    from combinations import enumerate_combs
    from bits import popcount
    words = [0]
    missed = []
    def adde(c,words=words,missed=missed):
        for w in words:
            if w&c==c:
                break
        else:
            p = popcount(words[-1]|c)
            if p <= 11 and not can_use_as_parity(words[-1]|c)[2]:
                words[-1] |= c
                if p == 11:
                    words.append(0)
            else:
                missed.append(c)
    import random
    r = list(enumerate_combs(23,3))
    def do(giveup=30):
        for i in r:
            adde(i)
            if len(words)>giveup:
                return giveup+1
        while len(missed):
            pl = 0
            while pl != len(missed):
                pl = len(missed)
                m = list(missed)
                random.shuffle(m)
                missed.clear()
                for i in m:
                    adde(i)
            if len(missed):
                words.append(0)
            #print(len(missed),end="   \r")
        return len(words)
    best = 1000000
    while 1:
        random.shuffle(r)
        b = do(best)
        if b <= best:# and can_do_checkset(words)[0]:
            best = b
            print("------------",b,"-------------")
            print(words)
        words.clear()
        words.append(0)
    return words,missed
            
# 27:
# [2843756, 832739, 6735656, 2694457, 2434026, 2685459, 8130932, 4429575, 4339838, 8325761, 1424436, 5419659, 1462224, 7553430, 6256642, 7262221, 6336340, 551645, 1836847, 1803736, 6406379, 3306913, 350964, 7638170, 1286359, 1764618, 4982229]
import numpy as np

from bits import popcount, unselect_bits


def to_bit_arr(l):
    return np.array([[(v>>i)&1 for i in range(23)] for v in l])

def from_bit_arr(a):
    return [sum(a[r,c]*(1<<c) for c in range(a.shape[1])) for r in range(a.shape[0])]
def pa(a):
    for y in range(a.shape[0]):
        for x in range(a.shape[1]):
            print(" #"[a[y,x]],end="")
        print("",end="\n")
        
    
#the same 27 but permuted rows/cols: [2047, 129055, 1129017, 1472385, 1633106, 1750556, 1911220, 1968239, 2219368, 2242983, 2904731, 3113763, 3222764, 3974283, 4660132, 4957482, 4970689, 5129893, 5560652, 5585987, 5751352, 6439389, 6704722, 6889072, 7398054, 7967122, 8265940]
checksets = [2047, 129055, 1129017, 1472385, 1633106, 1750556, 1911220, 1968239, 2219368, 2242983, 2904731, 3113763, 3222764, 3974283, 4660132, 4957482, 4970689, 5129893, 5560652, 5585987, 5751352, 6439389, 6704722, 6889072, 7398054, 7967122, 8265940]

checksets = [944922, 2726132, 4333229, 88942, 2290386, 5258697, 6670393, 2069377, 2397085, 5803349, 3519558, 7987984, 2876430, 3190578, 6047887, 3961881, 4628963, 5448428, 6952531, 3826145, 4685088, 7620322, 860797, 4932790, 3378377, 7733654, 5560929, 77575]

checksets = [944718, 7868989, 2836369, 8357024, 4508951, 5412451, 4656937, 4552803, 5465204, 7916000, 42493, 3384970, 2206163, 3513093, 2618392, 5263791, 1644766, 1991047, 3068500, 2637478, 1334244, 5026712, 1719745, 2001714, 7055875, 8193864, 6605322]


def decode_golay_l(checksets = checksets):
    from bits import popcount
    from bits import select_bits
    mats = []
    if 0:
      for cs in checksets:
        i = 0
        while popcount(cs) < 11:
            cs |= 1<<i
            i += 1
            
        mat = [(i>>1) for i in golay_mat]
        #run gaussian elimination on the cols for which cs is 0
        r = 0
        for i in range(22,-1,-1):
            if (cs>>i)&1:
                continue
            for j in range(len(mat)-1,r-1,-1):
                if (mat[j]>>i)&1:
                    mat[j],mat[r] = mat[r],mat[j]
                    break
            else:
                if (mat[r]>>i)&1 == 0:
                    print("noninvertible matrix:",cs,i,len(mats))
            for j in range(r+1,len(mat)):
                if j != r and (mat[j]>>i)&1:
                    mat[j] ^= mat[r]
            r += 1
        mats.append((cs,mat[::-1]))
    mats = [can_use_as_parity(cs) for cs in checksets]
    if any(v[2] for v in mats):
        print("noninvertable matrix")
    mats = [(a,transpose_bitmat(b),c) for a,b,c in mats]
    def decode_golay(word,d=0,mats=mats):
        if d:
            print("word:  ",bin((1<<23)|word)[3:]);j = 0
        for (cs,mat,f) in mats:
            if d:
                print("trying ",bin((1<<23)|cs)[3:].replace("0"," ").replace("1","C"),j,f);j+=1
            dat = cs^((1<<23)-1)
            #kw = select_bits(word,dat)
            #c = word
            #for i in range(12):
            #    if (kw>>i)&1:
            #i = -1
            #for b in range(23):
            #    if (cs>>b)&1:
            #        i += 1
            #        if (c>>b)&1:
            #            c ^= mat[i]
            c = gf2vm(word,mat)
            r = 0
            for b in range(23):
                if mat[b]&c==mat[b]:
                    r |= 1<<(22-b)
            if d:
                print("got    ",bin((1<<23)|r)[3:])
            if popcount(c) <= 3:
                return r^word
                        
    return decode_golay,mats

decode_golay,dgm = decode_golay_l()
def pdm(v):
    print(bin(v[0]|(1<<23))[3:]);pa(to_bit_arr(v[1])[:,::-1])


def full_code(cache=[]):
    if len(cache):
        return cache[0]
    r = []
    for i in range(1<<12):
        c = 0
        for b in range(12):
            if (i>>b)&1:
                c ^= golay_mat[b]
        r.append(c>>1)
    cache.append(r)
    return r
def full_code_slow(cache = []):
    if len(cache):
        return cache[0]
    from combinations import enumerate_combs
    c = {0}|set(enumerate_combs(23,1))|set(enumerate_combs(23,2))|set(enumerate_combs(23,3))|set(enumerate_combs(23,4))|set(enumerate_combs(23,5))|set(enumerate_combs(23,6))#|set(enumerate_combs(23,7))
    words = []
    for i in range(1<<23):
        for w in words:
            if w^i in c:
                break
        else:
            print(i,end="\r")
            words.append(i)
    cache.append(words)
    return words
        
    
def ref_decode(v,good=[]):
    if len(good) == 0:
        from combinations import enumerate_combs
        good.append({0}|set(enumerate_combs(23,1))|set(enumerate_combs(23,2))|set(enumerate_combs(23,3)))
    for w in full_code():
        if w^v in good[0]:
            return w
def ref_mat(w):
    bits = set()
    w ^= (1<<23)-1
    for i in range(23):
        if (1<<i)&w:
            bits.add(1<<i)
    print([bin(i) for i in bits])
    a = [v for v in full_code() if v&w in bits]
    a.sort(key=lambda v: v&w)
    return (w ^ ((1<<23)-1),a)
    

#another 27 solution
# [5097868, 6136432, 7365167, 385846, 6605030, 6724681, 5655189, 2344625, 1233890, 713883, 2299132, 1322603, 4750564, 8209688, 2446169, 1938726, 4406595, 2855213, 4202327, 3348890, 3979793, 4073346, 928118, 1811148, 1520073, 8069473, 5855778]


def prove():
    for i in range(0,1<<23,1<<12):
        for j in range(1<<12):
            if decode_golay(i+j) is None:
                print(i)
                return i
        print(i,end="\r")
        

    
def check_sol(words):
    from combinations import enumerate_combs
    for c in enumerate_combs(23,3):
        for w in words:
            if w&c==c:
                break
        else:
            return c
    
def can_do_checkset(checkset):
    n = 0
    for cs in checksets:
        i = 0
        while popcount(cs) < 11:
            cs |= 1<<i
            i += 1
        if can_use_as_parity(cs)[2]:
            return False,n
        n += 1
    return True,n

    
#best found so far
#------------ 27 -------------
#[7116081, 351216, 406253, 1719461, 7506822, 5000395, 5780366, 1557589, 273790, 4313501, 3649888, 4418323, 6113104, 3191994, 4020235, 2223668, 8151209, 6389320, 5882086, 5937768, 7538089, 5680645, 1019778, 2916975, 2773598, 2843089, 1763378]
#------------ 27 -------------
#[6972872, 4663501, 1941813, 4064978, 2921705, 4906005, 3232692, 2914823, 7594514, 315275, 243107, 7474477, 4275050, 1479262, 333782, 5776295, 1863771, 7675712, 863804, 7018673, 6865532, 478065, 1682762, 7169066, 5540344, 3857666, 43190]
#------------ 27 -------------
#[89718, 5123598, 4903251, 1073037, 2656233, 2358464, 467771, 8204433, 3524902, 5483045, 5642694, 3887690, 1880436, 2542172, 6645089, 4812986, 1455848, 6382238, 5096597, 1238682, 2758839, 8281344, 5343609, 3212607, 633810, 6542176, 433229]
#------------ 27 -------------
#[5418188, 390682, 4825715, 3350416, 1689830, 6510712, 3056790, 2230007, 3205466, 963940, 1330493, 7889453, 7094770, 4531149, 1874585, 6045220, 6343823, 2759463, 5322007, 4696355, 4145419, 2749513, 3562051, 2352932, 4354478, 5424040, 4886720]
#------------ 27 -------------
#[944718, 7868989, 2836369, 8357024, 4508951, 5412451, 4656937, 4552803, 5465204, 7916000, 42493, 3384970, 2206163, 3513093, 2618392, 5263791, 1644766, 1991047, 3068500, 2637478, 1334244, 5026712, 1719745, 2001714, 7055875, 8193864, 6605322]
