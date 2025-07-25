
def bit_reverse_pow2(w,b):
    if b <= 0:
        return w
    b -= 1
    m = ((1<<(1<<b))-1)
    while b > 0:
        w = ((w&m)<<(1<<b)) | ((w&~m)>>(1<<b))
        b -= 1
        m ^= m<<(1<<b)
    w = ((w&m)<<1) | ((w&~m)>>1)
    return w

def bit_length(w):
    w = abs(w)
    b = 1
    while w > 1<<b:
        b <<= 1
    b >>= 1
    c = b
    while c > 0:
        if w > 1<<b<<c:
            b |= c
        c >>= 1
    return b + (w >> b)

def bit_reverse(w,bits):
    if bits <= 1:
        return w
    l = bit_length(bits-1)
    return bit_reverse_pow2(w,l)>>((1<<l)-bits)

def popcount_pow2(w,b):
    for c in range(b):
        m = (((1<<(1<<c))-1)<<(1<<b))//((1<<(2<<c))-1)
        w = (w&m)+((w>>(1<<c))&m)
    return w

def popcount(w):
    if w < 0:
        return -popcount(-w)
    l = bit_length(w)
    b = bit_length(l-1)
    return popcount_pow2(w,b)

def expand_bits(w,sep=2):
    if sep == 1:
        return w
    b = bit_length(bit_length(w)-1)
    if b <= 0:
        return w
    b -= 1
    m = ((1<<(1<<b))-1)
    while b > 0:
        w = (w&m)|((w&~m)<<((sep-1)<<b))
        b -= 1
        m &= m>>(1<<b)
        m |= m<<(sep<<(b+1))
    w = (w&m)|((w&~m)<<(sep-1))
    return w

def contract_bits(w,sep=2):
    if sep == 1:
        return w
    m = (1<<(bit_length(w)//sep*sep+sep))//((1<<sep)-1)
    w &= m
    b = bit_length(bit_length(w)-1)
    for c in range(b):
        patbits = (sep<<c)<<1
        m = (((1<<(1<<c))-1)<<((1+(1<<b)//patbits)*patbits))//((1<<patbits)-1)
        if m == 0:
            break
        w = (w&m)|(((w>>(sep<<c))&m)<<(1+c))
    return w


def select_bits(w,mask):
    r = 0
    b = 0
    for i in range(bit_length(mask)):
        if (mask>>i)&1:
            r |= (w>>(i-b))&(1<<b)
            b += 1
    return r
def unselect_bits(w,mask):
    r = 0
    b = 0
    for i in range(bit_length(mask)):
        if (mask>>i)&1:
            r |= (w<<(i-b))&(1<<i)
            b += 1
    return r



def code(bits,dist,poplut=None):
    poplut = [popcount(i) for i in range(1<<bits)] if poplut is None else poplut
    c = []
    for i in range(1<<bits):
        for v in c:
            if (poplut[v^i] if poplut != 0 else popcount(v^i))<dist:
                break
        else:
            c.append(i)
    return c


        
def unique_dists(tot=2,length=4):
    if tot > 0:
        if length == 1:
            yield (tot,)
        else:
            for i in range(tot+1):
                for v in unique_dists(tot-i,length-1):
                    yield (i,)+v
    else:
        yield (0,)*length

