a method for correcting errors in golay code using covers
(came up with this several years ago (around 2020?) but forgot to put it online until now)


background theory:
The golay code encodes 12 data bits in 23 message bits in such a way that every pair of valid messages has a hamming distance of at least 7.
This allows either all possible errors of up to 7 message bits to be detected or all possible errors of up to 3 message bits to be corrected.
(in this context an error is when a 0 is received as a 1 or a 1 as a 0, (binary error channel, not binary erasure channel))

The extended golay code is golay with a parity bit appended. This lets it additionally detect (but not correct) 4 bit errors when used as forward-error-correction, or detect up to 8 bit errors when used for error detection.

extended golay code can be generated with this matrix:
                         1 1 1 1 1 1 1 1 1 1 2 2 2 2
     0_1_2_3_4_5_6_7_8_9_0_1_2_3_4_5_6_7_8_9_0_1_2_3
 0 :[1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 1 1 1 1 0 0 0 1]
 1 :[0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 1 1 1 1 0 1 0]
 2 :[0 0 1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 1 1 1 1 1 0 1]
 3 :[0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 1 0 0 1 1 1 1 1 0]
 4 :[0 0 0 0 1 0 0 0 0 0 0 0 1 1 0 0 1 0 0 1 1 1 0 1]
 5 :[0 0 0 0 0 1 0 0 0 0 0 0 1 1 1 0 0 1 0 0 1 1 1 0]
 6 :[0 0 0 0 0 0 1 0 0 0 0 0 1 1 1 1 0 0 1 0 0 1 0 1]
 7 :[0 0 0 0 0 0 0 1 0 0 0 0 1 1 1 1 1 0 0 1 0 0 1 0]
 8 :[0 0 0 0 0 0 0 0 1 0 0 0 0 1 1 1 1 1 0 0 1 0 0 1]
 9 :[0 0 0 0 0 0 0 0 0 1 0 0 0 0 1 1 1 1 1 0 0 1 1 0]
10 :[0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 1 0 1 0 1 0 1 1 1]
11 :[0 0 0 0 0 0 0 0 0 0 0 1 1 0 1 0 1 0 1 0 1 0 1 1]
multiply your 12 bit data vector with it in GF(2) to get your 24 bit message
example: data = 0x123 = [0 0 0 1 0 0 1 0 0 0 1 1]
 data * M  = M[3] + M[6] + M[10] + M[11]
           = [0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 1 0 0 1 1 1 1 1 0]
            +[0 0 0 0 0 0 1 0 0 0 0 0 1 1 1 1 0 0 1 0 0 1 0 1]
            +[0 0 0 0 0 0 0 0 0 0 1 0 0 1 0 1 0 1 0 1 0 1 1 1]
            +[0 0 0 0 0 0 0 0 0 0 0 1 1 0 1 0 1 0 1 0 1 0 1 1]
           = [0 0 0 1 0 0 1 0 0 0 1 1 3 2 2 3 1 1 3 2 2 3 3 3]
           = [0 0 0 1 0 0 1 0 0 0 1 1 1 0 0 1 1 1 1 0 0 1 1 1] in GF(2)
           = 0x1239E7

a message has 0 errors when M * message^t = 0
example: 0x1239E7
 M * m^t = [0 + 0 + 0 + 0 + 1 + 1 + 1 + 1 + 1 + 0 + 0 + 1]   [6]   [0]
           [0 + 0 + 0 + 0 + 0 + 0 + 1 + 1 + 1 + 0 + 1 + 0]   [4]   [0]
           [0 + 0 + 0 + 0 + 0 + 0 + 0 + 1 + 1 + 1 + 0 + 1]   [4]   [0]
           [1 + 0 + 0 + 0 + 1 + 1 + 0 + 0 + 1 + 1 + 1 + 0]   [6]   [0]
           [0 + 0 + 0 + 0 + 1 + 0 + 1 + 0 + 0 + 1 + 0 + 1]   [4]   [0]
           [0 + 0 + 0 + 0 + 1 + 0 + 0 + 1 + 0 + 1 + 1 + 0] = [4] = [0]
           [0 + 1 + 0 + 0 + 1 + 1 + 0 + 0 + 1 + 1 + 0 + 1]   [6]   [0]
           [0 + 0 + 0 + 0 + 1 + 1 + 1 + 0 + 0 + 0 + 1 + 0]   [4]   [0]
           [0 + 0 + 0 + 0 + 0 + 1 + 1 + 1 + 0 + 0 + 0 + 1]   [4]   [0]
           [0 + 0 + 0 + 0 + 0 + 1 + 1 + 1 + 1 + 1 + 1 + 0]   [6]   [0]
           [0 + 0 + 1 + 0 + 0 + 1 + 0 + 1 + 0 + 1 + 1 + 1]   [6]   [0]
           [0 + 0 + 0 + 1 + 1 + 0 + 1 + 0 + 1 + 0 + 1 + 1]   [6]   [0]
         = 0
however: 0x1239E6 has errors, so
 M * m^t = [0 + 0 + 0 + 0 + 1 + 1 + 1 + 1 + 1 + 0 + 0]   [5]   [1]
           [0 + 0 + 0 + 0 + 0 + 0 + 1 + 1 + 1 + 0 + 1]   [4]   [0]
           [0 + 0 + 0 + 0 + 0 + 0 + 0 + 1 + 1 + 1 + 0]   [3]   [1]
           [1 + 0 + 0 + 0 + 1 + 1 + 0 + 0 + 1 + 1 + 1]   [6]   [0]
           [0 + 0 + 0 + 0 + 1 + 0 + 1 + 0 + 0 + 1 + 0]   [3]   [1]
           [0 + 0 + 0 + 0 + 1 + 0 + 0 + 1 + 0 + 1 + 1] = [4] = [0]
           [0 + 1 + 0 + 0 + 1 + 1 + 0 + 0 + 1 + 1 + 0]   [5]   [1]
           [0 + 0 + 0 + 0 + 1 + 1 + 1 + 0 + 0 + 0 + 1]   [4]   [0]
           [0 + 0 + 0 + 0 + 0 + 1 + 1 + 1 + 0 + 0 + 0]   [3]   [1]
           [0 + 0 + 0 + 0 + 0 + 1 + 1 + 1 + 1 + 1 + 1]   [6]   [0]
           [0 + 0 + 1 + 0 + 0 + 1 + 0 + 1 + 0 + 1 + 1]   [5]   [1]
           [0 + 0 + 0 + 1 + 1 + 0 + 1 + 0 + 1 + 0 + 1]   [5]   [1]
         is not 0




working theory:
example: 0x1239E6
 if you assume the first 12 bits are correct and compute the checksum (0x9E7), you find it differes from the message by 3 or fewer bits, this means those differences are the errors.
 however, if you assume the last 12 bits are correct and compute the first part of the message from it, it will differ by more than 3 bits, this means an error is in the bits that were assumed correct.

So, the core principle of this technique is to choose sets of 12 bits to assume as correct, then compute the other bits from those 12, and check the hamming distance to the message.
 If it is 3 or less, the errors are found.
 So, what is needed is the smallest set of choices of 12 bits such that for any combination of 3 bits, there is a choice in the set which chooses none of them.
 So far I've found a set of 27 choices, but smaller sets may exist.
 





