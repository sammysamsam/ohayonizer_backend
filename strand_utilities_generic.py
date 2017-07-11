import random

    
full_restricted_seq = ['AAAA', 'TTTT', 'CCC', 'GGG', 'GGGG', 'CCCC']# + palindromic septameric units + palindromic pentameric units + alternating pyurine + pyrimidine


def get_restricted_sequences():
    pyrimadines = 'TC'
    purines = 'AG'

    lets = 'ACGT'
    #alternating purine - pyrimidines
    for pur1 in purines:
        for pyr1 in pyrimadines:
            for pur2 in purines:
                for pyr2 in pyrimadines:
                    for pur3 in purines:
                        for pyr3 in pyrimadines:
                            full_restricted_seq.append(pur1 + pyr1 + pur2 + pyr2 + pur3 + pyr3)
                            full_restricted_seq.append(pyr1 + pur1 + pyr2 + pur2 + pyr3 + pur3)


############################################

# take reverse complement of DNA
def reverse_complement(e_, s):
    rev = s[::-1]
    return e_.complement(rev)

def complement(e_, s):
    out = ''
    for let in s:
        if (let == 'A'):
            out+= 'T'
        elif (let == 'C'):
            out+= 'G'
        elif (let == 'G'):
            out+= 'C'
        elif (let == 'T'):
            out+= 'A'
        else:
            out+=let
    return out

#get one random base
def random_base(e_):
    baselist = ['A','T','C','G']
    return random.choice(baselist)


#get random sequence with given length
def random_sequence(e_, length):
    seq = ""
    for i in range(1,length):
        seq += random_base();
    return seq

def palindrome_check(e_, seq):
    #print("palindrome: " + seq + "  = "+ str((seq == e_.reverse_complement(seq))))
    if "o" in seq:
        return False
    return (seq == e_.reverse_complement(seq))

###################


"""
Based on a given blueprint, it returns an array of possible sequences that dont violate or worsen 
number of restricted sequences.

if palindromes are only found, it returns array with palindromes

"""

def get_new_restriction_score(e_, seq, new_base):

    score = 0
    new_seq = seq + new_base
    seq_len = len(new_seq)

    #restricted seq check
    for res in e_.full_restricted_seq:
        res_length = len(res)
        if(seq_len >= res_length):
            if new_seq[seq_len-res_length:] == res:
                score += 1

    #palindrome check
    if(seq_len >= 4):
        six = new_seq[-4:]
        if e_.palindrome_check(six):
            score += 2
            if six in seq:
                score += 3

    return score



###################



#TEST AREA 


# util = strand_utilities()

# print(util.get_new_restriction_score('AAA',"A"))
# print(util.palindrome_check("AAAT"))
# print(util.palindrome_check("AAA"))
