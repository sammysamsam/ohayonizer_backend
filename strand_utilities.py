import random
class strand_utilities:


    
    full_restricted_seq = ['AAAA', 'TTTT', 'CCC', 'GGG']# + palindromic septameric units + palindromic pentameric units + alternating pyurine + pyrimidine

    def __init__(self):
        self.gen_restricted_sequences()

    def gen_restricted_sequences(self):
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
                                self.full_restricted_seq.append(pur1 + pyr1 + pur2 + pyr2 + pur3 + pyr3)
                                self.full_restricted_seq.append(pyr1 + pur1 + pyr2 + pur2 + pyr3 + pur3)


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

    def get_new_restriction_score(e_, seq, new_base, complement_desired):

        score = 0
        new_seq = seq+new_base
        seq_len = len(new_seq)

        #restricted seq check
        for res in e_.full_restricted_seq:
            res_length = len(res)
            if(seq_len >= res_length):           
                if new_seq[seq_len-res_length:] == res:
                    score += 1

        #palindrome check
        if(seq_len > 4):
            six = new_seq[seq_len-4:]
            if e_.palindrome_check(six):
                score += 1
                if six in seq:
                    score += 1

        return score



###################


    def get_pentameric_possiblilities(self, blueprint_pent, complement_desired):
        init_score = strand_utilities.get_new_restriction_score(self, blueprint_pent[0:2] , blueprint_pent[2],complement_desired)
        init_score += strand_utilities.get_new_restriction_score(self, blueprint_pent[0:3] , blueprint_pent[3],complement_desired)
        init_score += strand_utilities.get_new_restriction_score(self, blueprint_pent[0:4] , blueprint_pent[4],complement_desired)


        results = []
        palim_results = []
        if(len(blueprint_pent) == 5):
            lets = []

            for i in range(0,5):
                if(blueprint_pent[i] == 'o'):
                    lets.append("ATCG")
                else:
                    lets.append(blueprint_pent[i])
            
            for i1 in lets[0]:
                for i2 in lets[1]:
                    for i3 in lets[2]:
                        for i4 in lets[3]:
                            for i5 in lets[4]:
                                
                                new_five = i1 + i2 + i3 + i4 + i5
                                new_score = strand_utilities.get_new_restriction_score(self, i1+i2, i3,complement_desired)
                                new_score += strand_utilities.get_new_restriction_score(self, i1+i2+i3, i4,complement_desired)
                                new_score += strand_utilities.get_new_restriction_score(self, i1+i2+i3+i4, i5,complement_desired)
                               
                                if (new_score == init_score):
                                    if (new_five == new_five[::-1]):
                                        palim_results.append(new_five)
                                    else:
                                        results.append(new_five)
        if(len(results) != 0):
            return results
        else:
            return palim_results




#TEST AREA 

"""
util = strand_utilities()

print(util.palindrome_check("TGGCCA"))
print(util.palindrome_check("AAAT"))
print(util.palindrome_check("AAA"))
"""
