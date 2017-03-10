import random
class strand_utilities:


    five_restricted_seq = ['AAAA', 'TTTT', 'CCC', 'GGG']
    full_restricted_seq = []

    def __init__(self):
        self.gen_full_restricted_sequences()


    def gen_full_restricted_sequences(self):
        restricted_sequences = ['AAAA', 'TTTT', 'CCC', 'GGG'] # smallest prohibited repeats of each category
        pyrimadines = 'TC'
        purines = 'AG'

        # Generate alternating sequences
        for pur1 in purines:
            for pyr1 in pyrimadines:
                for pur2 in purines:
                    for pyr2 in pyrimadines:
                        for pur3 in purines:
                            for pyr3 in pyrimadines:
                                restricted_sequences.append(pur1 + pyr1 + pur2 + pyr2 + pur3 + pyr3)
                                restricted_sequences.append(pyr1 + pur1 + pyr2 + pur2 + pyr3 + pur3)
        #print(restricted_sequences)
        self.full_restricted_seq = restricted_sequences


############################################

    # take reverse complement of DNA
    def reverse_complement(e, s):
        srev = s[::-1]
        return strand_utilities.complement(e, srev)

    def complement(e,s):
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
        return out       

    #get one random base
    def random_base():
        baselist = ['A','T','C','G']
        return random.choice(baselist)


    #get random sequence with given length
    def random_sequence(length):
        seq = ""
        for i in range(1,length):
            seq += random_base();
        return seq

    """
    Based on a given blueprint, it returns an array of possible sequences that dont violate or worsen 
    number of restricted sequences.

    if palindromes are only found, it returns array with palindromes

    """

    def get_new_restriction_score(e, seq, new_base, complement_desired):

        score = 0
        seq = seq+new_base
        seq_len = len(seq)

        #if(complement_desired):
        #    rev_seq = strand_utilities.complement(e,seq)

#        print(seq + " "+ rev_seq)
        for res in e.full_restricted_seq:
                res_length = len(res)

                if(seq_len >= res_length):           
                    if seq[seq_len-res_length:] == res:
                        #print('Restricted sequence found')
                        #print('Substring of current string ' + seq + ' is ' + res)
                        score += 1
        #            if complement_desired and rev_seq[seq_len-res_length:] == res:
        #                score += 1
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



#x = util.get_pentameric_possiblilities("GGGGG")
#print(x)
util = strand_utilities()
print(util.get_new_restriction_score("GACAGG",'G',True))



