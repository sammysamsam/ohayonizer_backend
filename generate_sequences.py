# This software designs strings that don't have reverse complements
# of length five with other strings and that don't have groups of 7 
# nucleotides with 6 in common.

# For now, it makes no restriction about strings being complementary
# to themselves nor does it forbid two strings from being
# very similar. Please tell me what those constraints should be.



"""







"""

import random

def gen_restricted_sequences():
    global restricted_sequences
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

# generates all five pentameric units
def genfives():
    global fives
    lets = 'ACGT'
    for i1 in lets:
        for i2 in lets:
            for i3 in lets:
                for i4 in lets:
                    for i5 in lets:
                        new_five = i1 + i2 + i3 + i4 + i5
                        restricted = False
                        for sequence in ['AAAA', 'TTTT', 'CCC', 'GGG'] :
                            if sequence in new_five:
                                restricted = True
                                #print(new_five + " is restricted.")
                                break
                        if not restricted:
                            fives[new_five] = False

    print("Length of fives: " + str(len(fives)) )

# generates all seven nucleotide units in order
def gensevens():
    global sevens
    lets = 'ACGT'
    for i1 in lets:
        for i2 in lets:
            for i3 in lets:
                for i4 in lets:
                    for i5 in lets:
                        for i6 in lets:
                            for i7 in lets:
                                new_seven = i1 + i2 + i3 + i4 + i5 + i6 + i7
                                restricted = False
                                for sequence in restricted_sequences:
                                    if sequence in new_seven:
                                        restricted = True
                                        #print(new_seven + " is restricted.")
                                        break
                                if not restricted:
                                    sevens[new_seven] = False
    print("Length of sevens: " + str(len(sevens)) + "\n")

# take reverse complement of DNA
def reverse_complement(s):
    srev = s[::-1]
    out = ''
    for let in srev:
        if (let == 'A'):
            out+= 'T'
        elif (let == 'C'):
            out+= 'G'
        elif (let == 'G'):
            out+= 'C'
        elif (let == 'T'):
            out+= 'A'
    return out




def get_next_base(prev_4,prev_6,complement_exists):
    next_possible_bases = []

    for base in 'ACGT':
        pentameric_unit = prev_4+base        # test pentameric unit
        septameric_unit = prev_6+base        # test septameric unit

        p_ = fives.get(pentameric_unit,'none')  # pentameric unit exists
        s_ = sevens.get(septameric_unit,'none') #septameric unit exists
 
        if (not p_) and ((not s_) or (len(septameric_unit) != 7 )): #if pentameric unit is not added yet && seven unit does not exist
            next_possible_bases.append(base)
    return next_possible_bases


# generates a new string of size n that doesn't intersect with existing strings
# those other strings are encoded in the fives

def gen_string(strand_length, complement_exists):
    global all_strings
    global fives
    global sevens

    #used for revert fives and sevens dictionary after failed attempts
    saved_fives = fives.copy()
    saved_sevens = sevens.copy()

    attempt = 1
    new_strand = []

    #ensure enough unadded units are available 
    while(attempt < 6) and fives.values().count(False) > strand_length:
        new_strand = []

        #randomize dictionary??


        #get first pentameric unit and add unit         
        for seq,exist in fives.iteritems():
            if(not exist):
                new_pentameric_unit = seq
                break
        fives[new_pentameric_unit] = True
        curr_length = 5

        #build rest of strand
        new_strand = new_pentameric_unit
        while (curr_length < strand_length):
            
            #get previous four bases for( _ _ _ _ + new base )
            prev_4 = new_strand[len(new_strand) - 4:]   

            #get previous six bases for( _ _ _ _ _ + new base )       
            if (curr_length >= 7):
                prev_6 = new_strand[len(new_strand) - 6:] 
            else:
                prev_6 = ''

            #get next possible base
            next_possible_bases = get_next_base(prev_4,prev_6,complement_exists)


            #CASE: add possible base (CONTINUE)
            if (len(next_possible_bases) > 0):
                chosen_unit = random.choice(next_possible_bases)
                new_strand += chosen_unit  
                update_fives(prev_4+chosen_unit, complement_exists)
                update_sevens(prev_6+chosen_unit, complement_exists)

            #CASE: no possible base (STOP)
            else:
                curr_length = strand_length 
                new_strand = []
            curr_length += 1
        
        
        if len(new_strand) == strand_length:
            all_strings.append(new_strand)
            return new_strand
        else:
            print("attempt failed: "+ str(attempt))     
            
            #undo changes to dictionary during failed attempt
            fives = saved_fives.copy()
            sevens = saved_sevens.copy()
            attempt +=1

    return "Could not generate a string in ", attempt, " attempts."
  



# updates the fives data structure with the new string
def update_fives(new_unit,complement_exists):
    if(len(new_unit) == 5):
        if(complement_exists):
            new_reverse_unit = reverse_complement(new_unit)
            fives[new_reverse_unit] = True
        fives[new_unit] = True
    
# updates the sevens data structure with the new string
def update_sevens(new_unit,complement_exists):
    if(len(new_unit) == 7):
        if(complement_exists):
            new_reverse_unit = reverse_complement(new_unit)
            update_sevens(new_reverse_unit,False)
    sevens[new_unit] = True
   
    #sevens approx added
    lets = "ACGT"
    i = 1
    while i < 6:
        for let in lets:
            new = new_unit[:i] + let + new_unit[i+1:]

            # if unit exists and has not been used yet
            if( sevens.get(new,'none') != 'none' and (not sevens[new])):
                sevens[new] = True
        i+=1




#MAIN METHOD

fives = {}
sevens = {}
all_strings = []
restricted_sequences = []

size_of_strand = 50

gen_restricted_sequences()
genfives()
gensevens()

for i in range(0, 10):
    print gen_string(size_of_strand,False)
    print("pentameric units remaining:" + str(fives.values().count(False)))
    print("septameric units remaining:" + str(sevens.values().count(False)))   
    print("\n")

