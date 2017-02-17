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
                    restricted_sequences.append(pur1 + pyr1 + pur2 + pyr2)
                    restricted_sequences.append(pyr1 + pur1 + pyr2 + pur2)

# generates all five pentameric units
def genfives():
    global fives
    global fives_present
    lets = 'ACGT'
    for i1 in lets:
        for i2 in lets:
            for i3 in lets:
                for i4 in lets:
                    for i5 in lets:
                        new_five = i1 + i2 + i3 + i4 + i5
                        restricted = False
                        fives.append(new_five)
                        for sequence in restricted_sequences:
                            if sequence in new_five:
                                restricted = True
                                fives_present.append(True)
                                #print(new_five + " is restricted.")
                                break
                        if not restricted:
                            fives_present.append(False)
    print("Length of fives: " + str(len(fives)) + "\nLength of fives_present: " + str(len(fives_present)))
    print("Number of restricted sequences in all possible pentamers: " + str(sum(fives_present)))

# generates all seven nucleotide units in order
def gensevens():
    global sevens
    global sevens_present
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
                                sevens.append(new_seven)
                                for sequence in restricted_sequences:
                                    if sequence in new_seven:
                                        restricted = True
                                        sevens_present.append(True)
                                        #print(new_seven + " is restricted.")
                                        break
                                if not restricted:
                                    sevens_present.append(False)
    print("Length of sevens: " + str(len(sevens)) + "\nLength of sevens_present: " + str(len(sevens_present)))
    print("Number of restricted sequences in all possible pentamers: " + str(sum(sevens_present)))

# sees whether the string of length seven is an approximate match
# to an already present string in the sevens
def sevensapprox(myseven):
    lets = "ACGT"
    j = sevens.index(myseven)
    if sevens_present[j]:
        return True
    i = 1
    while i < 6:
        for let in lets:
            new = myseven[:i] + let + myseven[i+1:]
            j = sevens.index(new)
            if sevens_present[j]:
                print "sevens violation of ", myseven, " with ", new
                return True
        i+=1
    return False

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
    

# generates a new string of size n that doesn't intersect with existing strings
# those other strings are encoded in the fives

def gen_string(strand_length):
    global all_strings
    attempt = 1
    new_strand = []

    while(attempt < 6) and (0 == len(new_strand)):
        new_strand = []
        #for i in range(strand_length): # this loop doesn't seem to be necessary; only the body does
              
        #find pentameric unit that hasnt been added
        current_index = random.randint(0, len(fives) - 1) # we don't need to generate a random int every time, probably
                                                    # just randomize list beforehand
        current_index = 0
        while (current_index < len(fives)) and (fives_present[current_index] or sevens_present[current_index]):
            current_index += 1
            # why does sevens_present need to be there? the indices wouldn't correspond to each other

        #Case: Found pentameric unit
        if (current_index < len(fives)): 
            new_strand = fives[current_index]
       
            curr_length = 5
            while (curr_length < strand_length):
                s = new_strand[len(new_strand) - 4:]     #get previous four bases for( _ _ _ _ + new base )
                if (curr_length >= 7):
                    t = new_strand[len(new_strand) - 6:]   #get previous six bases for( _ _ _ _ _ + new base )
                else:
                    t = ''
            
                good_ones = []
                for base in 'ACGT':
                    pentameric_unit = s+base        # test pentameric unit
                    septameric_unit = t+base        # test septameric unit
            
                    m = fives.index(pentameric_unit)  # find index of pentameric unit
                    if (not fives_present[m]) and ((len(septameric_unit) < 7) or (not sevensapprox(septameric_unit))): #if pentameric unit is not added yet && seven unit does not exist
                        good_ones.append(base)

                if (len(good_ones) > 0):
                    chosen_unit = random.choice(good_ones)
                    new_strand += chosen_unit  #ADD 
                else:
                    curr_length = strand_length # have to stop
                    new_strand = []
                    return "Could not generate a string"
                curr_length += 1

            new_rev_comp_strand = reverse_complement(new_strand)
            update_fives(new_rev_comp_strand)
            update_sevens(new_rev_comp_strand)
            update_fives(new_strand) # added these two lines so that new string would be added to the used fives and sevens
            update_sevens(new_strand)
            all_strings.append(new_strand)
            return new_strand
        else:
            new_strand = []
    return "Could not generate a string in ", attempt, " attempts."
  

# updates the fives data structure with the new string
def update_fives(newstring):
    global fives_present
    for i in range(len(newstring)-4):
        s = newstring[i:i+5]
        j = fives.index(s)
    # if fives_present[j]:
      # print "We have a fives problem at position ", i, " of newstring ", newstring, " with respect to ", fives[j]
      # print "The string in question is: ",  s
        fives_present[j] = True
    

# updates the sevens data structure with the new string
def update_sevens(newstring):
    global sevens_present
    for i in range(len(newstring)-6):
        s = newstring[i:i+7]
        j = sevens.index(s)
        #if sevens_present[j]:
            #print "We have a sevens problem at position ", i, " of newstring"
            #print "Letters are: ",  s
        sevens_present[j] = True


# DATA
fives = []
fives_present = []
sevens = []
sevens_present = []
all_strings = []
restricted_sequences = []

n = 7 # size of strings

# EXECUTION
gen_restricted_sequences()
genfives()
num_restricted_fives = sum(fives_present)
gensevens()
num_restricted_sevens = sum(sevens_present)
for i in range(0, 45):
    print gen_string(n)

x = sum(fives_present) - num_restricted_fives
print "number of fives that are present (including complements): ", x
x = sum(sevens_present) - num_restricted_sevens
print "number of sevens that are present (including complements): ", x
