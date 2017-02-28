"""

This software designs strand sequences, including its Zx[ bvcx complements, with unique pentameric (5) and 
septameric (7) units, i.e. AATCG & TAGGTCC.


in addition:
    - no polypurine (AAAA/TTTT) and  + CCC/GGG
    - avoid alternating CGs and other purines and pyrmadines (3-4)
    - no palindromic units
"""


import random
from strand_utilities import strand_utilities
util = strand_utilities()


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

                        # NO PALINDROMES, POLYPURINES    
                        if (not restricted) and (not new_five == new_five[::-1]) :
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
                                
                                # NO PALINDROMES, ALTERNATING PURINE-PYRIMIDINE (vice versa)
                                if (not restricted) and (not new_seven== new_seven[::-1]):
                                    sevens[new_seven] = False
    print("Length of sevens: " + str(len(sevens)) + "\n")


#####################


def get_next_base(prev_4,prev_6,blueprint,index):
    five_restricted = ['AAAA', 'TTTT', 'CCC', 'GGG','CCCC','GGGG','CCCCC','GGGGG','TTTTT','AAAAA']    
    global restricted_sequences
    next_possible_bases = []

    if(blueprint[index] != 'o'):
        next_possible_bases.append(blueprint[index])
        return next_possible_bases

    #blueprint contains specific bases in current pentameric unit

    if(blueprint[index - 4:index] != "oooo" ):
        init_score = 0
        if(len(prev_4) == 4):
            init_score = util.get_restriction_score(prev_4,five_restricted)
        if(len(prev_6) == 6):
            init_score += util.get_restriction_score(prev_6,restricted_sequences)

        for base in 'ACGT':
           
            pentameric_unit = prev_4+base  # test pentameric unit
            septameric_unit = prev_6+base  # test septameric unit
            new_score = 0

            if(len(prev_4) == 4):
                new_score = util.get_restriction_score(pentameric_unit,five_restricted)
            if(len(prev_6) == 6):
                new_score += util.get_restriction_score(septameric_unit,restricted_sequences)

            if(new_score == init_score):
                next_possible_bases.append(base)

    #blueprint does not contain specific bases in current pentameric unit

    else:
        for base in 'ACGT':
           
            pentameric_unit = prev_4+base           # test pentameric unit
            septameric_unit = prev_6+base           # test septameric unit

            p_ = fives.get(pentameric_unit,'none')  # pentameric unit exists
            s_ = sevens.get(septameric_unit,'none') # septameric unit exists
     
            if (not p_) and ((not s_) or (len(septameric_unit) != 7 )) : #if pentameric unit is not added yet && seven unit does not exist
                next_possible_bases.append(base)
    return next_possible_bases


def get_first_five_bases(blueprint):
    possibilities = util.get_pentameric_possiblilities(blueprint[0:5])

    if(len(possibilities)==0):
        print("no possibilities for first five bases")
        return "none"

    random.shuffle(possibilities)
    return possibilities[0]




def process_blueprint(strand_length, blueprint):
    if(len(blueprint) == 0):
        return blueprint.ljust(strand_length, 'o')
    else:
        return blueprint


#####################


# generates a new string of size n that doesn't intersect with existing strings
# those other strings are encoded in the fives

def gen_string(strand_length,blueprint, complement_exists):
    global all_strings
    global fives
    global sevens

    #used for revert fives and sevens dictionary after failed attempts
    saved_fives = fives.copy()
    saved_sevens = sevens.copy()
    blueprint = process_blueprint(strand_length, blueprint)
    new_strand = []

    #check if algorithm can start with building first five bases
    starting_bases = get_first_five_bases(blueprint)

    #ensure enough unadded units are available 
    attempt = 1
    while(attempt < 2500) and list(fives.values()).count(False) > strand_length and len(starting_bases) == 5:
        
        new_strand = starting_bases
        curr_length = 5
        while (curr_length < strand_length):
            
            #get previous four bases for( _ _ _ _ + new base )
            prev_4 = new_strand[len(new_strand) - 4:]   

            #get previous six bases for( _ _ _ _ _ + new base )       
            prev_6 = ''
            if (curr_length >= 7):
                prev_6 = new_strand[len(new_strand) - 6:] 

            #get next possible base
            next_possible_bases = get_next_base(prev_4,prev_6,blueprint,curr_length)

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
            print("order completed at attempt #: "+ str(attempt))     
            update_all_fives_sevens(new_strand,complement_exists)
            all_strings.append(new_strand)
            return new_strand
        else:
            #print("attempt failed: "+ str(attempt))     

            #undo changes to dictionary during failed attempt
            fives = saved_fives.copy()
            sevens = saved_sevens.copy()
            attempt +=1

    return "Could not generate a string in ", attempt, " attempts."
  



#####################


def update_all_fives_sevens(sequence,complement_exists):
    global fives 
    global sevens
    if(complement_exists):
        for i in range(len(sequence)-4):
            five = sequence[i:i+5]
            if(five in fives):
                fives[five] = True

        for d in range(len(sequence)-6):
            sev = sequence[d:d+7]
            if(sev in sevens):
                sevens[sev] = True




# updates the fives data structure with the new string
def update_fives(new_unit,complement_exists):
    global fives
    if(len(new_unit) == 5):
        if(complement_exists and (new_unit in fives)):
            fives[new_unit] = True

        new_reverse_unit = util.reverse_complement(new_unit)
        if( new_reverse_unit in fives):
            fives[new_reverse_unit] = True


# updates the sevens data structure with the new string
def update_sevens(new_unit,complement_exists):
    global sevens
    if(len(new_unit) == 7):
        if(complement_exists and (new_unit in sevens)):
            sevens[new_unit] = True
            add_sevens_approx(new_unit)

        new_reverse_unit = util.reverse_complement(new_unit)       
        if(new_reverse_unit in sevens):     
            sevens[new_reverse_unit]

        add_sevens_approx(new_reverse_unit)

def add_sevens_approx(new_unit):
        global sevens
        lets = "ACGT"
        i = 1
        while i < 6:
            for let in lets:
                new = new_unit[:i] + let + new_unit[i+1:]
                if( new in sevens and (not sevens[new])):
                    sevens[new] = True
            i+=1



########          MAIN METHOD        ############## 

fives = {}
sevens = {}
all_strings = []
restricted_sequences = []
size_of_strand = 50

gen_restricted_sequences()
genfives()
gensevens()

test_str = "CCoATooooAAAoooooooATTooGGGGGo"
test_str_2 = "oCCCooooooooooooooATTooGGGGGGo"
test = gen_string(30,test_str, True)
test2 = gen_string(30, test_str_2, True)
print(test_str)
print(test)
print(test_str_2)
print(test2)

"""
for i in range(0, 10):
    print gen_string(size_of_strand,"",True)
    print("pentameric units remaining:" + str(fives.values().count(False)))
    print("septameric units remaining:" + str(sevens.values().count(False)))   
"""
