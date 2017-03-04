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
                        else:
                            fives[new_five] = True
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
                                for sequence in util.full_restricted_seq:
                                    if sequence in new_seven:
                                        restricted = True
                                        break
                                
                                # NO PALINDROMES, ALTERNATING PURINE-PYRIMIDINE (vice versa)
                                if (not restricted) and (not new_seven== new_seven[::-1]):
                                    sevens[new_seven] = False
                                else:
                                    sevens[new_seven] = True
    print("Length of sevens: " + str(len(sevens)) + "\n")


#####################

def calc_five_restriction_score(sequence):
    total = 0;
    if(len(sequence) > 4):
        total += util.get_restriction_score(sequence, util.five_restricted_seq)
    return total

def calc_seven_restriction_score(sequence):
    total = 0;
    if(len(sequence) > 6):
        return util.get_restriction_score(sequence, util.full_restricted_seq)
    return total  


'''
Take the blueprint and keeps track of all the violations up to an index value iterating through the blueprint,
then checks if it increases when you add a base. If it does, use another base instead.

oooooooggooo
blueprint violation array: all 0s (no violations)
during string generation, if we put in:  oooooogggooo
then the number of violations will be greater than in the blueprint violation array
now the blueprint array has:
    ooogggooo
now the blueprint violation array will look like: [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]

'''

def get_blueprint_violation_array(blueprint):
    violation_array = []
    curr_blueprint_seq = ""

    recent_score = 0

    for base in blueprint:
        #print(curr_blueprint_seq)
        curr_blueprint_seq += base
        pent = ""
        sept = ""
        if(len(curr_blueprint_seq) > 4):
            pent = curr_blueprint_seq[len(curr_blueprint_seq) - 5:]
        if(len(curr_blueprint_seq) > 6):
            sept = curr_blueprint_seq[len(curr_blueprint_seq) - 7:]

        score = util.get_restriction_score(sept, util.full_restricted_seq)
        score += util.get_restriction_score(pent, util.five_restricted_seq)

        violation_array.append(score)

    return violation_array




def process_blueprint(strand_length, blueprint):
    if(len(blueprint) == 0):
        return blueprint.ljust(strand_length, 'o')
    else:
        return blueprint



#####################

def get_first_five_bases(blueprint):
    possibilities = util.get_pentameric_possiblilities(blueprint[0:5])

    if(len(possibilities)==0):
        print("no possibilities for first five bases")
        return "none"
    return random.choice(possibilities)


def get_next_base( prev_4, prev_6, blueprint, blueprint_violation_array,curr_length,complement_desired):
    global restricted_sequences
    global fives
    global sevens 

    blueprint_score = blueprint_violation_array[curr_length]
    blueprint_base = blueprint[curr_length]   
    next_possible_bases = []


    #print("next expected score : "+str(blueprint_score))
    #print('blueprint base: '+ str(blueprint_base)+ "  index: "+str(curr_length))

    #CASE: next base is blueprint base,
    if(blueprint_base != 'o'):
        #print("new 5 score : "+ str(calc_five_restriction_score(prev_4+blueprint_base))+"  new 7 score :" + str(calc_seven_restriction_score(prev_6+blueprint_base)))
        
        #CHECK: restriction score is the same as base violations array
        new_score = calc_five_restriction_score(prev_4+blueprint_base) + calc_seven_restriction_score(prev_6+blueprint_base)        
        if(new_score != blueprint_score):
            return []
        next_possible_bases.append(blueprint_base)
        return next_possible_bases

    #CASE: next base is NOT blueprint base but unit contain blueprint bases
    for base in 'ACGT':
       
        pentameric_unit = prev_4 + base         # test pentameric unit
        septameric_unit = prev_6 + base         # test septameric unit
        p_ = fives.get(pentameric_unit, False)   # pentameric unit exists
        s_ = sevens.get(septameric_unit, False)  # septameric unit exists


        # CHECK: new unit's reverse complement exists
        comp_bad = False
        if complement_desired:
            five_comp = util.reverse_complement(pentameric_unit)
            comp_bad = comp_bad or fives.get(five_comp, False)

            seven_comp = util.reverse_complement(septameric_unit)
            comp_bad = comp_bad or sevens.get(seven_comp, False)
        if comp_bad and blueprint_score==0: 
            print(" comp bad ")
            continue

        #CHECK: new base creates unwanted restricted sequence
        new_score = calc_five_restriction_score(pentameric_unit) + calc_seven_restriction_score(septameric_unit)
        
        # print("new five score : "+ str(calc_five_restriction_score(pentameric_unit)) + "  new seven score :" + str(calc_seven_restriction_score(septameric_unit)))

        if (new_score == blueprint_score) and ( ((not p_) and (not s_)) or blueprint_score != 0 ) : 
            next_possible_bases.append(base)

    return next_possible_bases





#####################


# generates a new string of size n that doesn't intersect with existing strings
# those other strings are encoded in the fives

def gen_string(strand_length, blueprint, complement_desired):
    global all_strings
    global fives
    global sevens

    #used for revert fives and sevens dictionary after failed attempts
    saved_fives = fives.copy()
    saved_sevens = sevens.copy()

    blueprint = process_blueprint(strand_length, blueprint)
    blueprint_violation_array = get_blueprint_violation_array(blueprint)

    new_strand = []

    backtrack_array = []

    #check if algorithm can start with building first five bases
    starting_bases = get_first_five_bases(blueprint)
    if len(starting_bases) == 0:
        print('No starting bases given constraints of blueprint.')
        return

    #ensure enough unadded units are available?
    
    attempt = 1
    while(attempt < 20):
        
        new_strand = starting_bases
        curr_length = 5
        while (curr_length < strand_length):
            #print("\n  * "+ str(new_strand) + " *  ")

            #get previous four bases for( _ _ _ _ + new base )
            prev_4 = new_strand[len(new_strand) - 4:]   

            #get previous six bases for( _ _ _ _ _ + new base )       
            prev_6 = ''
            if curr_length >= 7:
                prev_6 = new_strand[len(new_strand) - 6:] 

            #get next possible base
            next_possible_bases = get_next_base(prev_4, prev_6, blueprint, blueprint_violation_array, curr_length,complement_desired)
            #print("next bases: "+ str(next_possible_bases))

            #CASE: add possible base (CONTINUE)
            if len(next_possible_bases) > 0:
                chosen_unit = random.choice(next_possible_bases)
                new_strand += chosen_unit  

                #BACKTRACK CASE: add rest of base possiblilities to backtrack array
                next_possible_bases.remove(chosen_unit)
                backtrack_array.append(next_possible_bases)

                update_fives(prev_4 + chosen_unit, complement_desired)
                update_sevens(prev_6 + chosen_unit, complement_desired)

            #CASE: no possible base (STOP)
            else:
                curr_length = strand_length 
                new_strand = []
            curr_length += 1
        
        if len(new_strand) == strand_length:
            print("order completed at attempt #: " + str(attempt))     
            all_strings.append(new_strand)
            return new_strand
        else:
            #print("attempt failed: "+ str(attempt))     
            
            #undo changes to dictionary during failed attempt
            fives = saved_fives.copy()
            sevens = saved_sevens.copy()
            attempt += 1

    return "Could not generate a string in ", attempt, " attempts."
  



#####################


# updates the fives data structure with the new string
def update_fives(new_unit, complement_desired):
    global fives
    if len(new_unit) == 5:
        if complement_desired:
            fives[new_unit] = True

        new_reverse_unit = util.reverse_complement(new_unit)
        fives[new_reverse_unit] = True


# updates the sevens data structure with the new string
def update_sevens(new_unit, complement_desired):
    global sevens
    if len(new_unit) == 7:
        if complement_desired:
            add_sevens_approx(new_unit)

        new_reverse_unit = util.reverse_complement(new_unit)       
        add_sevens_approx(new_reverse_unit)

def add_sevens_approx(new_unit):
    global sevens
    lets = "ACGT"
    i = 1
    while i < 6:
        for let in lets:
            new = new_unit[:i] + let + new_unit[i + 1:]
            if new in sevens and (not sevens[new]):
                sevens[new] = True
        i += 1

#Include all restricted sequences in the dictionary from the start



########          MAIN METHOD        ############## 

fives = {}
sevens = {}
all_strings = []
size_of_strand = 50

genfives()
gensevens()

test_str = "GGoATooooAAAoooooooATTooGGGGGo"
test_str_2 = "oCCCooooooooooooooATTooGGGGGGo"
#print(get_blueprint_violation_array(test_str))

test = gen_string(30,test_str, True)
#test2 = gen_string(30, test_str_2, True)
print(test_str)
print(test)
#print(test_str_2)
#print(test2)


"""
for i in range(0, 10):
    print gen_string(size_of_strand,"",True)
    print("pentameric units remaining:" + str(fives.values().count(False)))
    print("septameric units remaining:" + str(sevens.values().count(False)))   
"""
