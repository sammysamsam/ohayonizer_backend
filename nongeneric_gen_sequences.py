"""
This software designs strand sequences, including its Zx[ bvcx complements, with unique pentameric (5) and 
septameric (7) units, i.e. AATCG & TAGGTCC.
in addition:
    - no polypurine (AAAA/TTTT) and  + CCC/GGG
    - avoid alternating CGs and other purines and pyrmadines (3-4)
    - no palindromic units
"""


import random
import pprint
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

                        if (util.palindrome_check(new_five)):
                            print(new_five)

                        for sequence in ['AAAA', 'TTTT', 'CCC', 'GGG'] :
                            if sequence in new_five:
                                restricted = True
                                ##print(new_five + " is restricted.")
                                break
                        # NO POLYPURINES    
                        if (not restricted):
                            fives[new_five] = False

    #print("Length of fives: " + str(len(fives)) )

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
                                if (util.palindrome_check(new_seven)):
                                    print(new_seven)

                                for sequence in util.full_restricted_seq:
                                    if sequence in new_seven:
                                        restricted = True
                                        break
                                
                                # NO ALTERNATING PURINE-PYRIMIDINE (vice versa)
                                if (not restricted):
                                    sevens[new_seven] = False



#####################




def get_blueprint_violation_array(blueprint,complement_desired):
    """
        Returns list where index = index of base in blueprint and value = sum of violation score of base with corresponding index and all new violation scores from previous bases.
        - Example: blueprint of "o o o g g g o o o" produces blueprint violation array of [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]

    """

    violation_array = []
    curr_blueprint_seq = ""
    total = 0
    for base in blueprint:
        ##print(curr_blueprint_seq)
        total += util.get_new_restriction_score(curr_blueprint_seq, base)
        curr_blueprint_seq += base

        violation_array.append(total)

    return violation_array

def process_blueprint(strand_length, blueprint):
    if(len(blueprint) == 0):
        return blueprint.ljust(strand_length, 'o')
    else:
        return blueprint



####################


def get_next_base(blueprint, blueprint_violation_array, curr_seq, complement_desired, front_edges, complement_front_edges):
    global fives
    global sevens

    global backtrack_array
    global backtrack_seq
    global added_units

    blueprint_score = blueprint_violation_array[len(curr_seq)]
    if len(curr_seq) > 0:
        blueprint_score -= blueprint_violation_array[len(curr_seq)-1]

    blueprint_base = blueprint[len(curr_seq)]   

    #CASE: specific base desired from blueprint
    if(blueprint_base != 'o'):
        new_score = util.get_new_restriction_score(curr_seq, blueprint_base)       
        if(new_score != blueprint_score):  
            return []
        return [blueprint_base]

    if len(curr_seq) > 6:
        front_edges = ['']
        complement_front_edges = []


    #CASE: next base is NOT blueprint base
    next_possible_bases = set(['A','T','C','G'])

    for e in front_edges:
        combined_edge = e + curr_seq
        prev_4 = combined_edge[-4:]
        prev_6 = combined_edge[-6:]

        for base in next_possible_bases.copy():    
            new_score = util.get_new_restriction_score(combined_edge, base) 

            p_ = fives.get(prev_4 + base, False)   # pentameric unit exists
            s_ = sevens.get(prev_6 + base, False)  # septameric unit exists

            # CHECK: new unit's reverse complement (if needed)
            if complement_desired and (fives.get(util.reverse_complement(prev_4 + base), False) or sevens.get(util.reverse_complement(prev_6 + base), False)):
                next_possible_bases.remove(base)
                continue
            # CHECK: new unit
            if not ((new_score == blueprint_score) and (not p_) and (not s_)): 
                next_possible_bases.remove(base)
                


    for e in complement_front_edges:
        combined_edge = e + util.reverse_complement(curr_seq)
        prev_4 = combined_edge[-4:]
        prev_6 = combined_edge[-6:]

        for base in next_possible_bases.copy():    
            new_score = util.get_new_restriction_score(combined_edge, base)  
            
            p_ = fives.get(prev_4 + base, False)   # pentameric unit exists
            s_ = sevens.get(prev_6 + base, False)  # septameric unit exists

            # CHECK: new unit's reverse complement exists
            if complement_desired and (fives.get(util.reverse_complement(prev_4 + base), False) or sevens.get(util.reverse_complement(prev_6 + base), False)):
                next_possible_bases.remove(base)
                continue
            if not ((new_score == blueprint_score) and (not p_) and (not s_)): 
                next_possible_bases.remove(base)



    next_possible_bases = list(next_possible_bases)
    if len(next_possible_bases) == 0:
        return ""

    chosen_base = random.choice(next_possible_bases)
    added_units[len(curr_seq)] = []

    #update backtrack variables
    for e in front_edges:
        if e == '':     
            added_units[len(curr_seq)] += update_fives( (curr_seq)[-4:] + chosen_base, complement_desired) 
            added_units[len(curr_seq)] += update_sevens( (curr_seq)[-6:] + chosen_base, complement_desired)  
        else:
            added_units[len(curr_seq)] += update_fives( (e + curr_seq)[-4:] + chosen_base, False) 
            added_units[len(curr_seq)] += update_sevens( (e + curr_seq)[-6:] + chosen_base, False)  

    for e in complement_front_edges:
        added_units[len(curr_seq)] += update_fives( ( e + util.reverse_complement(curr_seq))[-4:] + util.complement(chosen_base), False) 
        added_units[len(curr_seq)] += update_sevens( (e + util.reverse_complement(curr_seq))[-6:] + util.complement(chosen_base), False)  

    next_possible_bases.remove(chosen_base)
    if(len(next_possible_bases) > 0):
        backtrack_array.append(next_possible_bases)
        backtrack_seq.append(curr_seq)

    return chosen_base



#####################


#Backtrack Variables
backtrack_array = []
backtrack_seq = [] 
added_units = {}

def gen_string(strand_length, blueprint, complement_desired, front_edges=[], complement_front_edges=[], back_edges=[], complement_back_edges=[]):
    """
        Generates a new string of size n that doesn't intersect with existing strings. Those other strings are encoded in the fives

    """
    global fives
    global sevens

    global backtrack_array
    global backtrack_seq
    global added_units

    backtrack_array = []
    backtrack_seq = [] 
    added_units = {}

    blueprint = process_blueprint(strand_length, blueprint)
    blueprint_violation_array = get_blueprint_violation_array(blueprint, complement_desired)
    new_strand = ""
    
    attempt = 1
    while(attempt < 30):
        while (len(new_strand) < strand_length):
 
            next_possible_base = get_next_base(blueprint, blueprint_violation_array, new_strand, complement_desired, front_edges, complement_front_edges)

            #CASE: add possible base (CONTINUE)
            if next_possible_base != "":
                new_strand += next_possible_base

            #CASE: not possible base (BACKTRACK)
            else:
                backtrack_index = len(backtrack_array) - 1

                if(backtrack_index != -1):

                    revert_units(new_strand, backtrack_seq[backtrack_index], added_units ,complement_desired)
                    backtrack_bases = backtrack_array[backtrack_index]
                    backtrack_unit = random.choice(backtrack_bases)

                    #Backtrack Sequence              
                    new_strand = backtrack_seq[backtrack_index]
                    curr_length = len(new_strand)

                    #Replace removed base with alternate base
                    prev_4 = new_strand[len(new_strand) - 4:] 
                    prev_6 = ''                   
                    if curr_length >= 6:
                        prev_6 = new_strand[len(new_strand) - 6:] 
                    new_strand += backtrack_unit

                    #Update added units
                    added_units[curr_length] = update_fives(prev_4 + backtrack_unit, complement_desired) + update_sevens(prev_6 + backtrack_unit, complement_desired)

                    #Update backtrack variables
                    if(len(backtrack_bases) > 1):
                        backtrack_bases.remove(backtrack_unit)
                        backtrack_array[backtrack_index] = backtrack_bases
                    else:
                        backtrack_seq.pop()
                        backtrack_array.pop()
                else:
                    #print("stoped at: "+new_strand + "\n***ATTEMPT #"+str(attempt)+"***\n")       
                    new_strand = ""
                    break

        if len(new_strand) == strand_length:
            #print("\norder completed at attempt #: " + str(attempt))                
            return str(new_strand)

        revert_units(new_strand, "", added_units, complement_desired)             
        attempt += 1

    #print( "Could not generate a string in "+ str(attempt) + " attempts." )
    return ""
  

#####################

def revert_units(curr_seq, rev_seq, added_bases, complement_desired):
    global fives
    global sevens

    curr_len = len(curr_seq)
    rev_len = len(rev_seq)

    #print("backtrack to sequence:\n# "+rev_seq+"\n")
    for curr_index in range( rev_len, curr_len):
        for unit in added_bases[curr_index]:

            #print("remove: "+unit)
            if len(unit) == 5:
                fives[unit] = False
            else:
                sevens[unit] = False



# updates the fives data structure with the new string
def update_fives(new_unit, complement_desired):
    global fives
    unit_list = []
    if len(new_unit) == 5:
        if complement_desired and (not fives.get(new_unit, True)): #not exist = True, added = True, not added = False
            fives[new_unit] = True
            unit_list.append(new_unit)

        new_reverse_unit = util.reverse_complement(new_unit)
        if (not fives.get(new_reverse_unit, True)):
            fives[new_reverse_unit] = True
            unit_list.append(new_reverse_unit)

    return unit_list

# updates the sevens data structure with the new string
def update_sevens(new_unit, complement_desired):
    global sevens
    unit_list = []
    if len(new_unit) == 7:
        if complement_desired and (not sevens.get(new_unit, True)):
            unit_list.extend(update_sevens_approx(new_unit))

        new_reverse_unit = util.reverse_complement(new_unit) 
        if (not sevens.get(new_reverse_unit, True)):      
            unit_list.extend(update_sevens_approx(new_reverse_unit))
    return unit_list

def update_sevens_approx(new_unit):
    global sevens
    unit_list = []
    lets = "ACGT"
    i = 1
    while i < 6:
        for let in lets:
            new = new_unit[:i] + let + new_unit[i + 1:]
            if (not sevens.get(new, True)):
                sevens[new] = True
                unit_list.append(new)
        i += 1
    return unit_list

#Include all restricted sequences in the dictionary from the start



########          MAIN METHOD        ############## 

fives = {}
sevens = {}
genfives()
gensevens()

# print ("total fives: "+ str(len(fives)))
# print ("total sevens: "+ str(len(sevens)) + "\n")
"""
#test1
test_str = "GGoATooooAAAAooooAoAToGGoGGGoGooooATGCoo"
test_str_2 = "oCCCooooooooooooooATTooGGGoooGGGooooooooooo"
print("blueprint violation array: " + str(get_blueprint_violation_array(test_str)))
print ("total fives: "+ str(len(fives)))
test2 = gen_string(40, test_str, True)
print(test2)
print("\n\nlength of fives used: "+ str(sum(fives.values())))
print ("total fives: "+ str(len(fives)))
test = gen_string(43, test_str_2, True)
print(test)
print("\n\nlength of fives used: "+ str(sum(fives.values())))
"""


#test2
test2 = gen_string(77, "", True ,["ATAAGC",' ATATA','GGGGGG'])
print(test2)

test3 = gen_string(17, "", True)
print(test3)


# print("\nresult: " + test2)
# print("length of fives used: "+ str(sum(fives.values())))
# print("length of sevens used: "+ str(sum(sevens.values())))

# test3 = gen_string(100, "", True)
# print("\nresult: " + test3)
# print("length of fives used: "+ str(sum(fives.values())))
# print("length of sevens used: "+ str(sum(sevens.values())))

# test4 = gen_string(100, "", True)
# print("\nresult: " + test4)
# print("length of fives used: "+ str(sum(fives.values())))
# print("length of sevens used: "+ str(sum(sevens.values())))

# test5 = gen_string(100, "", True)
# print("\nresult: " + test5)
# print("length of fives used: "+ str(sum(fives.values())))
# print("length of sevens used: "+ str(sum(sevens.values())))

# print("\n___________\n1: " + test2)
# print("2: " + test3)
# print("3: " + test4)
# print("4: " + test5)

"""
#test3
for i in range(0, 10):
    #print gen_string(size_of_strand,"",True)
    #print("pentameric units remaining:" + str(fives.values().count(False)))
    #print("septameric units remaining:" + str(sevens.values().count(False)))   
"""
