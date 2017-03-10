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
                                ##print(new_five + " is restricted.")
                                break
                        # NO PALINDROMES, POLYPURINES    
                        if (not restricted) and (not new_five == new_five[::-1]) :
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
                                for sequence in util.full_restricted_seq:
                                    if sequence in new_seven:
                                        restricted = True
                                        break
                                
                                # NO PALINDROMES, ALTERNATING PURINE-PYRIMIDINE (vice versa)
                                if (not restricted) and (not new_seven== new_seven[::-1]):
                                    sevens[new_seven] = False
    #print("Length of sevens: " + str(len(sevens)) + "\n")


#####################


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

def get_blueprint_violation_array(blueprint,complement_desired):
    violation_array = []
    curr_blueprint_seq = ""
    total = 0
    for base in blueprint:
        ##print(curr_blueprint_seq)

        total += util.get_new_restriction_score(curr_blueprint_seq, base, complement_desired)
        curr_blueprint_seq += base

        violation_array.append(total)

    return violation_array



def process_blueprint(strand_length, blueprint):
    if(len(blueprint) == 0):
        return blueprint.ljust(strand_length, 'o')
    else:
        return blueprint



#####################

def get_first_five_bases(blueprint, complement_desired):
    possibilities = []
    if len(blueprint) > 4:
        possibilities = util.get_pentameric_possiblilities(blueprint[0:5], complement_desired)

    if(len(possibilities)==0):
        #print("no possibilities for first five bases")
        return "none"
    return random.choice(possibilities)


def get_next_base(prev_4, prev_6, blueprint, blueprint_violation_array, curr_seq, curr_length, complement_desired):
    global restricted_sequences
    global fives
    global sevens 

    blueprint_score = blueprint_violation_array[curr_length]-blueprint_violation_array[curr_length-1]
    blueprint_base = blueprint[curr_length]   
    #print("next expected score : "+str(blueprint_score))
    #print('blueprint base: '+ str(blueprint_base)+ "  index: "+str(curr_length))

    #CASE: specific base desired from blueprint
    if(blueprint_base != 'o'):
        new_score = util.get_new_restriction_score(curr_seq, blueprint_base, complement_desired)       
        
        if(new_score != blueprint_score):  
            return []
        return [blueprint_base]

    #CASE: next base is NOT blueprint base but unit contain blueprint bases
    next_possible_bases = []
    for base in 'ACGT':    
        pentameric_unit = prev_4 + base         # test pentameric unit
        septameric_unit = prev_6 + base         # test septameric unit
        
        p_ = fives.get(pentameric_unit, False)   # pentameric unit exists
        s_ = sevens.get(septameric_unit, False)  # septameric unit exists

        # CHECK: new unit's reverse complement exists
        if complement_desired:
            five_comp = util.reverse_complement(pentameric_unit)
            seven_comp = util.reverse_complement(septameric_unit)
            if fives.get(five_comp, False) or sevens.get(seven_comp, False):
                #print(" comp bad ")    
                continue

        #CHECK: new base creates unwanted restricted sequence
        new_score = util.get_new_restriction_score(curr_seq, base, complement_desired)
        
        #print("new score : " + str(new_score) + "|  " + pentameric_unit +":"+ str(p_) + "  || "+ septameric_unit + ":" + str(s_))

        if (new_score == blueprint_score) and (not p_) and (not s_): 
            next_possible_bases.append(base)

    #print(str(next_possible_bases)+"\n")
    return next_possible_bases





#####################


# generates a new string of size n that doesn't intersect with existing strings
# those other strings are encoded in the fives

def gen_string(strand_length, blueprint, complement_desired):
    global all_strings
    global fives
    global sevens

    blueprint = process_blueprint(strand_length, blueprint)
    blueprint_violation_array = get_blueprint_violation_array(blueprint, complement_desired)
    #print("\n order processed\nblueprint: "+str(blueprint)+"\nblueprint violation array: "+str(blueprint_violation_array)+"\n\n")

    new_strand = []

    #ensure enough unadded units are available?
    
    attempt = 1
    while(attempt < 20):

        #BACKTRACK VARIABLES
        backtrack_array = []
        backtrack_seq = []      

        new_strand = get_first_five_bases(blueprint, complement_desired)
        curr_length = len(new_strand)

        #check if algorithm can start with building first five bases       
        if curr_length == 0:
            #print('No starting bases given constraints of blueprint.')
            return

        while (curr_length < strand_length):

            #get previous four/six bases for( _ _ _ _ + new base )
            prev_4 = new_strand[len(new_strand) - 4:] 
            prev_6 = ''                   
            if curr_length >= 6:
                prev_6 = new_strand[len(new_strand) - 6:] 

            next_possible_bases = get_next_base(prev_4, prev_6, blueprint, blueprint_violation_array, new_strand, curr_length,complement_desired)

            #CASE: add possible base (CONTINUE)
            if len(next_possible_bases) > 0:
                chosen_unit = random.choice(next_possible_bases)

                #BACKTRACK: add rest of base possiblilities to backtrack array
                next_possible_bases.remove(chosen_unit)
                if(len(next_possible_bases) > 0):
                    backtrack_array.append(next_possible_bases)
                    backtrack_seq.append(new_strand)
                
                #add new strand
                new_strand += chosen_unit 
                update_fives(prev_4 + chosen_unit, complement_desired)
                update_sevens(prev_6 + chosen_unit, complement_desired)
                
                #print("# "+ str(new_strand))

            #CASE: not possible base (BACKTRACK)
            else:
                backtrack_index = len(backtrack_array)-1
                if(backtrack_index != -1):
                    #print("\n\nBACKTRACK \n go back " + str(len(new_strand) - len(backtrack_seq[backtrack_index])) + " sequences")
                    #print("UNITS REVERTED FROM " + str( sum(fives.values())))
                    
                    #backtrack 5/7 unit dictionary
                    revert_units(new_strand, backtrack_seq[backtrack_index], blueprint ,complement_desired)

                    #backtrack sequence
                    new_strand = backtrack_seq[backtrack_index]

                    #pick new base and update backtrack_array
                    backtrack_bases = backtrack_array[backtrack_index]
                    backtrack_unit = random.choice(backtrack_bases)

                    if(len(backtrack_bases) > 1):
                        backtrack_bases.remove(backtrack_unit)
                        backtrack_array[backtrack_index] = backtrack_bases
                    else:
                        backtrack_seq.pop()
                        backtrack_array.pop()

                    #update newly created five and seven units
                    prev_4 = new_strand[len(new_strand) - 4:] 
                    prev_6 = ''                   
                    if curr_length >= 6:
                        prev_6 = new_strand[len(new_strand) - 6:] 
                    update_fives(prev_4 + backtrack_unit, complement_desired)
                    update_sevens(prev_6 + backtrack_unit, complement_desired)

                    #add new base to new strand
                    new_strand += backtrack_unit
                    curr_length = len(new_strand)
                   
                    #print('TO '+str(sum(fives.values())))
                    #print("\n# "+ str(new_strand)+'\n')

                    continue
                else:
                    #backtrack added 5/7 units
                    revert_units(new_strand, "", blueprint, complement_desired)         
                    new_strand = []
                    break

            curr_length += 1
        
        if len(new_strand) == strand_length:
            #print("\norder completed at attempt #: " + str(attempt))     
            all_strings.append(new_strand)
            return new_strand
        else:
            attempt += 1

    return "Could not generate a string in ", attempt, " attempts."
  

#####################

def revert_units(curr_seq, rev_seq, blueprint, complement_desired):
    global fives
    global sevens

    curr_len = len(curr_seq)
    rev_len = len(rev_seq)
    #print("backtrack seq : "+rev_seq)
    #print("curr seq:       "+curr_seq)
    for curr_index in range( rev_len, curr_len):
        #print("\nremove base: "+ curr_seq[curr_index])
        
        #base is blueprint base
        if(blueprint[curr_index] != 'o'):
            continue 

        if(curr_index > 4):
            pent = curr_seq[curr_index - 5:curr_index]
            if complement_desired and fives.get(pent, False) :
                fives[pent] = False
            
            p_reverse_unit = util.reverse_complement(pent)   
            if fives.get(p_reverse_unit, False):
                fives[p_reverse_unit] = False
            #print("pent: "+pent)
        if(curr_index > 6):
            sept = curr_seq[curr_index - 7:curr_index]   
            if complement_desired and sevens.get(sept, False) :      
                sevens[sept] = False

            s_reverse_unit = util.reverse_complement(sept)
            if sevens.get(s_reverse_unit, False):
                sevens[s_reverse_unit] = False 
            #print("sept: "+sept)

# updates the fives data structure with the new string
def update_fives(new_unit, complement_desired):
    global fives
    if len(new_unit) == 5:
        if complement_desired and (not fives.get(new_unit, True)): #not exist = False, added = True, not added = False
            fives[new_unit] = True

        new_reverse_unit = util.reverse_complement(new_unit)
        if (not fives.get(new_reverse_unit, True)):
            fives[new_reverse_unit] = True


# updates the sevens data structure with the new string
def update_sevens(new_unit, complement_desired):
    global sevens
    if len(new_unit) == 7:
        if complement_desired and (not sevens.get(new_unit, True)):
            update_sevens_approx(new_unit)

        new_reverse_unit = util.reverse_complement(new_unit) 
        if (not sevens.get(new_reverse_unit, True)):      
            update_sevens_approx(new_reverse_unit)

def update_sevens_approx(new_unit):
    global sevens
    lets = "ACGT"
    i = 1
    while i < 6:
        for let in lets:
            new = new_unit[:i] + let + new_unit[i + 1:]
            if new in sevens and (not sevens.get(new, True) ):
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

test_str = "GGoATooooAAAAooooAoAToGGoGGGoGooooATGCoo"
test_str_2 = "oCCCooooooooooooooATTooGGGoooGGGooooooooooo"
##print(get_blueprint_violation_array(test_str))



test3 = gen_string(100, "",True)
test3 = gen_string(100, "",True)
test3 = gen_string(100, "",True)
#test = gen_string(40,test_str, True)
test2 = gen_string(43, test_str_2, True)
#test3 = gen_string(100,"",True)

##print(test_str)

##print(str(test2)+"\n")
##print(test_str_2)



#test4 = gen_string(100,"",True)


print(test2)
print(test3)
##print(test4)


"""
for i in range(0, 10):
    #print gen_string(size_of_strand,"",True)
    #print("pentameric units remaining:" + str(fives.values().count(False)))
    #print("septameric units remaining:" + str(sevens.values().count(False)))   
"""
