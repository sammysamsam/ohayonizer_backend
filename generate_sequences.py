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

                        if (util.palindrome_check(new_five)):
                            print(new_five)

                        for sequence in ['AAAA', 'TTTT', 'CCC', 'GGG']:
                            if sequence in new_five:
                                restricted = True
                                ##print(new_five + " is restricted.")
                                break
                        # NO POLYPURINES    
                        if (not restricted):
                            fives[new_five] = False



    #print("Length of fives: " + str(len(fives)) )

def gen_nmers(prev_nmers):
    ''' Returns the dictionary of n+1mers using the nmer list passed in '''
    letters = 'ACGT'
    new_nmers = {}
    for nmer in prev_nmers:
        for letter in letters:
            new_nmer = nmer + letter
            restricted = False

            if (util.palindrome_check(new_nmer)):
                print(new_nmer)

            for sequence in ['AAAA', 'TTTT', 'CCC', 'GGG']:
                if sequence in new_nmer:
                    restricted = True
                    ##print(new_five + " is restricted.")
                    break
            # NO POLYPURINES    
            if (not restricted):
                new_nmers[new_nmer] = False
    return new_nmers


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

def get_first_six_bases(blueprint, blueprint_violation_array, complement_desired, front_edges):
    possibilities = []
    if len(blueprint) > 6:
        possibilities = util.get_pentameric_possiblilities(blueprint[0:5], complement_desired)

    if(len(possibilities)==0):
        #print("no possibilities for first five bases")
        return "none"
    return random.choice(possibilities)

def get_last_six_bases(sequence, blueprint, blueprint_violation_array, complement_desired, back_edges):




def get_next_base(prev_4, prev_6, blueprint, blueprint_violation_array, curr_seq, curr_length, complement_desired, nmers, mmers):
    global restricted_sequences

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
        n_unit = prev_4 + base         # test pentameric unit
        m_unit = prev_6 + base         # test septameric unit
        
        #CHECK: new base creates unwanted restricted sequence
        new_score = util.get_new_restriction_score(curr_seq, base, complement_desired)        

        n_ = nmers.get(n_unit, False)   # pentameric unit exists
        m_ = mmers.get(m_unit, False)  # septameric unit exists

        # CHECK: new unit's reverse complement exists
        if complement_desired:
            n_comp = util.reverse_complement(n_unit)
            m_comp = util.reverse_complement(m_unit)
            if nmers.get(five_comp, False) or mmerss.get(seven_comp, False):
                #print(" comp bad ")    
                continue
        
        #print("new score : " + str(new_score) + "|  " + n_unit +":"+ str(p_) + "  || "+ m_unit + ":" + str(s_))

        if (new_score == blueprint_score) and (not n_) and (not m_): 
            next_possible_bases.append(base)

    #print(str(next_possible_bases)+"\n")
    return next_possible_bases


#####################


# generates a new string of size n that doesn't intersect with existing strings
# those other strings are encoded in the fives

def gen_string(strand_length, blueprint, complement_desired, strand_list, front_edges = [], back_edges = []):
    global all_strings
    global fives
    global sevens

    blueprint = process_blueprint(strand_length, blueprint)
    blueprint_violation_array = get_blueprint_violation_array(blueprint, complement_desired)
    #print("\n order processed\nblueprint: "+str(blueprint)+"\nblueprint violation array: "+str(blueprint_violation_array)+"\n\n")

    new_strand = []

    #ensure enough unadded units are available?
    for n in range(5, 10): # start at 5s and then if after 200 attempts doesn't work go forth
        m = n + 2
        attempt = 1
        while(attempt < 200):

            #BACKTRACK VARIABLES
            backtrack_array = []
            backtrack_seq = [] 

            
    ###########
            added_units = {0:[],1:[],2:[],3:[],4:[]}
            new_strand = get_first_six_bases(blueprint, blueprint_violation_array , complement_desired)
    ###########

            curr_length = len(new_strand)
            if curr_length == 0:
                return 'No starting bases given constraints of blueprint.'

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

                    #update backtrack variables
                    next_possible_bases.remove(chosen_unit)
                    if(len(next_possible_bases) > 0):
                        backtrack_array.append(next_possible_bases)
                        backtrack_seq.append(new_strand)
                    
                    #add new strand
                    new_strand += chosen_unit 
                    
                    #update added_bases list
                    added_units[curr_length] = update_nmers(prev_4 + chosen_unit, complement_desired, nmers, n) + update_nmers(prev_6 + chosen_unit, complement_desired, mmers, m)

                #CASE: not possible base (BACKTRACK)
                else:
                    #print("# "+ str(new_strand))

                    backtrack_index = len(backtrack_array)-1
                    if(backtrack_index != -1):
                       # print("backtrack length: "+str(len(backtrack_array)))
                        
                        revert_units(new_strand, backtrack_seq[backtrack_index], added_units ,complement_desired)
                        backtrack_bases = backtrack_array[backtrack_index]
                        backtrack_unit = random.choice(backtrack_bases)

                        #BACKTRACK SEQUENCE                 
                        new_strand = backtrack_seq[backtrack_index]
                        curr_length = len(new_strand)

                        #replace removed base with alternate base
                        prev_4 = new_strand[len(new_strand) - 4:] 
                        prev_6 = ''                   
                        if curr_length >= 6:
                            prev_6 = new_strand[len(new_strand) - 6:] 
                        new_strand += backtrack_unit

                        #update added units
                        added_units[curr_length] = update_nmers(prev_4 + backtrack_unit, complement_desired, nmers, n) + update_nmers(prev_6 + backtrack_unit, complement_desired, mmers, m)

                        #update backtrack variables
                        if(len(backtrack_bases) > 1):
                            backtrack_bases.remove(backtrack_unit)
                            backtrack_array[backtrack_index] = backtrack_bases
                        else:
                            backtrack_seq.pop()
                            backtrack_array.pop()

                        #print("\nBacktrack to: " + str(new_strand))
                    else:
                        #backtrack added 5/7 units
                        print("stoped at: "+new_strand)
                        print("***ATTEMPT #"+str(attempt)+"***\n")
                        revert_units(new_strand, "", added_units, complement_desired, nmers, mmers, n)
                        new_strand = []
                        break

                curr_length += 1
            
            if len(new_strand) == strand_length:
                #print("\norder completed at attempt #: " + str(attempt))
                all_strings.append(new_strand)
                #revert_units(new_strand, "", added_units, complement_desired)
                return str(new_strand)
            else:
                attempt += 1

    return "Could not generate a string in "+ str(attempt) + " attempts."
  

#####################

def revert_units(curr_seq, rev_seq, added_bases, complement_desired, nmers, mmers, n):

    curr_len = len(curr_seq)
    rev_len = len(rev_seq)

    #print("backtrack to sequence:\n# "+rev_seq+"\n")
    for curr_index in range( rev_len, curr_len):
        for unit in added_bases[curr_index]:

            #print("remove: "+unit)
            if len(unit) == n:
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

def update_nmers(new_unit, complement_desired, nmers, n):
    ''' n is length of nmers '''
    unit_list = []
    if len(new_unit) == n:
        if complement_desired and (not nmers.get(new_unit, True)):
            unit_list.extend(update_nmers_approx(new_unit, nmers, n))

        new_reverse_unit = util.reverse_complement(new_unit)
        if (not nmers.get(new_reverse_unit, True)):
            unit_list.extend(update_nmers_approx(new_reverse_unit, nmers, n))
    return unit_list


def update_nmers_approx(new_unit, nmers, n):
    unit_list = []
    lets = "ACGT"
    i = 1
    while i < n:
        for let in lets:
            new = new_unit[:i] + let + new_unit[i + 1:]
            if (not nmers.get(new, True)):
                nmers[new] = True
                unit_list.append(new)
        i += 1
    return unit_list


#Include all restricted sequences in the dictionary from the start



########          MAIN METHOD        ############## 

fives = {}
sevens = {}
all_strings = []
size_of_strand = 50

genfives()
gensevens()
print ("total fives: "+ str(len(fives)))
print ("total sevens: "+ str(len(sevens)) + "\n")

nmers = gen_nmers(fives)
mmers = gen_nmers(sevens)

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

test2 = gen_string(100, "", True)
print("\nresult: " + test2)
print("length of fives used: "+ str(sum(fives.values())))
print("length of sevens used: "+ str(sum(sevens.values())))

test3 = gen_string(100, "", True)
print("\nresult: " + test3)
print("length of fives used: "+ str(sum(fives.values())))
print("length of sevens used: "+ str(sum(sevens.values())))

test4 = gen_string(100, "", True)
print("\nresult: " + test4)
print("length of fives used: "+ str(sum(fives.values())))
print("length of sevens used: "+ str(sum(sevens.values())))

test5 = gen_string(100, "", True)
print("\nresult: " + test5)
print("length of fives used: "+ str(sum(fives.values())))
print("length of sevens used: "+ str(sum(sevens.values())))

print("\n___________\n1: " + test2)
print("2: " + test3)
print("3: " + test4)
print("4: " + test5)

"""
#test3
for i in range(0, 10):
    #print gen_string(size_of_strand,"",True)
    #print("pentameric units remaining:" + str(fives.values().count(False)))
    #print("septameric units remaining:" + str(sevens.values().count(False)))   
"""
