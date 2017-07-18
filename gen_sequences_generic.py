"""

This code (python3) designs strand sequences

Notes:
    - no polypurine (AAAA/TTTT) and  + CCC/GGG
    - avoid alternating CGs and other purines and pyrmadines (3-4)
    - no palindromic units
"""


import random
import itertools
import strand_utilities as util


def gen_nmers(n_size):
    """
    generates all possible n_size nucleotide units as key and values TRUE (if restricted) or FALSE (if valid unit)

    :param n_size: int
    :return:
    """

    bases = ['A', 'T', 'C', 'G']
    nmers = []
    for nmer in itertools.combinations_with_replacement(bases, n_size):
        nmers.append(nmer)

    result = {}
    for new_nmer in nmers:
        if (util.palindrome_check(new_nmer)) or new_nmer in {'AAAA', 'TTTT', 'CCC', 'GGG'}:
            result[new_nmer] = True
        else:
            result[new_nmer] = False
    return result


#####################


def get_blueprint_violation_array(blueprint, complement_desired):
    """
    Returns list where index = index of base in blueprint and value = sum of violation score of base with corresponding
    index and all new violation scores from previous bases.

    :param blueprint: str
    :param complement_desired: bool
    :return: int[]

    """

    violation_score_at_all_positions = []
    curr_blueprint_seq = ""
    total = 0
    for base in blueprint:
        # score at each position is the cumulative score as each base is added
        total += util.get_new_restriction_score(curr_blueprint_seq, base, complement_desired)
        curr_blueprint_seq += base

        violation_score_at_all_positions.append(total)
    return violation_score_at_all_positions


def process_blueprint(strand_length, blueprint):
    """
    Prelim check on blueprint to make sure it is valid.

    :param strand_length: int
    :param blueprint: string[]
    :return: string
    """
    bases = {"A", "T", "C", "G", "o"}
    if any((b not in bases) for b in blueprint):
        raise ValueError('blueprint: ' + blueprint + ' contains characters other than A,T,C,G,o')

    if len(blueprint) == 0:
        return blueprint.ljust(strand_length, 'o')
    else:
        return blueprint


# ===================================================

def get_next_base(n, m, blueprint, blueprint_violation_array, curr_seq, curr_length, complement_desired):
    """
    Calculates the possible bases that can be added without adding unwanted restricted sequences.

    :param n_mer: string
    :param m_mer: string
    :param blueprint: string
    :param blueprint_violation_array: int[]
    :param curr_seq:
    :param curr_length:
    :param complement_desired:
    :return: possible bases (string[])
    """

    blueprint_score = blueprint_violation_array[curr_length] - blueprint_violation_array[curr_length-1]
    blueprint_base = blueprint[curr_length]

    # print("next expected score : "+str(blueprint_score))
    # print('blueprint base: '+ str(blueprint_base)+ "  index: "+str(curr_length))

    # case: specific base desired from blueprint
    if blueprint_base is not 'o':
        new_score = util.get_new_restriction_score(curr_seq, blueprint_base, complement_desired)       
        if new_score is not blueprint_score:
            return []
        return [blueprint_base]

    # case: next base is NOT blueprint base but unit contain blueprint bases
    next_possible_bases = []
    for base in 'ACGT':    
        n_unit = curr_seq[n-1:] + base         # test new n unit
        m_unit = curr_seq[m-1:] + base         # test new m unit
        
        # new base creates unwanted restricted sequence
        new_score = util.get_new_restriction_score(curr_seq, base, complement_desired)        

        if nmers.get(n_unit, False) or mmers.get(m_unit, False):
            continue

        # CHECK: new unit's reverse complement exists
        if complement_desired:
            n_comp = util.reverse_complement(n_unit)
            m_comp = util.reverse_complement(m_unit)
            if nmers.get(n_comp, False) or mmers.get(m_comp, False):
                continue
        
        #print("new score : " + str(new_score) + "|  " + n_unit +":"+ str(p_) + "  || "+ m_unit + ":" + str(s_))

        if (new_score == blueprint_score) and (not n_) and (not m_): 
            next_possible_bases.append(base)

    #print(str(next_possible_bases)+"\n")
    return next_possible_bases


#####################




def gen_strand(strand_length, blueprint, complement_desired, n, m):
    """
    generates a new string of size n that doesn't intersect with existing strings those other strings are encoded in
    the fives

    :param strand_length: int
    :param blueprint: str[]
    :param complement_desired: bool
    :param n: int
    :param m: int
    :return: str
    """
    global all_strings
    global nmers
    global sevens

    blueprint = process_blueprint(strand_length, blueprint)
    blueprint_violation_array = get_blueprint_violation_array(blueprint, complement_desired)

    new_strand = []
    for n_size in range(4, 10): # start at 5s and then if after 200 attempts doesn't work go forth
        m_size = n_size + 2

        for attempt in range(1,100):
            backtrack_array = []
            backtrack_seq = [] 
            added_units = {}


            next_possible_base = get_next_base(n_size, m_size, blueprint, blueprint_violation_array,
                                               new_strand, complement_desired, front_edges, complement_front_edges)

            # CASE: add possible base (CONTINUE)
            if next_possible_base != "":
                new_strand += next_possible_base

            # CASE: not possible base (BACKTRACK)
            else:
                backtrack_index = len(backtrack_array) - 1

                if(backtrack_index != -1):

                    revert_units(new_strand, backtrack_seq[backtrack_index], added_units ,complement_desired)
                    backtrack_bases = backtrack_array[backtrack_index]
                    backtrack_unit = random.choice(backtrack_bases)

                    # Backtrack Sequence
                    new_strand = backtrack_seq[backtrack_index]
                    curr_length = len(new_strand)

                    # Replace removed base with alternate base
                    prev_4 = new_strand[len(new_strand) - 4:] 
                    prev_6 = ''                   
                    if curr_length >= 6:
                        prev_6 = new_strand[len(new_strand) - 6:] 
                    new_strand += backtrack_unit

                    # Update added units
                    added_units[curr_length] = update_fives(prev_4 + backtrack_unit, complement_desired) + update_sevens(prev_6 + backtrack_unit, complement_desired)

                    # Update backtrack variables
                    if (len(backtrack_bases) > 1):
                        backtrack_bases.remove(backtrack_unit)
                        backtrack_array[backtrack_index] = backtrack_bases
                    else:
                        backtrack_seq.pop()
                        backtrack_array.pop()
                else:
                    # print("stoped at: "+new_strand + "\n***ATTEMPT #"+str(attempt)+"***\n")
                    new_strand = ""
                    break

        if len(new_strand) == strand_length:
            # and (( backedge_check(new_strand, back_edges, complement_back_edges)) or not 'o' in blueprint):
            return str(new_strand)

        revert_units(new_strand, "", added_units, complement_desired)             
        attempt += 1

    # print( "Could not generate a string in "+ str(attempt) + " attempts." )
    return ""
  

#####################

def revert_units(curr_seq, rev_seq, added_bases, complement_desired):
    """
    remove the added nmer and mmer units that were added to rev_seq to create curr_seq. (this
    :param curr_seq:
    :param rev_seq:
    :param added_bases:
    :param complement_desired:
    :return:
    """
    global fives
    global sevens

    curr_len = len(curr_seq)
    rev_len = len(rev_seq)

    # print("backtrack to sequence:\n# "+rev_seq+"\n")
    for curr_index in range( rev_len, curr_len):
        for unit in added_bases[curr_index]:

            # print("remove: "+unit)
            if len(unit) == 5:
                fives[unit] = False
            else:
                sevens[unit] = False


def update_nmers(new_unit, complement_desired, nmers, n):

    unit_list = []
    if len(new_unit) == n:
        if complement_desired and (not nmers.get(new_unit, True)):
            unit_list.extend(update_nmers_approx(new_unit, nmers, n))

        new_reverse_unit = util.reverse_complement(new_unit)
        if not nmers.get(new_reverse_unit, True):
            unit_list.extend(update_nmers_approx(new_reverse_unit, nmers, n))
    return unit_list


def update_nmers_approx(new_unit, nmers, n):
    unit_list = []
    lets = "ACGT"
    i = 1
    while i < n:
        for let in lets:
            new = new_unit[:i] + let + new_unit[i + 1:]
            if not nmers.get(new, True):
                nmers[new] = True
                unit_list.append(new)
        i += 1
    return unit_list


#Include all restricted sequences in the dictionary from the start



########          MAIN METHOD        ############## 
#
# fives = {}
# sevens = {}
# all_strings = []
# size_of_strand = 50
#
# genfives()
# gensevens()
#
# print ("total fives: "+ str(len(fives)))
# print ("total sevens: "+ str(len(sevens)) + "\n")
#
# nmers = gen_nmers(fives)
# mmers = gen_nmers(sevens)


# # test1
# test_str = "GGoATooooAAAAooooAoAToGGoGGGoGooooATGCoo"
# test_str_2 = "oCCCooooooooooooooATTooGGGoooGGGooooooooooo"
#
# print("blueprint violation array: " + str(get_blueprint_violation_array(test_str)))
#
# print ("total fives: "+ str(len(fives)))
# test2 = gen_string(40, test_str, True)
# print(test2)
# print("\n\nlength of fives used: "+ str(sum(fives.values())))
#
# print ("total fives: "+ str(len(fives)))
# test = gen_string(43, test_str_2, True)
# print(test)
# print("\n\nlength of fives used: "+ str(sum(fives.values())))



#
# #test2
#
# test2 = gen_string(100, "", True)
# print("\nresult: " + test2)
# print("length of fives used: "+ str(sum(fives.values())))
# print("length of sevens used: "+ str(sum(sevens.values())))
#
# test3 = gen_string(100, "", True)
# print("\nresult: " + test3)
# print("length of fives used: "+ str(sum(fives.values())))
# print("length of sevens used: "+ str(sum(sevens.values())))
#
# test4 = gen_string(100, "", True)
# print("\nresult: " + test4)
# print("length of fives used: "+ str(sum(fives.values())))
# print("length of sevens used: "+ str(sum(sevens.values())))
#
# test5 = gen_string(100, "", True)
# print("\nresult: " + test5)
# print("length of fives used: "+ str(sum(fives.values())))
# print("length of sevens used: "+ str(sum(sevens.values())))
#
# print("\n___________\n1: " + test2)
# print("2: " + test3)
# print("3: " + test4)
# print("4: " + test5)


# test3
# for i in range(0, 10):
#     #print gen_string(size_of_strand,"",True)
#     #print("pentameric units remaining:" + str(fives.values().count(False)))
#     #print("septameric units remaining:" + str(sevens.values().count(False)))
