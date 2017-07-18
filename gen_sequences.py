"""
This code (python3) designs strand sequences, including its Zx[ bvcx complements, with unique pentameric (5) and
septameric (7) units, i.e. AATCG & TAGGTCC.

in addition:
    - no polypurine (AAAA/TTTT) and  + CCC/GGG
    - avoid alternating CGs and other purines and pyrmadines (3-4)
    - no palindromic units
"""

import random
import strand_utilities as util

BACKTRACK_ARRAY = []
BACKTRACK_SEQ = []
ADDED_UNITS = {}
FIVES = {}
SEVENS = {}


def gen_restricted_sequences():
    """
    generate sequences that are known to cause issues in specific
    :return:
    """
    pyrimadines = 'TC'
    purines = 'AG'

    result = ['AAAA', 'TTTT', 'CCC', 'GGG', 'GGGG', 'CCCC']

    # alternating purine - pyrimidines
    for pur1 in purines:
        for pyr1 in pyrimadines:
            for pur2 in purines:
                for pyr2 in pyrimadines:
                    for pur3 in purines:
                        for pyr3 in pyrimadines:
                            result.append(pur1 + pyr1 + pur2 + pyr2 + pur3 + pyr3)
                            result.append(pyr1 + pur1 + pyr2 + pur2 + pyr3 + pur3)

    return result

RESTRICTED_SEQ = gen_restricted_sequences()


def gen_fives():
    """
    generates all possible non restricted five nucleotide units as key and values = FALSE into
    global FIVES variable

    :return: None
    """

    global FIVES
    lets = 'ACGT'
    for i1 in lets:
        for i2 in lets:
            for i3 in lets:
                for i4 in lets:
                    for i5 in lets:
                        new_five = i1 + i2 + i3 + i4 + i5
                        restricted = False

                        if util.palindrome_check(new_five):
                            print(new_five)

                        for sequence in {'AAAA', 'TTTT', 'CCC', 'GGG'}:
                            if sequence in new_five:
                                restricted = True
                                # print(new_five + " is restricted.")
                                break

                        # only add if not restricted
                        if not restricted:
                            FIVES[new_five] = False

    # print("Length of FIVES: " + str(len(FIVES)) )


def gen_sevens():
    """
    generates all possible non restricted seven nucleotide units as key and value = FALSE into global SEVENS variable

    :return: None
    """

    global SEVENS
    global RESTRICTED_SEQ
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

                                if util.palindrome_check(new_seven):
                                    print(new_seven)

                                for rs in set(RESTRICTED_SEQ):
                                    if rs in new_seven:
                                        restricted = True
                                        break
                                
                                if not restricted:
                                    SEVENS[new_seven] = False


# ===================================================

def get_blueprint_violation_array(blueprint):
    """
    Returns list where index = index of base in blueprint and value = sum of violation score of base with corresponding
    index and all new violation scores from previous bases.

    - Example: blueprint of "o o o g g g o o o" produces blueprint violation array of [0, 0, 0, 0, 0, 1, 1, 1, 1, 1]

    :param blueprint: str
    :return: int[]


    """

    violation_array = []
    curr_blueprint_seq = ""
    total = 0
    for base in blueprint:
        # print(curr_blueprint_seq)
        total += util.get_new_restriction_score(curr_blueprint_seq, base, RESTRICTED_SEQ)
        curr_blueprint_seq += base

        violation_array.append(total)

    return violation_array


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


def get_next_base(blueprint, blueprint_violation_array, curr_seq, complement_desired, front_edges, complement_front_edges):
    """
    Calculates the possible bases that can be added without adding unwanted restricted sequences.

    :param blueprint: string
    :param blueprint_violation_array: int[]
    :param curr_seq: string
    :param complement_desired: bool
    :param front_edges: string[]
    :param complement_front_edges: string[]
    :return: string[]
    """

    global FIVES
    global SEVENS
    global RESTRICTED_SEQ
    global ADDED_UNITS

    # remove front edge consideration when necessary (only matters in the first 6 bases)
    if len(curr_seq) > 6:
        front_edges = ['']
        complement_front_edges = []
    else:
        front_edges.append('')

    # calculate expected restricted score based on blueprint
    blueprint_score = 0
    if len(curr_seq) > 0:
        blueprint_score -= blueprint_violation_array[len(curr_seq)-1]

    blueprint_base = blueprint[len(curr_seq)]   

    # if specific base desired from blueprint
    if blueprint_base != 'o':
        new_score = util.get_new_restriction_score(curr_seq, blueprint_base, RESTRICTED_SEQ)
        if new_score is not blueprint_score:
            return ""
        else:
            return blueprint_base

    # if next base is NOT blueprint base
    else:
        next_possible_bases = {'A', 'T', 'C', 'G'}
        new_score = 0

        for e in front_edges:

            prev_4 = (e + curr_seq)[-4:]
            prev_6 = (e + curr_seq)[-6:]
            prev_blueprint_6 = (e + curr_seq)[-6:]

            for base in next_possible_bases.copy():

                # check if unit has been added before
                if (prev_4 + base) in FIVES or (prev_6 + base) in SEVENS:
                    next_possible_bases.remove(base)
                    continue

                # check if edge + base creates new unwanted restricted sequence
                if util.get_new_restriction_score(prev_6, base, RESTRICTED_SEQ) is not util.get_new_restriction_score(prev_blueprint_6, base, RESTRICTED_SEQ):
                    next_possible_bases.remove(base)
                    continue

        for e in complement_front_edges:
            combined_edge = e + util.reverse_complement(curr_seq)
            prev_4 = combined_edge[-4:]
            prev_6 = combined_edge[-6:]
            prev_blueprint_6 = combined_edge[:-6]

            for base in next_possible_bases.copy():

                if util.reverse_complement(prev_4 + base) in FIVES or util.reverse_complement(prev_6 + base) in SEVENS:
                    next_possible_bases.remove(base)
                    continue

                # check if edge + base creates new unwanted restricted sequence
                if util.get_new_restriction_score(prev_6, base, RESTRICTED_SEQ) is not util.get_new_restriction_score(prev_blueprint_6, base, RESTRICTED_SEQ):
                    next_possible_bases.remove(base)
                    continue

        next_possible_bases = list(next_possible_bases)
        if len(next_possible_bases) == 0:
            return ""
        chosen_base = random.choice(next_possible_bases)
    
        ADDED_UNITS[len(curr_seq)] = []

        # update backtrack variables
        for e in front_edges:
            if e == '':
                ADDED_UNITS[len(curr_seq)] += update_fives( (curr_seq)[-4:] + chosen_base, complement_desired)
                ADDED_UNITS[len(curr_seq)] += update_sevens( (curr_seq)[-6:] + chosen_base, complement_desired)
            else:
                ADDED_UNITS[len(curr_seq)] += update_fives( (e + curr_seq)[-4:] + chosen_base, False)
                ADDED_UNITS[len(curr_seq)] += update_sevens( (e + curr_seq)[-6:] + chosen_base, False)
        for e in complement_front_edges:
            ADDED_UNITS[len(curr_seq)] += update_fives( ( e + util.reverse_complement(curr_seq))[-4:] + util.complement(chosen_base), False)
            ADDED_UNITS[len(curr_seq)] += update_sevens( (e + util.reverse_complement(curr_seq))[-6:] + util.complement(chosen_base), False)

        if(len(next_possible_bases) > 1):
            next_possible_bases.remove(chosen_base)
            BACKTRACK_ARRAY.append(next_possible_bases)
            BACKTRACK_SEQ.append(curr_seq)

    return chosen_base


def backedge_check(curr_seq, back_edges, complement_back_edges):
    global FIVES
    global SEVENS

    for e in back_edges:
        combined_edge = curr_seq + e 
        if FIVES.get(combined_edge[-5:], False) or SEVENS.get(combined_edge[-7:], False):
            return False

        for pos in range(1,len(e)):

            p_ = FIVES.get(combined_edge[:-pos][-5:], False)   # pentameric unit exists
            s_ = SEVENS.get(combined_edge[:-pos][-7:], False)  # septameric unit exists

            if p_ or s_: 
                return False
                
    for e in complement_back_edges:
        combined_edge = util.reverse_complement(curr_seq) + e
        if FIVES.get(combined_edge[-5:], False) or SEVENS.get(combined_edge[-7:], False):
             return False       

        for pos in range(0,len(e)):

            p_ = FIVES.get(combined_edge[:-pos][-5:], False)   # pentameric unit exists
            s_ = SEVENS.get(combined_edge[:-pos][-7:], False)  # septameric unit exists

            if p_ or s_: 
                return False
    return True


def gen_string(strand_length, blueprint, complement_desired, front_edges=[], complement_front_edges=[], back_edges=[],
               complement_back_edges=[]):
    """
    Generates a new string of size n that doesn't intersect with existing strings. Those other strings are encoded in the FIVES

    :param strand_length: int
    :param blueprint: str
    :param complement_desired: bool
    :param front_edges: str[]
    :param complement_front_edges: str[]
    :param back_edges: str[]
    :param complement_back_edges: str[]

    :return: str
    """

    global FIVES
    global SEVENS
    global BACKTRACK_ARRAY
    global ADDED_UNITS

    blueprint = process_blueprint(strand_length, blueprint)
    blueprint_violation_array = get_blueprint_violation_array(blueprint)
    new_strand = ""
    
    attempt = 1

    while attempt < 30:
        while len(new_strand) < strand_length:
 
            next_possible_base = get_next_base(blueprint, blueprint_violation_array, new_strand, complement_desired,
                                               front_edges, complement_front_edges)

            # Case: add possible base (CONTINUE)
            if next_possible_base != "":
                new_strand += next_possible_base

            # Case: not possible base (BACKTRACK)
            else:
                backtrack_index = len(BACKTRACK_ARRAY) - 1

                if backtrack_index >= 0:

                    revert_units(new_strand, BACKTRACK_SEQ[backtrack_index], ADDED_UNITS)
                    backtrack_bases = BACKTRACK_ARRAY[backtrack_index]
                    backtrack_unit = random.choice(backtrack_bases)

                    # Backtrack Sequence
                    new_strand = BACKTRACK_SEQ[backtrack_index]
                    curr_length = len(new_strand)

                    # Replace removed base with alternate base
                    prev_4 = new_strand[len(new_strand) - 4:] 
                    prev_6 = ''                   
                    if curr_length >= 6:
                        prev_6 = new_strand[len(new_strand) - 6:] 
                    new_strand += backtrack_unit

                    # Update added units
                    ADDED_UNITS[curr_length] = update_fives(prev_4 + backtrack_unit, complement_desired) + \
                                               update_sevens(prev_6 + backtrack_unit, complement_desired)

                    # Update backtrack variables
                    if len(backtrack_bases) > 1:
                        backtrack_bases.remove(backtrack_unit)
                        BACKTRACK_ARRAY[backtrack_index] = backtrack_bases
                    else:
                        BACKTRACK_SEQ.pop()
                        BACKTRACK_ARRAY.pop()
                else:
                    # print("stoped at: "+new_strand + "\n***ATTEMPT #"+str(attempt)+"***\n")
                    new_strand = ""
                    break

        if len(new_strand) == strand_length and (backedge_check(new_strand, back_edges, complement_back_edges)
                                                 or not 'o' in blueprint):
            return str(new_strand)

        revert_units(new_strand, "", ADDED_UNITS, complement_desired)             
        attempt += 1

    # print( "Could not generate a string in "+ str(attempt) + " attempts." )
    return ""
  

#####################

def revert_units(curr_seq, rev_seq, added_bases):
    """
    Finds the 5 and 7 units added to reverted sequence to form current sequence. These units (which are keys in global
    variables FIVES and SEVENS are then reverted back to having FALSE as their value (they were switched to TRUE when
    added to current_sequence).

    :param curr_seq: str
    :param rev_seq: str
    :param added_bases: str[]
    :return:
    """
    global FIVES
    global SEVENS

    curr_len = len(curr_seq)
    rev_len = len(rev_seq)

    # print("backtrack to sequence:\n# "+rev_seq+"\n")
    for curr_index in range(rev_len, curr_len):
        for unit in added_bases[curr_index]:

            # print("remove: "+unit)
            if len(unit) == 5:
                FIVES[unit] = False
            else:
                SEVENS[unit] = False


def update_fives(new_unit, complement_desired):
    """
    Updates the FIVES data structure with the new string
    :param new_unit:
    :param complement_desired:
    :return:
    """
    global FIVES
    unit_list = []
    if len(new_unit) == 5:
        if complement_desired and (not FIVES.get(new_unit, True)): # not exist = True, added = True, not added = False
            FIVES[new_unit] = True
            unit_list.append(new_unit)

        new_reverse_unit = util.reverse_complement(new_unit)
        if not FIVES.get(new_reverse_unit, True):
            FIVES[new_reverse_unit] = True
            unit_list.append(new_reverse_unit)

    return unit_list


def update_sevens(new_unit, complement_desired):

    """
    Updates the SEVENS data structure with the new string
    :param new_unit:
    :param complement_desired:
    :return:
    """

    global SEVENS
    unit_list = []
    if len(new_unit) == 7:
        if complement_desired and (not SEVENS.get(new_unit, True)):
            unit_list.extend(update_sevens_approx(new_unit))

        new_reverse_unit = util.reverse_complement(new_unit) 
        if not SEVENS.get(new_reverse_unit, True):
            unit_list.extend(update_sevens_approx(new_reverse_unit))
    return unit_list


def update_sevens_approx(new_unit):
    """
    Updates the SEVENS aprox data structure with the new string
    :param new_unit:
    :return:
    """

    global SEVENS
    unit_list = []
    lets = "ACGT"
    i = 1
    while i < 6:
        for let in lets:
            new = new_unit[:i] + let + new_unit[i + 1:]
            if not SEVENS.get(new, True):
                SEVENS[new] = True
                unit_list.append(new)
        i += 1
    return unit_list


if __name__ == '__main__':

    gen_fives()
    gen_sevens()

    # test1
    test_str = "GGoATooooAAAAooooAoAToGGoGGGoGooooATGCoo"

    print("blueprint violation array: " + str(get_blueprint_violation_array(test_str)))
    print("total FIVES: "+ str(len(FIVES)))
    test2 = gen_string(40, test_str, True)


    # test_str_2 = "oCCCooooooooooooooATTooGGGoooGGGooooooooooo"
    # print("\n\nlength of FIVES used: "+ str(sum(FIVES.values())))
    # print ("total FIVES: "+ str(len(FIVES)))
    # test = gen_string(43, test_str_2, True)
    # print(test)
    # print("\n\nlength of FIVES used: "+ str(sum(FIVES.values())))

    # test2
    # print("\nresult: " + test2)
    # print("length of FIVES used: "+ str(sum(FIVES.values())))
    # print("length of SEVENS used: "+ str(sum(SEVENS.values())))

    # test3 = gen_string(100, "", True)
    # print("\nresult: " + test3)
    # print("length of FIVES used: "+ str(sum(FIVES.values())))
    # print("length of SEVENS used: "+ str(sum(SEVENS.values())))

    # test4 = gen_string(100, "", True)
    # print("\nresult: " + test4)
    # print("length of FIVES used: "+ str(sum(FIVES.values())))
    # print("length of SEVENS used: "+ str(sum(SEVENS.values())))

    # test5 = gen_string(100, "", True)
    # print("\nresult: " + test5)
    # print("length of FIVES used: "+ str(sum(FIVES.values())))
    # print("length of SEVENS used: "+ str(sum(SEVENS.values())))

    # print("\n___________\n1: " + test2)
    # print("2: " + test3)
    # print("3: " + test4)
    # print("4: " + test5)

    # #test3
    # for i in range(0, 10):
    #     #print gen_string(size_of_strand,"",True)
    #     #print("pentameric units remaining:" + str(FIVES.values().count(False)))
    #     #print("septameric units remaining:" + str(SEVENS.values().count(False)))
    #
