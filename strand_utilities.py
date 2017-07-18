import random


def reverse_complement(s):
    """
    get reverse complement of given sequence
    :param s:
    :return:
    """

    rev = s[::-1]
    return complement(rev)


def complement(s):
    """
    get complement of given sequence
    :param s:
    :return:
    """

    out = ''
    for let in s:
        if let is 'A':
            out += 'T'
        elif let is 'C':
            out += 'G'
        elif let is 'G':
            out += 'C'
        elif let is 'T':
            out += 'A'
        else:
            out += let
    return out


def random_base():
    """
    get one random base
    :return: char
    """

    baselist = ['A', 'T', 'C', 'G']
    return random.choice(baselist)


def random_sequence(length):
    """
    get random sequence with given length
    :param length:
    :return: str
    """

    seq = ""
    for i in range(1, length):
        seq += random_base()
    return seq


def palindrome_check(seq):
    """
    check if given sequence is palindrome
    :param seq:
    :return: bool
    """

    # print("palindrome: " + seq + "  = "+ str((seq == reverse_complement(seq))))
    if "o" in seq:
        return False
    return seq is reverse_complement(seq)


def get_new_restriction_score(seq, new_base, restricted_seqs):
    """
    based on a given blueprint, it returns an array of possible sequences that dont violate or worsen
    number of restricted sequences.

    if palindromes are only found, it returns array with palindromes
    :param seq:
    :param new_base:
    :param restricted_seqs:
    :return: int
    """

    score = 0
    new_seq = seq + new_base
    seq_len = len(new_seq)

    # restricted seq check
    for res in restricted_seqs:
        res_length = len(res)
        if seq_len >= res_length:
            if new_seq[seq_len-res_length:] == res:
                score += 1

    # palindrome check
    if seq_len >= 4:
        six = new_seq[-4:]
        if palindrome_check(six):
            score += 2
            if six in seq:
                score += 3

    return score

