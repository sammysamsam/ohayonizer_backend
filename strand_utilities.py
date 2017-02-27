import random
class strand_utilities:

	# take reverse complement of DNA
	def reverse_complement(e, s):
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

	#get one random base
	def random_base():
		baselist = ['A','T','C','G']
		return random.choice(baselist)


	#get random sequence with given length
	def random_sequence(length):
		seq = ""
		for i in range(1,length):
			seq += random_base();
		return seq

	"""
	Based on a given blueprint, it returns an array of possible sequences that dont violate or worsen 
	number of restricted sequences.

	if palindromes are only found, it returns array with palindromes

	"""
	def get_pentameric_possiblilities(e, blueprint):
		restricted_seq = ['AAAA', 'TTTT', 'CCC', 'GGG','CCCC','GGGG','TTTTT','AAAAA']
		init_score = 0
		for res in restricted_seq:
			if res in blueprint:
				init_score += 1

		results = []
		palim_results = []
		if(len(blueprint) == 5):
			lets = []

			for i in range(0,5):
				if(blueprint[i] == 'o'):
					lets.append("ATCG")
				else:
					lets.append(blueprint[i])

			for i1 in lets[0]:
				for i2 in lets[1]:
					for i3 in lets[2]:
						for i4 in lets[3]:
							for i5 in lets[4]:
								new_five = i1 + i2 + i3 + i4 + i5
								new_score = 0
								for sequence in restricted_seq:
									if sequence in new_five:
										new_score+=1
								
								if (new_score == init_score):
									if (new_five == new_five[::-1]):
										palim_results.append(new_five)
									else:
										results.append(new_five)
			if(len(results) != 0):
				return results
			else:
				return palim_results

"""

#TEST AREA 

util = strand_utilities()

x = util.get_pentameric_possiblilities("GGGGG")
print(x)



"""

