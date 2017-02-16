# This software designs strings that don't have reverse complements
# of length five with other strings and that don't have groups of 7 
# nucleotides with 6 in common.

# For now, it makes no restriction about strings being complementary
# to themselves nor does it forbid two strings from being
# very similar. Please tell me what those constraints should be.



"""







"""

import random

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
            fives.append(i1+i2+i3+i4+i5)
            fives_present.append(False)

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
                sevens.append(i1+i2+i3+i4+i5+i6+i7)
                sevens_present.append(False)

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
def reversecomplement(s):
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

def genstring(strand_length):
  attempt = 1
  newstring = []
  
  while(attempt < 6) and (0 == len(newstring)):
    newstring = []
    for i in range(strand_length):
              
      #find pentameric unit that hasnt been added
      current_index = random.randint(0, len(fives) - 50)      
      while (fives_present[current_index] or sevens_present[current_index]) and (current_index < len(fives)): 
        current_index+=1

      #Case: Found pentameric unit
      if (current_index < len(fives)): 
        newstring = fives[current_index]
       
        k = 5
        while (k < strand_length):       #while k < strand length 
          
          s = newstring[len(newstring) - 4:]     #get previous four bases for( _ _ _ _ + new base )
          if (k >= 7):
            t = newstring[len(newstring) - 6:]   #get previous six bases for( _ _ _ _ _ + new base )
          else:
            t = ''
            
          goodones = []
          for base in 'ACGT':
            pentameric_unit = s+base        # test pentameric unit
            septameric_unit = t+base        # test septameric unit
            
            m = fives.index(pentameric_unit)  #find index of pentameric unit
            if (not fives_present[m]) and ((7 > len(septameric_unit)) or (not sevensapprox(septameric_unit))) : #if pentameric unit is not added yet && seven unit does not exist
              goodones.append(base)

          if (0 < len(goodones)):
            chosen_unit = random.choice(goodones)
            newstring += chosen_unit  #ADD 
          else:
            k = strand_length # have to stop
            newstring = []
          k+= 1

        newrev = reversecomplement(newstring)
        updatefives(newrev)
        updatesevens(newrev)
        allstrings.append(newstring)
        return newstring
      else:
        newstring = []
  return "Could not generate a string in ", attempt, " attempts."

  

# updates the fives data structure with the new string
def updatefives(newstring):
  global fives_present
  for i in range(len(newstring)-4):
    s = newstring[i:i+5]
    j = fives.index(s)
    # if fives_present[j]:
      # print "We have a fives problem at position ", i, " of newstring ", newstring, " with respect to ", fives[j]
      # print "The string in question is: ",  s
    fives_present[j] = True
    

# updates the sevens data structure with the new string
def updatesevens(newstring):
  global sevens_present
  for i in range(len(newstring)-6):
    s = newstring[i:i+7]
    j = sevens.index(s)
    if sevens_present[j]:
      print "We have a sevens problem at position ", i, " of newstring"
      print "Letters are: ",  s
    sevens_present[j] = True
    
  



# DATA
fives = []
fives_present = []
sevens = []
sevens_present = []
allstrings = []

n = 25 # size of strings

# EXECUTION
genfives()
gensevens()
print genstring(n)
print genstring(n)
print genstring(n)
print genstring(n)
print genstring(n)
print genstring(n)
print genstring(n)



x = [i for i in range(len(fives_present)) if  fives_present[i] ]
print "number of fives that are present: ", len(x)
x = [i for i in range(len(sevens_present)) if  sevens_present[i] ]
print "number of sevens that are present: ", len(x)
