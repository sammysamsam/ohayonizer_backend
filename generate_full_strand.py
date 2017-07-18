from strand_utilities import strand_utilities
from gen_sequences import *
import pprint
from itertools import permutations

util = strand_utilities()


def get_front_edges(strandname, component_list, full_strand_recipe, n_size):
    edges = []
    complement_edges = []
    for recipe in full_strand_recipe:
        prev_edge = ""

        for comp_name in recipe:
            if prev_edge != "":
                if strandname == comp_name:
                    edges.append(prev_edge)
                if strandname + "'" == comp_name:
                    complement_edges.append(prev_edge)

            for c in component_list:
                if len(c['sequence']) > n_size:
                    if c['name'] == comp_name:
                        prev_edge = c['sequence'][-n_size:]
                        break
                    if c['name'] + "'" == comp_name:            
                        prev_edge = util.reverse_complement(c['sequence'])[-n_size:]    
                        break

    return [edges, complement_edges]

def get_back_edges(strandname, component_list, full_strand_recipe, n_size):
    edges = []
    complement_edges = []
    for recipe in full_strand_recipe:
        prev_is_strand = False
        prev_is_complement_strand = False

        for comp_name in recipe:
            if strandname == comp_name:
                prev_is_strand = True
            if strandname + "'" == comp_name:
                prev_is_complement_strand = True
            
            seq = ""
            for c in component_list:
                if len(seq) >= n_size:
                    if c['name'] == comp_name:
                        edges.append(c['sequence'][0 : n_size])
                        break      
                    if c['name'] + "'" == comp_name:            
                        seq = util.reverse_complement(c['sequence']) 
                        complement_edges.append(seq[0 : n_size])
                        break
    return [edges,complement_edges]
 
def get_nmers(size, component_list, full_strand_recipe):

    ''' 
    Returns the dictionary of n+1mers using the nmer list passed in. If contains restricted sequence, value = True
    
    '''

    #step 1: get full nmer list (before strand nmers are added)
    letters = (' A C G T'*size).split()
    nmers = {}
    for new_nmer in [''.join(s) for s in permutations(letters,size)]:
        #check if nmer unit can base pair with self or contains repeated units
        if (util.palindrome_check(new_nmer)) or new_nmer in set(['AAAA', 'TTTT', 'CCC', 'GGG'])::
            nmers[new_nmer] = True
        else:
            nmers[new_nmer] = False  
    
    #step 2: add full strand nmers
    for full_strand in full_strand_recipe:


    return nmers


def generate_strands(component_list, full_strand_recipe):
    """
    Generate strands brah

    Params:
        component_list (dict): dictionary of representing attributes of desired components 
            i.e. {'name':'a', 'length':15, 'blueprint':"", 'complement':True}
        full_strand_recipe (list):  list of desired full strands 
            i.e.{ 'full strand name':['a','b','c'] , 'full strand name2': ['a','b'] }

    """

    attempts = 0
    while attempts < 2000:
        fail = False
        index = 0
        while index < len(component_list):
            component = component_list[index]
            
            for nmer_size in range(4,8):
                front_edges = get_front_edges(nmer_size, component['name'], component_list, full_strand_recipe)
                back_edges = get_back_edges(nmer_size, component['name'], component_list, full_strand_recipe)
                current_nmers = get_nmers(nmer_size)
                current_mmers = get_nmers(nmer_size+2)
                
                seq = gen_string(component['length'], component['blueprint'], component['complement'], front_edges[0], front_edges[1], back_edges[0], back_edges[1])
            
                if seq == '':
                    if index == 0:
                        fail = True
                        break
                    index -= 1
                    print("backtrack to "+component_list[index-1]['name'])
                    print("backtrack to index: "+str(index))
                else:
                    component['sequence'] = seq
                    index += 1

        if fail:
            print('fail at '+ component['name'])
            attempts +=1

            for component in component_list:
                component['sequence'] = ''

            fives = {}
            sevens = {}
            genfives()
            gensevens()

            random.shuffle(component_list)
        
        else:
            return [component_list, assembleFullStrands(component_list, full_strand_recipe) ]

    print("FAIL \n\n")



def assembleFullStrands(component_list, full_strand_recipe):
    fullstrands = {}

    for recipe in full_strand_recipe:
        seq = '' 
        name = ''

        for comp_name in recipe:
            name += comp_name + " "

            found = False
            for c in component_list:
                if c['name'] == comp_name:
                    seq += c['sequence']
                    found = True   
                    break 
                elif c['name'] + "'" == comp_name:            
                    seq += util.reverse_complement(c['sequence'])
                    found = True
                    break
            if not found:
                print(' COULD NOT FIND '+ comp_name) 

        fullstrands[name] = seq

    #pprint.pprint(fullstrands)
    return fullstrands

#C=8; A/E/Z/Y = 5; B/D = 15; S = 6 or 7

#---------------------EXPERIMENT-----------------------------------------------------------------------------------------

# d = {'name':'d', 'length':15, 'blueprint':"ATCAAGCATAGTATC", 'complement':True, 'sequence':''}
# g = {'name':'g', 'length':15, 'blueprint':"TCACAGAACCATCCA", 'complement':True, 'sequence':''}
# e = {'name':'e', 'length':5, 'blueprint':"GTAAG", 'complement':True, 'sequence':''}
# z = {'name':'z', 'length':5, 'blueprint':"AGGTT", 'complement':True, 'sequence':''}
# y = {'name':'y', 'length':5, 'blueprint':"TATGA", 'complement':True, 'sequence':''}
# a = {'name':'a', 'length':5, 'blueprint':"TCCAG", 'complement':True, 'sequence':''}
# b = {'name':'b', 'length':12, 'blueprint':"AATCTAATAACC", 'complement':True, 'sequence':''}
# c = {'name':'c', 'length':8, 'blueprint':"GCGTGGTG", 'complement':True, 'sequence':''}
# f = {'name':'f', 'length':5, 'blueprint':"GCTTC", 'complement':True, 'sequence':''}
# x = {'name':'x', 'length':5, 'blueprint':"CCTTA", 'complement':True, 'sequence':''}
# w = {'name':'w', 'length':5, 'blueprint':"TAACA", 'complement':True, 'sequence':''}
# j = {'name':'j', 'length':8, 'blueprint':"TGTGTTCG", 'complement':True, 'sequence':''}
# s = {'name':'s', 'length':8, 'blueprint':"", 'complement':True, 'sequence':''}

# component_list = [y,b,e,w,d,g,z,a,c,f,x,j,s]

# full_strand_recipe = [ ['b'] , ['d'] ]
# full_strand_recipe += [ ["d'", "z'"] , ["c'", 'a'] , ["a'", "c"] , ["d'", 'y'] , ["b'", "a'"] , ["j","y"], ["j'","y'"]]
# full_strand_recipe += [ ["e'", "y'", 'd'] ]
# full_strand_recipe += [ ['a', 'b', "e'", 'z', 'd', "s'"] , ['s', "d'", 'y', 'e', "b'"] ]  
# full_strand_recipe += [ ["g'", 'x', "f","d'"] , ["y'", 'd','f', "w", 'g'] ]
# print (full_strand_recipe)
# #print(full_strand_recipe) D, d’z’, c’a, a’c, abe’zd, d’yeb’, b’a’, b, e’y’d, d’y, y’dfwg, g’xf’d’

# results = generate_strands(component_list, full_strand_recipe)
# component_list = results[0]
# full_strand_list = results[1]

# for x in component_list:
#     del x['blueprint']
#     del x['length']


# pprint.pprint(component_list)

# count = 0
# print('{"CL":[')
# for temp in full_strand_list:
#     count +=1
#     print('{"name": "' + temp.replace("'","`") + '", "sequence": "'+full_strand_list[temp] + '", "complement": "false"},')
# print('],"FSL":[]}')

# d = {'name':'d', 'length':15, 'blueprint':"", 'complement':True, 'sequence':''}
# g = {'name':'g', 'length':15, 'blueprint':"", 'complement':True, 'sequence':''}
# e = {'name':'e', 'length':5, 'blueprint':"", 'complement':True, 'sequence':''}
# z = {'name':'z', 'length':5, 'blueprint':"", 'complement':True, 'sequence':''}
# y = {'name':'y', 'length':5, 'blueprint':"", 'complement':True, 'sequence':''}
# a = {'name':'a', 'length':5, 'blueprint':"", 'complement':True, 'sequence':''}
# b = {'name':'b', 'length':12, 'blueprint':"", 'complement':True, 'sequence':''}
# c = {'name':'c', 'length':8, 'blueprint':"", 'complement':True, 'sequence':''}
# f = {'name':'f', 'length':5, 'blueprint':"", 'complement':True, 'sequence':''}
# x = {'name':'x', 'length':5, 'blueprint':"", 'complement':True, 'sequence':''}
# w = {'name':'w', 'length':5, 'blueprint':"", 'complement':True, 'sequence':''}
# j = {'name':'j', 'length':8, 'blueprint':"", 'complement':True, 'sequence':''}
# s = {'name':'s', 'length':7, 'blueprint':"", 'complement':True, 'sequence':''}

# component_list = [y,b,e,w,d,g,z,a,c,f,x,j,s]

# full_strand_recipe = [ ['b'] , ['d'] ]
# full_strand_recipe += [ ["d'", "z'"] , ["c'", 'a'] , ["a'", "c"] , ["d'", 'y'] , ["b'", "a'"] , ["j","y"], ["j'","y'"]]
# full_strand_recipe += [ ["e'", "y'", 'd'] ]
# full_strand_recipe += [ ['a', 'b', "e'", 'z', 'd', "s'"] , ['s', "d'", 'y', 'e', "b'"] ]  
# full_strand_recipe += [ ["g'", 'x', "f'","d'"] , ["y'", 'd','f', "w", 'g'] ]
# print (full_strand_recipe)
# #print(full_strand_recipe) D, d’z’, c’a, a’c, abe’zd, d’yeb’, b’a’, b, e’y’d, d’y, y’dfwg, g’xf’d’

# results = generate_strands(component_list, full_strand_recipe)
# component_list = results[0]
# full_strand_list = results[1]

# for x in component_list:
#     del x['blueprint']
#     del x['length']


# pprint.pprint(component_list)

# count = 0
# print('{"CL":[')
# for temp in full_strand_list:
#     count +=1
#     print('{"name": "' + temp.replace("'","`") + '", "sequence": "'+full_strand_list[temp] + '", "complement": "false"},')
# print('],"FSL":[]}')


=======
# d = {'name':'d', 'length':15, 'blueprint':"ATCAAGCATAGTATC", 'complement':True, 'sequence':''}
# g = {'name':'g', 'length':15, 'blueprint':"TCACAGAACCATCCA", 'complement':True, 'sequence':''}
# e = {'name':'e', 'length':5, 'blueprint':"GTAAG", 'complement':True, 'sequence':''}
# z = {'name':'z', 'length':5, 'blueprint':"AGGTT", 'complement':True, 'sequence':''}
# y = {'name':'y', 'length':5, 'blueprint':"TATGA", 'complement':True, 'sequence':''}
# a = {'name':'a', 'length':5, 'blueprint':"TCCAG", 'complement':True, 'sequence':''}
# b = {'name':'b', 'length':12, 'blueprint':"AATCTAATAACC", 'complement':True, 'sequence':''}
# c = {'name':'c', 'length':8, 'blueprint':"GCGTGGTG", 'complement':True, 'sequence':''}
# f = {'name':'f', 'length':5, 'blueprint':"GCTTC", 'complement':True, 'sequence':''}
# x = {'name':'x', 'length':5, 'blueprint':"CCTTA", 'complement':True, 'sequence':''}
# w = {'name':'w', 'length':5, 'blueprint':"TAACA", 'complement':True, 'sequence':''}
# j = {'name':'j', 'length':8, 'blueprint':"TGTGTTCG", 'complement':True, 'sequence':''}
# s = {'name':'s', 'length':7, 'blueprint':"", 'complement':True, 'sequence':''}

# component_list = [y,b,e,w,d,g,z,a,c,f,x,j]

# full_strand_recipe = [ ['b'] , ['d'] ]
# full_strand_recipe += [ ["d'", "z'"] , ["c'", 'a'] , ["a'", "c"] , ["d'", 'y'] , ["b'", "a'"] , ["j","y"], ["j'","y'"]]
# full_strand_recipe += [ ["e'", "y'", 'd'] ]
# full_strand_recipe += [ ['a', 'b', "e'", 'z', 'd', "s'"] , ['s', "d'", 'y', 'e', "b'"] ]  
# full_strand_recipe += [ ["g'", 'x', "f","d'"] , ["y'", 'd','f', "w", 'g'] ]
# print (full_strand_recipe)
# #print(full_strand_recipe) D, d’z’, c’a, a’c, abe’zd, d’yeb’, b’a’, b, e’y’d, d’y, y’dfwg, g’xf’d’

# results = generate_strands(component_list, full_strand_recipe)
# component_list = results[0]
# full_strand_list = results[1]

# for x in component_list:
#     del x['blueprint']
#     del x['length']


# pprint.pprint(component_list)

# count = 0
# print('{"CL":[')
# for temp in full_strand_list:
#     count +=1
#     print('{"name": "' + temp.replace("'","`") + '", "sequence": "'+full_strand_list[temp] + '", "complement": "false"},')
# print('],"FSL":[]}')


#-------------------------TEST-------------------------------------------------------------------------------------

# a = {'name':'a', 'length':25, 'blueprint':"", 'complement':True, 'sequence':''}
# b = {'name':'b', 'length':25, 'blueprint':"", 'complement':True, 'sequence':''}
# c = {'name':'c', 'length':15, 'blueprint':"", 'complement':False, 'sequence':''}
# d = {'name':'d', 'length':15, 'blueprint':"", 'complement':False, 'sequence':''}
# e = {'name':'e', 'length':15, 'blueprint':"", 'complement':False, 'sequence':''}
# f = {'name':'f', 'length':15, 'blueprint':"", 'complement':False, 'sequence':''}

# component_list = [b,c,d,e,f,a]
# full_strand_recipe = [ ['a','b','c'], ["a'","b"], ['d','e','f',"a'"] ] 
# print(full_strand_recipe)

# results = generate_strands(component_list, full_strand_recipe)
# component_list = results[0]
# full_strand_list = results[1]
# #print(get_front_edges('b', component_list, full_strand_recipe))
# pprint.pprint(component_list)

# count = 0
# print('{"CL":[')
# for temp in full_strand_list:
#     count +=1
#     print('{"name": "' + temp.replace("'","`") + '", "sequence": "'+full_strand_list[temp] + '", "complement": "false"},')
# print('],"FSL":[]}')







