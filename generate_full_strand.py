from strand_utilities import strand_utilities
from nongeneric_gen_sequences import *
import pprint


util = strand_utilities()


def get_front_edges(strandname, component_list, full_strand_recipe):
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

            seq = ""
            for c in component_list:
                if c['name'] == comp_name:
                    seq = c['sequence']
                if c['name'] + "'" == comp_name:            
                    seq = util.reverse_complement(c['sequence'])    

            prev_edge = seq[-6:]


    return [edges, complement_edges]


def generate_strands(component_list, full_strand_recipe):
    """
    Generate strands brah

    Params:
        component_list (dict): dictionary of representing attributes of desired components 
            i.e. {'name':'a', 'length':15, 'blueprint':"", 'complement_desired':True}
        full_strand_recipe (list):  list of desired full strands 
            i.e.{ 'full strand name':['a','b','c'] , 'full strand name2': ['a','b'] }

    """
    attempts = 0

    while attempts < 100:
        
        fail = False
        success_count = 0

        for index in range(0, len(component_list)):
            component = component_list[index]
            front_edges = get_front_edges(component['name'], component_list, full_strand_recipe)

            seq = gen_string(component['length'], component['blueprint'], component['complement_desired'], front_edges[0], front_edges[1])
            
            if seq == '':
                if index == 0:
                    attempts +=1
                    fail = True
                    break

                index -= 1
                #print("backtrack to "+component_list[index-1]['name'])
                print("backtrack to index: "+str(index))
            else:
                component['sequence'] = seq

        if fail:
            print('fail at '+ component['name'])
            attempts += 1
            success_count = 0
            for component in component_list:
                component['sequence'] = ''
            random.shuffle(component_list)
        else:
            success_count += 1
            print('success count '+ str(success_count))
            if success_count == 2:
                break

    if attempts < 100:
        pprint.pprint(component_list)
        assembleFullStrands(component_list, full_strand_recipe)
    else:
        print("FAIL")
 

def assembleFullStrands(component_list, full_strand_recipe):
    fullstrands = {}

    for recipe in full_strand_recipe:
        seq = ''       
        name = ''
        for comp_name in recipe:
            name += comp_name +" "

            found = False
            for c in component_list:
                if c['name'] == comp_name:
                    seq += c['sequence']
                    found = True
                elif c['name'] + "'" == comp_name:            
                    seq += util.reverse_complement(c['sequence'])
                    found = True
            if not found:
                print(' COULD NOT FIND '+ comp_name) 


            fullstrands[name] = seq

    pprint.pprint(fullstrands)
    return fullstrands


#---------------------EXPERIMENT-----------------------------------------------------------------------------------------

e = {'name':'e', 'length':5, 'blueprint':"", 'complement_desired':True, 'sequence':''}
z = {'name':'z', 'length':5, 'blueprint':"", 'complement_desired':True, 'sequence':''}
y = {'name':'y', 'length':5, 'blueprint':"", 'complement_desired':True, 'sequence':''}
a = {'name':'a', 'length':5, 'blueprint':"", 'complement_desired':True, 'sequence':''}
b = {'name':'b', 'length':12, 'blueprint':"", 'complement_desired':True, 'sequence':''}
c = {'name':'c', 'length':8, 'blueprint':"", 'complement_desired':True, 'sequence':''}
d = {'name':'d', 'length':15, 'blueprint':"", 'complement_desired':True, 'sequence':''}
f = {'name':'f', 'length':5, 'blueprint':"", 'complement_desired':True, 'sequence':''}
w = {'name':'w', 'length':5, 'blueprint':"", 'complement_desired':True, 'sequence':''}
x = {'name':'x', 'length':5, 'blueprint':"", 'complement_desired':True, 'sequence':''}
g = {'name':'g', 'length':15, 'blueprint':"", 'complement_desired':True, 'sequence':''}

component_list = [y,z,e,a,b,c,d,f,w,x,g]

# full_strand_recipe = {}
# full_strand_recipe['1'] = ['a','b''c']
# full_strand_recipe['2'] = ["a'","b"] 
# full_strand_recipe['3'] = ['d','e','f',"a'"]

full_strand_recipe = [ ['b'] , ['d'] ]
full_strand_recipe += [ ["d'", "z'"] , ["c'", 'a'] , ["a'", 'c'] , ["d'", 'y'] , ["b'", "a'"] ]
full_strand_recipe += [ ["e'", "y'", 'd'] ]
full_strand_recipe += [ ['a', 'b', "e'", 'z', 'd'] , ["d'", 'y', 'e', "b'"] ]  
full_strand_recipe += [ ["g'", 'x', "f'","d'"] , ["y'", 'x','f', "d'"] ]

print(full_strand_recipe)

generate_strands(component_list, full_strand_recipe)
#print(get_front_edges('b', component_list, full_strand_recipe))


#-------------------------TEST-------------------------------------------------------------------------------------

"""
a = {'name':'a', 'length':15, 'blueprint':"", 'complement_desired':True, 'sequence':''}
b = {'name':'b', 'length':15, 'blueprint':"", 'complement_desired':True, 'sequence':''}
c = {'name':'c', 'length':15, 'blueprint':"", 'complement_desired':False, 'sequence':''}
d = {'name':'d', 'length':15, 'blueprint':"", 'complement_desired':False, 'sequence':''}
e = {'name':'e', 'length':15, 'blueprint':"", 'complement_desired':False, 'sequence':''}
f = {'name':'f', 'length':15, 'blueprint':"", 'complement_desired':False, 'sequence':''}

component_list = [b,c,d,e,f,a]
full_strand_recipe = [ ['a','b','c'], ["a'","b"], ['d','e','f',"a'"] ] 
print(full_strand_recipe)

generate_strands(component_list, full_strand_recipe)
#print(get_front_edges('b', component_list, full_strand_recipe))


"""







