from strand_utilities import strand_utilities
util = strand_utilities()

full_strand_recipe_list = [] # [ ['a','b','c'] , ['a','b'] ] -> ab'c , ab 
component_list = {}

def get_front_edges(strandname):
    global full_strand_recipe_list
    edges = []
    complement_edges = []
    for recipe in full_strand_recipe_list:
        prev_edge = ""

        for comp_name in recipe:
            if prev_edge != "":
                if strandname == comp_name:
                    edges.append(prev_edge)
                if strandname + "'" == comp_name:
                    complement_edges.append(prev_edge)
            
            seq = component_list[comp_name]
            if len(seq) >= 4:
                prev_edge = seq[len(seq)-4 : len(seq)]
                if "'" in comp_name:            
                    prev_edge = util.reverse_complement(prev_edge)
            else:
                prev_edge = ""
    return [edges,complement_edges]

#test:
# full_strand_recipe_list = [['a','b','c','b']]
# component_list =  {'a':'ATCGAT','b':"",'c':"AAAA"}
# x = get_front_edges('b')
# print(str(x))
# print("\n.............\n")

def get_back_edges(strandname):
    global full_strand_recipe_list
    edges = []
    complement_edges = []
    for recipe in full_strand_recipe_list:
        prev_is_strand = False
        prev_is_complement_strand = False

        for comp_name in recipe:
            if strandname == comp_name:
                prev_is_strand = True
            if strandname + "'" == comp_name:
                prev_is_complement_strand = True
            
            seq = component_list[comp_name]
            if len(seq) >= 4:
                if prev_is_strand:
                    edges.append(seq[0 : 4])
                if prev_is_complement_strand:
                    complement_edges.append(util.reverse_complement(seq[0 : 4]))

    return [edges,complement_edges]

#test:
# full_strand_recipe_list = [['a','b','c','b','a']]
# component_list =  {'a':'ATCGAT','b':"",'c':"AAAA"}
# x = get_back_edges('b')
# print(str(x))

def generate_empty_component_list(name_list):
    component_list = {}
    for name in name_list:
        component_list[name] = ""
    return component_list

def generate_strands(name_list, full_strand_recipe_list):
    comp_list = generate_empty_component_list(name_list)
    index = 0
    while index < len(name_list):

    	index +=1





