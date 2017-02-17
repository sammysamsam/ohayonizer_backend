
"""
    
factors:
    - no polypurine (AAAA/TTTT) and  + CCC/GGG
    - no branch migration in junctions
    - no symmetry
    - at least 5 nucleotides between junctions (per arm)
    - avoid alternating CGs and other purines and pyrmadines (3-4 )
    - having CG on the ends of strands is optimal to prevent fraying

approach:
    - using unique units (eliminates symmetry)
    - junction flanking pairs (eliminates branch migration)




user input:

list of desired strand components + length + "blueprint/partially filled in sequence"

list of full strands ( linked components + marked points of junction???)

percentage of CG/AT per component?? -> melting point per component???



big picture: 

1. strand components + blueprint 
2. full strands (no junctions) 
3. full strands + junctions 
4. melting point per full strand or component ???

"""

print("");
