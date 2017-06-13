

        first_bases_len = 6
        if strand_length < 6:
            first_bases_len = strand_length

        for num in range(0 , first_bases_len):
            good_base_count = {'A':0, 'T':0, 'C':0, 'G':0}

            for e in front_edges:
                combined_edge = e + new_strand

                prev_4 = combined_edge[-4:]
                prev_6 = combined_edge[-6:]

                bases = get_next_base(prev_4, prev_6, blueprint, blueprint_violation_array, combined_edge, len(new_strand), complement_desired)
                
                for base in bases:
                    good_base_count[base] += 1
            
            for e in complement_front_edges:
                combined_edge = e + util.reverse_complement(new_strand)

                prev_4 = combined_edge[-4:]
                prev_6 = combined_edge[-6:]

                bases = get_next_base(prev_4, prev_6, blueprint, blueprint_violation_array, combined_edge, len(new_strand), complement_desired)
                
                for base in bases:
                    good_base_count[util.complement(base)] += 1


            base_choices = []
            for b in good_base_count:
                if good_base_count[b] == len(front_edges) + len(complement_front_edges):
                    base_choices.append(b)

            if len(base_choices) == 0:
                new_strand = ""
                break
            else:
                new_strand += random.choice(base_choices)



def get_first_six_bases(strand_length, blueprint, blueprint_violation_array, complement_desired, front_edges, complement_front_edges):
    new_strand = ""

    first_bases_len = 6
    if strand_length < 6:
        first_bases_len = strand_length

    for num in range(0 , first_bases_len):
        good_base_count = {'A':0, 'T':0, 'C':0, 'G':0}

        for e in front_edges:
            combined_edge = e + new_strand

            prev_4 = combined_edge[-4:]
            prev_6 = combined_edge[-6:]

            bases = get_next_base(prev_4, prev_6, blueprint, blueprint_violation_array, combined_edge, len(new_strand), complement_desired)
            
            for base in bases:
                good_base_count[base] += 1
        
        for e in complement_front_edges:
            combined_edge = e + util.reverse_complement(new_strand)

            prev_4 = combined_edge[-4:]
            prev_6 = combined_edge[-6:]

            bases = get_next_base(prev_4, prev_6, blueprint, blueprint_violation_array, combined_edge, len(new_strand), complement_desired)
            
            for base in bases:
                good_base_count[util.complement(base)] += 1


        base_choices = []
        for b in good_base_count:
            if good_base_count[b] == len(front_edges) + len(complement_front_edges):
                base_choices.append(b)

        if len(base_choices) == 0:
            return ""
        else:
            new_strand += random.choice(base_choices)
    return new_strand