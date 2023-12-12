## purpose of this function is unclear. 
## creates a list and then return a dict with indices mapping to elements
## couldn't we just return the list? the dict is redundant
def make_db_mapping(db): 
    name_list = [x for x in db.proteins]

    #edit this after testing
    num_to_prot_name = dict()
    for i, name in enumerate(name_list):
        num_to_prot_name[i] = name_list[i]
    
    return num_to_prot_name

def make_db_mapping_by_key(db):
    name_list = [x for x in db.proteins]

    # edit this after testing
    num_to_prot_name = dict()
    for i, x in enumerate(name_list):
        prot_description = x[0]
        prot_arr = prot_description.split("|")
        prot_arr2 = prot_arr[2].split(" ")
        prot_name=prot_arr2[0]
        num_to_prot_name[prot_name] = name_list[i]
    
    return num_to_prot_name

def label_alignments(alignment):
    #b_mass, b_start, b_end, 1, b_charge, b_pid, b_score
    type = alignment[3]
    return "Hybrid" if type else "Native"

def get_scores(alignment_info):
    b_scores, y_scores = set(), set()
    for merge in alignment_info:
        bs, ys = merge[2][0][7], merge[2][1][7]
        b_scores.add(bs)
        y_scores.add(ys)
        
    return b_scores, y_scores

def get_precursor_dist(alignment):
    return alignment[1]

def get_sequence(alignment_info, label):
    if label == "Hybrid":
        b_seq, y_seq = alignment_info[0][2][0][6], alignment_info[0][2][1][6]
        total_sequence = f'{b_seq}-{y_seq}'
    else:
        b_seq = alignment_info[0][2][0][6]
        total_sequence = b_seq

    return total_sequence

def find_parent_proteins(alignment_info, db_mapping):
    get_name = lambda x: x.split('|')[-1].split()[0]
    left_parents, right_parents = set(), set()
    for merge in alignment_info:
        left_num, right_num = merge[2][0][5], merge[2][1][5]
        left_parent, right_parent = get_name(db_mapping[left_num].description), get_name(db_mapping[right_num].description)
        left_parents.add(left_parent)
        right_parents.add(right_parent)
    return left_parents, right_parents

def get_extensions(alignment_info, protein_list):
    left_extensions, right_extensions = [], []
    for merge in alignment_info:
        b_side, y_side = merge[2][0], merge[2][1]
        b_parent_seq, y_parent_seq = protein_list[b_side[5]][1], protein_list[y_side[5]][1]
        b_extension = b_parent_seq[max(b_side[1]-1, 0):b_side[1]]
        y_extension = y_parent_seq[y_side[2]:min(y_side[2]+1,len(y_parent_seq))]
        left_extensions.append(b_extension) if b_extension != '' else '-'
        right_extensions.append(y_extension) if y_extension != '' else '-'
        
    return left_extensions, right_extensions

# spec_num, non_hybrid, proteins, sequence, b_score, y_score, total_score, precursor_distance, total_mass_error
def postprocessing(alignments, db, input_spectrum, num_hybrids, num_natives):
    postprocessed_alignments = []
    db_mapping = make_db_mapping(db)
    sorted_alignments = sorted(alignments, key = lambda x: (x[0], x[1]), reverse=True)
    i = 0
    if (len(sorted_alignments) > 0):
        max = num_hybrids if label_alignments(sorted_alignments[0]) == "Hybrid" else num_natives
        
        while i < len(sorted_alignments) and i < max:
            alignment = sorted_alignments[i]
            alignment_info = alignments[alignment]
            label = label_alignments(alignment)
            left_proteins, right_proteins = find_parent_proteins(alignment_info, db_mapping)
            sequence = get_sequence(alignment_info, label)
            extended_sequence = get_extensions(alignment_info, db.proteins)
            b_scores, y_scores = get_scores(alignment_info)
            total_score = alignment[0]
            gaussian_score = alignment[1]
            mass_error_sum = alignments[alignment][0][3]

            postprocessed_alignments.append((label, left_proteins, right_proteins, sequence, b_scores, y_scores, total_score/len(sequence.replace("-","")), gaussian_score, extended_sequence, input_spectrum.precursor_mass, input_spectrum.precursor_charge, mass_error_sum))
            i += 1
    return postprocessed_alignments