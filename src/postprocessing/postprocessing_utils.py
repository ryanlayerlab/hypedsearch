from preprocessing import clustering


def make_db_mapping(sqllite_database):
    name_list = []
    proteins = sqllite_database.query_proteins()
    [name_list.append(x) for x in proteins]

    num_to_prot_name = dict()
    for i, x in enumerate(name_list):
        num_to_prot_name[i] = name_list[i]
    
    return num_to_prot_name

def make_db_mapping_by_key(db):
    name_list = []
    [name_list.append(x) for x in db.proteins]

    num_to_prot_name = dict()
    for i, x in enumerate(name_list):
        prot_description = x[0]
        prot_arr = prot_description.split("|")
        prot_arr2 = prot_arr[2].split(" ")
        prot_name=prot_arr2[0]
        num_to_prot_name[prot_name] = name_list[i]
    
    return num_to_prot_name

def label_alignments(alignment):
    type = alignment[3]
    if type:
        label = "Hybrid"
    else:
        label = "Native"
    return label

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
        total_sequence = b_seq + '-' + y_seq
    else:
        b_seq = alignment_info[0][2][0][6]
        total_sequence = b_seq

    return total_sequence

def find_parent_proteins(alignment_info, sqllite_database):
    left_parents, right_parents = set(), set()
    for merge in alignment_info:
        left_num, right_num = merge[2][0][5], merge[2][1][5]
        get_name = lambda x: x.split('|')[-1].split()[0]
        left_parent = sqllite_database.get_protein(left_num)[2]
        right_parent = sqllite_database.get_protein(right_num)[2]
        left_parents.add(left_parent)
        right_parents.add(right_parent)
    return left_parents, right_parents

def get_extensions(alignment_info, sqllite_database):
    left_extensions, right_extensions = [], []
    for merge in alignment_info:
        b_side, y_side = merge[2][0], merge[2][1]
        b_parent_seq = sqllite_database.get_protein(b_side[5])[2]
        y_parent_seq = sqllite_database.get_protein(y_side[5])[2]
        b_extension = b_parent_seq[max(b_side[1]-1, 0):b_side[1]]
        y_extension = y_parent_seq[y_side[2]:min(y_side[2]+1,len(y_parent_seq))]
        left_extensions.append(b_extension) if b_extension != '' else '-'
        right_extensions.append(y_extension) if y_extension != '' else '-'
        
    return left_extensions, right_extensions

def postprocessing(rescored_alignments, sqllite_database, spectrum, number_hybrids, number_natives):
    postprocessed_alignments = []
    sorted_alignments = sorted(rescored_alignments, key = lambda x: (x[0], x[1]), reverse=True)
    i = 0
    if (len(sorted_alignments) > 0):
        max = number_hybrids if label_alignments(sorted_alignments[0]) == "Hybrid" else number_natives
        
        while i < max and i < len(sorted_alignments):
            alignment = sorted_alignments[i]
            alignment_info = rescored_alignments[alignment]
            label = label_alignments(alignment)
            left_proteins, right_proteins = find_parent_proteins(alignment_info, sqllite_database)
            sequence = get_sequence(alignment_info, label)
            extended_sequence = get_extensions(alignment_info, sqllite_database)
            b_scores, y_scores = get_scores(alignment_info)
            total_score = alignment[0]
            gaussian_score = alignment[1]
            mass_error_sum = rescored_alignments[alignment][0][3]

            postprocessed_alignments.append((label, left_proteins, right_proteins, sequence, b_scores, y_scores, total_score/len(sequence.replace("-","")), gaussian_score, extended_sequence, spectrum.precursor_mass, spectrum.precursor_charge, mass_error_sum))
            i += 1
    return postprocessed_alignments