from preprocessing import clustering


def make_db_mapping(db):
    name_list = []
    [name_list.append(x) for x in db.proteins]

    num_to_prot_name = dict()
    for i, x in enumerate(name_list):
        num_to_prot_name[i] = name_list[i]
    
    return num_to_prot_name

def label_alignments(alignment):
    #b_mass, b_start, b_end, 1, b_charge, b_pid, b_score
    type = alignment[3]
    if type:
        label = "Hybrid"
    else:
        label = "Natural"
    return label

def get_score(alignment):
    b_c, y_c = alignment[2][0], alignment[2][1]
    return b_c[6], y_c[6]

def get_precursor_dist(alignment):
    return alignment[1]

def get_sequence(alignment, proteins):
    hybrid = alignment[3]
    b_c, y_c = alignment[2][0], alignment[2][1]
    b_pid, y_pid = b_c[5], y_c[5]
    b_start, y_start = b_c[1], y_c[1]
    b_end, y_end = b_c[2], y_c[2]
    if hybrid:
        total_sequence = clustering.find_sequence(b_pid, b_start, b_end, proteins) + '-' + clustering.find_sequence(y_pid, y_start, y_end, proteins)
    else:
        total_sequence = clustering.find_sequence(b_pid, b_start, y_end, proteins)

    return total_sequence

def find_parent_protein(alignment, db_mapping):
    b_c, y_c = alignment[2][0], alignment[2][1]
    left_num, right_num = b_c[5], y_c[5]
    get_name = lambda x: x.split('|')[-1].split()[0]
    left_parent, right_parent = get_name(db_mapping[left_num].description), get_name(db_mapping[right_num].description)
    return left_parent, right_parent

def get_extended_sequence(alignment, protein_list, label, sequence):
    if label != "Hybrid":
        return sequence
    else:
        b_start, y_end = alignment[2][0][1], alignment[2][1][2]
        b_parent_id, y_parent_id = alignment[2][0][5], alignment[2][1][5]
        b_parent_seq, y_parent_seq = protein_list[b_parent_id][1], protein_list[y_parent_id][1]
        b_extensions = b_parent_seq[max(b_start-5, 0):b_start]
        y_extensions = y_parent_seq[y_end:min(y_end+5,len(y_parent_seq))]
        return b_extensions + sequence + y_extensions

# spec_num, non_hybrid, proteins, sequence, b_score, y_score, total_score, precursor_distance, total_mass_error
def postprocessing(alignments, db):
    postprocessed_alignments = []
    db_mapping = make_db_mapping(db)
    for alignment in alignments[:min(10, len(alignments))]:
        label = label_alignments(alignment)
        left_protein, right_protein = find_parent_protein(alignment, db_mapping)
        sequence = get_sequence(alignment, db.proteins)
        extended_sequence = get_extended_sequence(alignment, db.proteins, label, sequence)
        b_score, y_score = get_score(alignment)
        total_score = alignment[0]
        precursor_distance = get_precursor_dist(alignment)
        postprocessed_alignments.append((label, left_protein, right_protein, sequence, b_score, y_score, total_score, precursor_distance, extended_sequence, alignment))
    return postprocessed_alignments