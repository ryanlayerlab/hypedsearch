
def make_db_mapping(db):
    name_list = []
    [name_list.append(x) for x in db.proteins]

    num_to_prot_name = dict()
    for i, x in enumerate(name_list):
        num_to_prot_name[i] = name_list[i]
    
    return num_to_prot_name

def label_alignments(alignment):
    b_c, y_c = alignment[2][3], alignment[2][4]
    if b_c[2]==y_c[1]-1 or b_c[4] == y_c[4]:
        label = "Natural"
    else:
        label = "Hybrid"
    return label

def get_score(alignment):
    b_c, y_c = alignment[2][3], alignment[2][4]
    return b_c[3], y_c[3]

def get_precursor_dist(alignment):
    return 1/alignment[1]

def get_sequence(alignment):
    b_c, y_c = alignment[2][3], alignment[2][4]
    b_seq = b_c[4]
    y_seq = y_c[4]
    if b_seq == y_seq:
        total_seq = b_seq
    else:
        total_seq = b_seq +'-'+ y_seq
    return total_seq

def find_parent_protein(alignment, db_mapping):
    b_c, y_c = alignment[2][3], alignment[2][4]
    left_num, right_num = b_c[0], y_c[0]
    get_name = lambda x: x.split('|')[-1].split()[0]
    left_parent, right_parent = get_name(db_mapping[left_num].description), get_name(db_mapping[right_num].description)
    return left_parent, right_parent

# spec_num, non_hybrid, proteins, sequence, b_score, y_score, total_score, precursor_distance, total_mass_error
def postprocessing(alignments, db):
    postprocessed_alignments = []
    db_mapping = make_db_mapping(db)
    for alignment in alignments[:min(10, len(alignments))]:
        label = label_alignments(alignment)
        left_protein, right_protein = find_parent_protein(alignment, db_mapping)
        sequence = get_sequence(alignment)
        b_score, y_score = get_score(alignment)
        total_score = alignment[0]
        precursor_distance = get_precursor_dist(alignment)
        postprocessed_alignments.append((label, left_protein, right_protein, sequence, b_score, y_score, total_score, precursor_distance, alignment))
    return postprocessed_alignments