def generate_truth_set(path):
    correct_sequences = []
    with open(path, 'r') as truth_set:
        for line in truth_set[1:]:
            split_line = line.split(';')
            correct_sequences.append(split_line[9])

    return correct_sequences

def check_truth_for_hybrid(b_side, y_side, truth):
    # TODO
    # Make sure this is working
    if (b_side == truth[:-(len(truth)-len(b_side))]):
        if (y_side == truth[(len(truth) - len(y_side)):]):
    # if (b_side == truth[:-(len(truth)-len(b_side))]) and (y_side == truth[(len(truth) - len(y_side)):]):
            return True
    return False

def check_truth(tup, truth):
    seq = tup[0]
    if tup[2]:
        return seq == truth
    b_side, y_side = seq.rstrip().split("-")
    return check_truth_for_hybrid(b_side, y_side, truth)

def evaluate_initial_hits(merged_seqs, SpecMill_ans, spec_num):
    for i, tup in enumerate(merged_seqs):
        if check_truth(tup, SpecMill_ans):
            print("Correct sequence found at:", i, " for spectrum ", spec_num)
            break
            