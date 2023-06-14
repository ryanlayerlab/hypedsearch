import gen_spectra
from utils import ppm_to_da

def overlap_scoring(sequence, input_masses, ppm_tolerance):
    total_score = 0
    spectrum = gen_spectra.gen_spectrum(sequence)
    masses = sorted(spectrum['spectrum'])
    input_masses = sorted(input_masses)
    o_ctr, t_ctr = 0, 0
    observed = input_masses[o_ctr]
    theoretical = masses[t_ctr]
    while (o_ctr < len(input_masses) and t_ctr < len(masses)):
        tol = ppm_to_da(observed, ppm_tolerance)
        if theoretical < observed - tol:
            t_ctr = t_ctr + 1
            if t_ctr < len(masses):
                theoretical = masses[t_ctr]
        elif observed + tol < theoretical:
            o_ctr = o_ctr + 1
            if o_ctr < len(input_masses):
                observed = input_masses[o_ctr]
        elif abs(observed-theoretical) <= tol:
            total_score = total_score + 1
            o_ctr = o_ctr + 1
            t_ctr = t_ctr + 1
            if o_ctr < len(input_masses) and t_ctr < len(masses):
                observed = input_masses[o_ctr]
                theoretical = masses[t_ctr]
                
    return(total_score)

def get_target_data(target_seq, proteins, input_masses, ppm_tolerance):
    target_left_pids, target_right_pids = [], []
    target_left_indices, target_right_indices = [], []
    protein_list = proteins.proteins
    b_part, y_part = target_seq.split("-")
    for i, protein in enumerate(protein_list):
        protein_seq = protein[1]
        if b_part in protein_seq:
            target_left_pids.append(i)
            start = protein_seq.find(b_part)
            end = start + len(b_part)
            print(protein_seq[start:end])
            target_left_indices.append((start, end))
        if y_part in protein_seq:
            target_right_pids.append(i)
            start = protein_seq.find(y_part)
            end = start + len(y_part)
            print(protein_seq[start:end])
            target_right_indices.append((start, end))
            
    target_score = overlap_scoring(target_seq.replace("-", ""), input_masses, ppm_tolerance)
    
    return target_seq, target_left_pids, target_right_pids, target_left_indices, target_right_indices, target_score

def check_in_matched_masses(matched_masses_b, matched_masses_y, left_pids, left_indices, right_pids, right_indices):
    #left part
    good_b_entries, good_y_entries = [], []
    for i, pid in enumerate(left_pids):
        for key in matched_masses_b.keys():
            for entry in matched_masses_b[key]:
                entry_pid = entry[5]
                entry_left_start = entry[1]
                if pid == entry_pid:
                    if left_indices[i][0] == entry_left_start:
                        good_b_entries.append(entry)
    #right part
    for i, pid in enumerate(right_pids):
        for key in matched_masses_y.keys():
            for entry in matched_masses_y[key]:
                entry_pid = entry[5]
                entry_right_end = entry[2]
                if pid == entry_pid and right_indices[i][1] == entry_right_end:
                    good_y_entries.append(entry)
            
    if len(good_b_entries) > 0 and len(good_y_entries) > 0:
        print("Good hits in matched_masses_b and matched_masses_y")
    if len(good_b_entries) > 0 and len(good_y_entries) == 0:
        print("Good hits in matched_masses_b but not matched_masses_y")
    if len(good_b_entries) == 0 and len(good_y_entries) > 0:
        print("No good hits in matched_masses_b but in matched_masses_y")
    if len(good_b_entries) == 0 and len(good_y_entries) == 0:
        print("No good hits in matched_masses_b nor in matched_masses_y")
    return good_b_entries, good_y_entries
    
def check_in_sorted_clusters(b_sorted_clusters, y_sorted_clusters, good_b_hits, good_y_hits, target_seq):
    good_maxb_hit, good_maxy_hit = [], []
    good_maxb_dict, good_maxy_dict = dict(), dict()
    good_b_clusters, good_y_clusters = [], []
    b_part, y_part = target_seq.split("-")
    for hit in good_b_hits:
        pid = hit[5]
        if pid not in good_maxb_dict.keys():
            good_maxb_dict[pid] = []
        good_maxb_dict[pid].append(hit)
            
    for hit in good_y_hits:
        pid = hit[5]
        if pid not in good_maxy_dict.keys():
            good_maxy_dict[pid] = []
        good_maxy_dict[pid].append(hit)
    
    for key in good_maxb_dict.keys():
        best_hit = max(good_maxb_dict[key], key=lambda x: x[2])
        good_maxb_hit.append(best_hit)
        
    for key in good_maxy_dict.keys():
        best_hit = min(good_maxy_dict[key], key=lambda x: x[1])
        good_maxy_hit.append(best_hit)
    
    for hit in good_maxb_hit:
        for mz in b_sorted_clusters.keys():
            for cluster in b_sorted_clusters[mz]:
                if cluster.pid == hit[5]:
                    if cluster.start == hit[1]:
                        if cluster.end > hit[1]+len(b_part):
                            print("check this")
                        else:
                            good_b_clusters.append((mz, cluster))

    for hit in good_maxy_hit:
        for mz in y_sorted_clusters.keys():
            for cluster in y_sorted_clusters[mz]:
                if cluster.pid == hit[5]:
                    if cluster.end == hit[2]:
                        if cluster.start < hit[2]-len(y_part):
                            print("check this")
                        else:
                            good_y_clusters.append((mz, cluster))
    # want to return good clusters and print whether these good clusters are in our sorted list

    if len(good_b_clusters) > 0 and len(good_y_clusters) > 0:
        print("Good hits in matched_masses_b and matched_masses_y")
    if len(good_b_clusters) > 0 and len(good_y_clusters) == 0:
        print("Good hits in matched_masses_b but not matched_masses_y")
    if len(good_b_clusters) == 0 and len(good_y_clusters) > 0:
        print("No good hits in matched_masses_b but in matched_masses_y")
    if len(good_b_clusters) == 0 and len(good_y_clusters) == 0:
        print("No good hits in matched_masses_b nor in matched_masses_y")
    return good_b_clusters, good_y_clusters

def check_in_natives(natives, good_b, good_y):
    good_natives = []
    for native in natives:
        for b in good_b:
            for y in good_y:
                if b == native[1] and y == native[2]:
                    good_natives.append(native)
    
    if len(good_natives) > 0:
        print("good natives were found")
        
    else:
        print("no good natives were found")
    
    return good_natives

def check_in_hybrids(hybrid_dict, good_b, good_y, all_b, all_y):
            
    for x in all_b:
        for b in all_b[x]:
            if all_b[x].count(b) > 1:
                print("b_sorted_clusters are non unique")
                break
        
    for x in all_y:
        for y in all_y[x]:
            if all_y[x].count(y) > 1:
                print("y_sorted_clusters are non unique")
                break
    
    good_hybrids = []
    for b in good_b:
        conv_b = b[0]
        for y in good_y:
            conv_y = y[0]
            if conv_b in hybrid_dict.keys():
                if conv_y in hybrid_dict[conv_b]:
                    good_hybrids.append((b[1].score + y[1].score, b[1], y[1]))
            
           
    
    if len(good_hybrids) > 0:
        print("good hybrids were found")
        
    else:
        print("no good hybrids were found")
    
    return good_hybrids

def check_in_combined_hybrids(good_hybrids, b_merged_hybrids):
    good_combined = []
    for hybrid in good_hybrids:
        if hybrid in b_merged_hybrids:
            good_combined.append(hybrid)
            
    if len(good_combined) > 0:
        print("good combined hybrids were found")
        
    else:
        print("no good combined hybrids were found")
    
    return good_combined

def check_in_alignments(target_left_pids, target_left_indices, target_right_pids, target_right_indices, native_alignments, hybrid_alignments):
    good_natives, good_hybrids = [],[]
    for left_native, right_native in native_alignments:
        left_pid = left_native[5]
        left_start, left_end = left_native[1], left_native[2]
        right_pid = right_native[5]
        right_start, right_end = right_native[1], right_native[2]
        if left_pid in target_left_pids and (left_start, left_end) in target_left_indices and right_pid in target_right_pids and (right_start, right_end) in target_right_indices:
            good_natives.append((left_native, right_native))
    
    for left_hybrid, right_hybrid in hybrid_alignments:
        left_pid = left_hybrid[5]
        left_start, left_end = left_hybrid[1], left_hybrid[2]
        right_pid = right_hybrid[5]
        right_start, right_end = right_hybrid[1], right_hybrid[2]
        if left_pid in target_left_pids and (left_start, left_end) in target_left_indices and right_pid in target_right_pids and (right_start, right_end) in target_right_indices:
            good_hybrids.append((left_hybrid, right_hybrid))
            
    if len(good_natives) > 0:
        print("Good native alignments were found")
        
    if len(good_hybrids) > 0:
        print("Good hybrid alignments were found")
        
    else:
        print("No good native or hybrid alignments were found")
    
    return good_natives, good_hybrids

def check_score(good_natives, good_hybrids, target_score):
    for native in good_natives:
        native_score = native[0]
        if native_score != target_score:
            print("Scores were not equal. target score was:", target_score, "while Hypedsearch scored:", native_score)
            
    for hybrid in good_hybrids:
        hybrid_score = hybrid[0]
        if hybrid_score != target_score:
            print("Scores were not equal. target score was:", target_score, "while Hypedsearch scored:", hybrid_score)