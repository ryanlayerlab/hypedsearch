import computational_pipeline.gen_spectra
from lookups.utils import ppm_to_da
import lookups.objects

def overlap_scoring(sequence, input_masses, ppm_tolerance):
    total_score = 0
    spectrum = computational_pipeline.gen_spectra.generate_spectrum(sequence)
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

def get_target_data(aligned_spectrum_params):
    spectrum = aligned_spectrum_params.spectrum
    base_alignment_params = aligned_spectrum_params.base_alignment_params
    precursor_charge = spectrum.precursor_charge
    precursor_mass = spectrum.precursor_mass
    target_seq = base_alignment_params.target_seq
    sqllite_database = base_alignment_params.sqllite_database
    ppm_tolerance = base_alignment_params.ppm_tolerance
    precursor_tolerance = base_alignment_params.precursor_tolerance

    protein_list = sqllite_database.query_proteins()
    input_masses = spectrum.mz_values
    theoretical_prec = computational_pipeline.gen_spectra.get_precursor(target_seq.replace("-", ""), precursor_charge)
    if abs(theoretical_prec - precursor_mass) <= precursor_tolerance:
        print(target_seq,"is a valid hybrid to find")
    else:
        print(target_seq,"is NOT a valid hybrid to find")
    
    target_left_pids, target_right_pids = [], []
    target_left_indices, target_right_indices = [], []
    if "-" in target_seq:
        b_part, y_part = target_seq.split("-")
    else:
        b_part, y_part = target_seq, target_seq
    for i, protein in enumerate(protein_list):
        protein_seq = protein[1]
        if b_part in protein_seq:
            start = protein_seq.find(b_part)
            end = start + len(b_part)
            print(protein_seq[start:end])
            target_left_pids.append(i)
            target_left_indices.append((start, end))
        if y_part in protein_seq:
            start = protein_seq.find(y_part)
            end = start + len(y_part)
            print(protein_seq[start:end])
            target_right_pids.append(i)
            target_right_indices.append((start, end))
            
    target_score = overlap_scoring(target_seq.replace("-", ""), input_masses, ppm_tolerance)

    target_data = lookups.objects.TargetData(target_seq=target_seq,
        target_left_pids=target_left_pids,
        target_right_pids=target_right_pids,
        target_left_indices=target_left_indices,
        target_right_indices=target_right_indices,
        target_score=target_score)

    return target_data

def check_in_matched_masses(alignment_data,target_data):
    converted_precursor_b = alignment_data.converted_precursor_b
    converted_precursor_y = alignment_data.converted_precursor_y
    matched_masses_b = alignment_data.matched_masses_b
    matched_masses_y = alignment_data.matched_masses_y
    left_pids = target_data.target_left_pids
    left_indices = target_data.target_left_indices
    right_pids = target_data.target_right_pids
    right_indices = target_data.target_right_indices
    good_b_entries, good_y_entries = [], []
    for i, pid in enumerate(left_pids):
        for key in matched_masses_b.keys():
            for entry in matched_masses_b[key]:
                entry_pid = entry[5]
                entry_left_start = entry[1]
                if pid == entry_pid:
                    if left_indices[i][0] == entry_left_start:
                        good_b_entries.append(entry)
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
    target_alignment_data = lookups.objects.AlignmentData(converted_precursor_b=converted_precursor_b,converted_precursor_y=converted_precursor_y,matched_masses_b=good_b_entries,matched_masses_y=good_y_entries)
    return target_alignment_data
    
def check_in_sorted_clusters(target_alignment_data,sorted_clusters):
    b_sorted_clusters = sorted_clusters.b_sorted_clusters
    y_sorted_clusters = sorted_clusters.y_sorted_clusters
    good_b_hits = target_alignment_data.matched_masses_b
    good_y_hits = target_alignment_data.matched_masses_y
    good_maxb_hit, good_maxy_hit = [], []
    good_maxb_dict, good_maxy_dict = dict(), dict()
    good_b_clusters, good_y_clusters = [], []
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
                        good_b_clusters.append((mz, cluster))

    for hit in good_maxy_hit:
        for mz in y_sorted_clusters.keys():
            for cluster in y_sorted_clusters[mz]:
                if cluster.pid == hit[5]:
                    if cluster.end == hit[2]:
                        good_y_clusters.append((mz, cluster))

    if len(good_b_clusters) > 0 and len(good_y_clusters) > 0:
        print("Good hits in b_clusters and y_clusters")
    if len(good_b_clusters) > 0 and len(good_y_clusters) == 0:
        print("Good hits in b_clusters but not y_clusters")
    if len(good_b_clusters) == 0 and len(good_y_clusters) > 0:
        print("No good hits in b_clusters but in y_clusters")
    if len(good_b_clusters) == 0 and len(good_y_clusters) == 0:
        print("No good hits in b_clusters nor in y_clusters")
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
    for score, left_native, right_native in native_alignments:
        left_pid = left_native.pid
        left_start, left_end = left_native.start, left_native.end
        right_pid = right_native.pid
        right_start, right_end = right_native.start, right_native.end
        if left_pid in target_left_pids and (left_start, left_end) in target_left_indices and right_pid in target_right_pids and (right_start, right_end) in target_right_indices:
            good_natives.append((score, left_native, right_native))
    
    for score, left_hybrid, right_hybrid in hybrid_alignments:
        left_pid = left_hybrid.pid
        left_start, left_end = left_hybrid.start, left_hybrid.end
        right_pid = right_hybrid.pid
        right_start, right_end = right_hybrid.start, right_hybrid.end
        if left_pid in target_left_pids and (left_start, left_end) in target_left_indices and right_pid in target_right_pids and (right_start, right_end) in target_right_indices:
            good_hybrids.append((score, left_hybrid, right_hybrid))
            
    if len(good_natives) > 0:
        print("Good native alignments were found")
        
    if len(good_hybrids) > 0:
        print("Good hybrid alignments were found")
        
    else:
        print("No good native or hybrid alignments were found")
    
    return good_natives, good_hybrids

def check_score(rescored_natives, rescored_hybrids, good_natives, good_hybrids, target_score):
    good_scored = []
    for native in good_natives:
        for renative in rescored_natives:
            if native[1] == renative[2] and native[2] == renative[3]:
                native_score = native[0]
                if native_score != target_score:
                    print("Scores were not equal for this native. target score was:", target_score, "while Hypedsearch scored:", native_score)
                good_scored.append((native_score, renative[2], renative[3]))
            
    for hybrid in good_hybrids:
        for rehybrid in rescored_hybrids:
            if hybrid[1] == rehybrid[2] and hybrid[2] == rehybrid[3]:
                hybrid_score = rehybrid[0]
                if hybrid_score != target_score:
                    print("Scores were not equal for this hybrid. target score was:", target_score, "while Hypedsearch scored:", hybrid_score)
                good_scored.append((hybrid_score, rehybrid[2], rehybrid[3]))
    
    return good_scored

def check_in_rescored(rescored, good_scored):
    good_rescored = []
    for score, left_alignment, right_alignment in good_scored:
        for new_score, _, rescored_left, rescored_right, _ in rescored:
            if score == new_score and left_alignment == rescored_left and right_alignment == rescored_right:
                good_rescored.append((new_score, rescored_left, rescored_right))
                
    if len(good_rescored) > 0:
        print("Good rescored alignments were found")
    
    else:
        print("Lost in final filtering")
        
        
def check_in_searches(target_data,target_alignment_data,ppm_tolerance,precursor_charge,target_seq):
    b_searches =  target_alignment_data.matched_masses_b
    y_searches = target_alignment_data.matched_masses_y
    target_left_pids = target_data.target_left_pids
    target_right_pids = target_data.target_right_pids
    target_left_indices = target_data.target_left_indices
    target_right_indices = target_data.target_right_indices
    good_b_searches, good_y_searches = [],[]
    if "-" in target_seq:
        left_part, right_part = target_seq.split("-")
    else:
        left_part, right_part = target_seq, target_seq
    left_prec = computational_pipeline.gen_spectra.get_precursor(left_part, precursor_charge)
    right_prec = computational_pipeline.gen_spectra.get_precursor(right_part, precursor_charge)
    b_tol = ppm_to_da(left_prec, ppm_tolerance)
    y_tol = ppm_to_da(right_prec, ppm_tolerance)
    
    for prec in b_searches.keys():
        if abs(prec - left_prec) < b_tol:
            for cluster in b_searches[prec]:
                cluster_pid = cluster[5]
                cluster_start = cluster[1]
                cluster_end = cluster[2]
                if cluster_pid in target_left_pids:
                    if (cluster_start, cluster_end) in target_left_indices:
                        good_b_searches.append((prec, cluster))
                        
    for prec in y_searches.keys():
        if abs(prec - right_prec) < y_tol:
            for cluster in y_searches[prec]:
                cluster_pid = cluster[5]
                cluster_start = cluster[1]
                cluster_end = cluster[2]
                if cluster_pid in target_right_pids:
                    if (cluster_start, cluster_end) in target_right_indices:
                        good_y_searches.append((prec, cluster))
    
    if len(good_b_searches) > 0 and len(good_y_searches) > 0:
        print("Good hits in b_searches and y_searches")
    if len(good_b_searches) > 0 and len(good_y_searches) == 0:
        print("Good hits in b_searches but not y_searches")
    if len(good_b_searches) == 0 and len(good_y_searches) > 0:
        print("No good hits in b_searches but in y_searches")
    if len(good_b_searches) == 0 and len(good_y_searches) == 0:
        print("No good hits in b_searches nor in y_searches")
    return good_b_searches, good_y_searches

def check_in_merges(good_searches):
    hybrid_merges = check_merges_params.hybrid_merges
    native_merges = check_merges_params.native_merges
    good_b_searches = good_searches.good_b_searches
    good_y_searches = good_searches.good_y_searches
    good_merges = []
    for key in hybrid_merges:
        for merge in hybrid_merges[key]:
            left_part, right_part = merge[0], merge[1]
            for good_b in good_b_searches:
                if left_part == good_b[1]:
                    for good_y in good_y_searches:
                        if right_part == good_y[1]:
                            good_merges.append((left_part, right_part, 1))
    
    for key in native_merges:
        for merge in native_merges[key]:
            left_part, right_part = merge[0], merge[1]
            for good_b in good_b_searches:
                if left_part == good_b[1]:
                    good_merges.append((left_part, right_part, 0))
                        
    if len(good_merges) > 0:
        print("Good merges were found")
    else:
        print("No good merges were found")
    
    return good_merges

def check_in_rescored_merges(rescored_merges,good_merges):
    good_rescored = []
    sorted_rescored = sorted(rescored_merges, key=lambda x: (x[0], x[1]), reverse = True)
    for merge in good_merges:
        for i, remerged in enumerate(sorted_rescored):
            if merge[0][6] != merge[1][6]:
                if merge[0][6] + merge[1][6] == remerged[2]:
                    good_rescored.append((i, merge))
                    break
            else: 
                if merge[0][6] == remerged[2]:
                    good_rescored.append((i, merge))
                    break
    good = False
    for index, merge in good_rescored:
        if index < 10:
            good = True
            break
        
    if len(good_rescored) == 0:
        print("There were no good rescored")
    elif not good and len(good_rescored) > 0:
        print("There were good_rescored but they all score outside of the top 10. Best scores in top", good_rescored[0][0])
    else:
        print("There were good rescored and it scores in the top 10")
    
    return good_rescored

def check_in_unique(unique_merges, good_merges):
    good_unique = []
    for merge in good_merges:
        if merge[2]: #hybrid
            good_seq = merge[0][6] + merge[1][6]
        else:
            good_seq = merge[0][6]
        for (seq, score) in unique_merges:
            if seq == good_seq:
                good_unique.append(((seq, score), unique_merges[(seq, score)]))
    
    if len(good_unique) > 0:
        print("After uniqueness, there are good merges")
    else:
        print("There were no good merges after uniqueness")
    
    return good_unique