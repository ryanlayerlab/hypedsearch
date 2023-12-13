from collections import namedtuple
import operator
from lookups.utils import ppm_to_da
import computational_pipeline.gen_spectra
from lookups.constants import WATER_MASS, PROTON_MASS
from computational_pipeline.sqlite import database_file
from lookups.constants import AMINO_ACIDS

def create_clusters_from_foos(foos):
    if len(foos) == 0: 
        return
    clusters = [len(foos), foos[0].pid]
    max_len, max_hit = 0, None
    for foo in foos:
        length = foo.end - foo.start + 1
        if length > max_len:
            max_len = length
            max_hit = foo
    clusters.extend([max_hit.mz, max_hit.start, max_hit.end, max_hit.charge, foos])
    return clusters

def parse_foos(Hit, all_hits):
    hits = []
    #Matched masses data is of form (mass, start, end, ion_int, charge, protein_num)
    for A in all_hits:
        pid = int(A[2][5])
        start = int(A[2][1])
        end = int(A[2][2])
        mz = A[1]
        charge = A[2][4]

        hits.append( Hit(pid=pid, start=start, end=end, mz=mz, charge=charge) )
    return hits

def create_clusters(ion, b_hits, y_hits):
    all_clusters = []
    Foo = namedtuple('Foo', 'pid start end mz charge')
    if ion == 'b':
        foos = parse_foos(Foo, b_hits)
        sorted_foos = sorted(foos, key=operator.attrgetter('pid', 'start', 'end'))
        last_pid = None
        last_start = None
        foos = []
        for foo in sorted_foos:
            if last_pid == foo.pid and last_start == foo.start:
                foos.append(foo)
            else:
                if foos != []:
                    clusters = create_clusters_from_foos(foos)
                    all_clusters.append(clusters)
                foos = [foo]
            last_pid = foo.pid
            last_start = foo.start
        if foos != []:
            clusters = create_clusters_from_foos(foos)
            all_clusters.append(clusters)
    else:
        foos = parse_foos(Foo, y_hits)
        sorted_foos = sorted(foos, key=operator.attrgetter('pid', 'end', 'start'))
        last_pid = None
        last_start = None
        foos = []
        for foo in sorted_foos:
            if last_pid == foo.pid and last_end == foo.end:
                foos.append(foo)
            else:
                if foos != []:
                    clusters = create_clusters_from_foos(foos)
                    all_clusters.append(clusters)
                foos = [foo]
            last_pid = foo.pid
            last_end = foo.end
        if foos != []:
            clusters = create_clusters_from_foos(foos)
            all_clusters.append(clusters)
    return all_clusters

def parse_indices(index_set):
    indices = []
    for index in index_set:
        string = str(index)
        A = string.rstrip().split(',')
        start, end, seq, mz = A[:4]
        disallowed_characters = " ()\'"
        for character in disallowed_characters:
            start = start.replace(character, "")
            end = end.replace(character, "")
            seq = seq.replace(character, "")
            mz = mz.replace(character, "")
        target_tuple = (int(start), int(end), seq, float(mz))
        indices.append(target_tuple)
    return indices

def find_sequence(pid, start_ind, end_ind, proteins):
    protein = proteins[pid]
    prot_seq = protein[1]
    target = prot_seq[start_ind:end_ind]
    return target

def append_AA(next_AA, current_mass, ion, charge):
    raw_current_mass = computational_pipeline.gen_spectra.get_raw_mass(current_mass, ion, charge)
    new_raw = raw_current_mass + AMINO_ACIDS[next_AA]
    normalized_raw = computational_pipeline.gen_spectra.calc_combined_mass(new_raw, ion)
    return normalized_raw

def test_digest_match(protein_list, pid, digest, start, end, ion):
    parent_seq = protein_list[pid][1]
    if ion == 0:
        if parent_seq[start] == digest[0]: #if we correctly cut to the left. Means sequence starts with D
            return True
        target_position = start - 1
        if not target_position < 0:
            if parent_seq[target_position] == digest[1]: #if we correctly cut to the right. Means sequence starts after an R
                return True
    else:
        if parent_seq[end-1] == digest[1]: #if we correctly cut to the right. Means sequence ends in an R
            return True
        target_position = end
        if target_position < len(parent_seq):
            if parent_seq[target_position] == digest[0]: #if we correctly cut to the left. Means sequence ends to the left of a D
                return True
    
    return False

# insane amount of nested if else in this one. needs to be changed
def find_extensions(conv_prec,current_mass,ion,charge,pid,protein_list,start,end,ppm_tolerance,seq,score):
    #goal is to get a list of all database extensions leading up to the precursor mass
    prot_seq = protein_list[pid][1]
    bad_chars = ['B', 'X', 'U', 'Z', 'O', 'J']
    extensions = []
    repeat = True
    current_seq = seq
    
    if ion == 0:
        current_position = end
        while repeat:
            if current_position < len(prot_seq):
                next_AA = prot_seq[current_position]
                if next_AA not in bad_chars:
                    current_mass = append_AA(next_AA, current_mass, ion, charge)
                    charge = 2
                    tol = ppm_to_da(current_mass, ppm_tolerance)
                    if current_mass < conv_prec + tol:
                        current_seq = current_seq + next_AA
                        tup = (current_mass, start, current_position+1, 0, 2, pid, current_seq, score)
                        extensions.append(tup)
                        current_position += 1
                    else:
                        repeat = False
                else:
                    repeat = False
            else:
                repeat = False
    
    else:
        current_position = start - 1
        while repeat:
            if current_position >= 0:
                next_AA = prot_seq[current_position]
                if next_AA not in bad_chars:
                    current_mass = append_AA(next_AA, current_mass, ion, charge)
                    charge = 2
                    tol = ppm_to_da(current_mass, ppm_tolerance)
                    if current_mass < conv_prec + tol:
                        current_seq = next_AA + current_seq
                        tup = (current_mass, current_position, end, 1, 2, pid, current_seq, score)
                        extensions.append(tup)
                        current_position = current_position - 1
                    else:
                        repeat = False
                else:
                    repeat = False
            else:
                repeat = False
    
    return extensions

def check_unique(merges):
    for m in merges:
        if merges.count(m) > 1:
            print(m, "has a count of", merges.count(m))
            return True
    return False                

def convert_components(component_arr, ion, score, seq):
    Foo = namedtuple('Foo', 'pid start end mz charge')
    converted_components = []
    if ion == 0:
        prev_component = component_arr[-1]
        for component in reversed(component_arr):
            if component.end == prev_component.end:
                converted_components.append((component.mz, component.start, component.end, ion, component.charge, component.pid, seq[:len(seq)-i], score-i))
                prev_component = component
            i = 0
            while component.end < prev_component.end:
                i +=1
                new_component = Foo(
                    prev_component.pid, prev_component.start, prev_component.end - 1, 
                    computational_pipeline.gen_spectra.max_mass(seq[:len(seq) - i], 'b', prev_component.charge), 
                    prev_component.charge
                )
                converted_components.append((new_component.mz, new_component.start, new_component.end, ion, new_component.charge, new_component.pid, seq[:len(seq)-i], score-i))
                prev_component = new_component
    else:
        prev_component = component_arr[0]
        for component in component_arr:
            if component.start == prev_component.start:
                converted_components.append((component.mz, component.start, component.end, ion, component.charge, component.pid, seq[i:], score-i))
                prev_component = component
            i = 0
            while component.start > prev_component.start + 1:
                i +=1
                new_component = Foo(
                    prev_component.pid, prev_component.start + 1, prev_component.end, 
                    computational_pipeline.gen_spectra.max_mass(seq[i:], 'y', prev_component.charge), 
                    prev_component.charge
                )
                converted_components.append((new_component.mz, new_component.start, new_component.end, ion, new_component.charge, new_component.pid, seq[i:], score-i))
                prev_component = new_component

    return converted_components

def old_score_clusters(ion, clusters, conv_prec, protein_list, prec_charge, ppm_tol):
    sorted_cluster = namedtuple('sorted_cluster', 'score pid start end mz charge components seq')
    cluster_dict = dict()
    for cluster in clusters:
        score, pid, mz, start, end, charge = cluster[:6]
        seq = find_sequence(pid, start, end, protein_list)
        components = convert_components(cluster[6], ion, score, seq)
        extensions = find_extensions(conv_prec,mz,ion,charge,pid,protein_list,start,end,ppm_tol,seq,score)
        target_cluster = sorted_cluster(score, pid, start, end, mz, charge, components + extensions, seq)
        converted_precursor = computational_pipeline.gen_spectra.convert_ion_to_precursor(mz, ion, charge, prec_charge)
        if converted_precursor not in cluster_dict:
            cluster_dict[converted_precursor] = []
        cluster_dict[converted_precursor].append(target_cluster)
        
    return cluster_dict
    
def min_info(cluster):
    return (cluster.pid, cluster.start, cluster.end, cluster.score, cluster.mz, cluster.charge, cluster.components, cluster.seq)

def bsearch(key, Y):
    low, mid, high = -1, -1, len(Y)
    while (high - low > 1):
        mid = (high + low) // 2
        if Y[mid].start < key:
            low = mid
        else:
            high = mid
    return high

def ysearch(key, B):
    low, mid, high = -1, -1, len(B)
    while (high - low > 1):
        mid = (high + low) // 2
        if B[mid].start < key:
            low = mid
        else:
            high = mid
    return high

def Ryan_merge(b_sorted_clusters, y_sorted_clusters):
    merge_seqs = list()

    B = {}
    for mz in b_sorted_clusters.keys():
        for c in b_sorted_clusters[mz]:
            if c.pid not in B:
                B[c.pid] = []
            B[c.pid].append(c)

    Y = {}
    for mz in y_sorted_clusters.keys():
        for c in y_sorted_clusters[mz]:
            if c.pid not in Y:
                Y[c.pid] = []
            Y[c.pid].append(c)

    for pid in B:
        if pid not in Y:
            continue

        sorted_B = sorted(B[pid], key=operator.attrgetter('pid', 'start', 'end'))
        sorted_Y = sorted(Y[pid], key=operator.attrgetter('pid', 'start', 'end'))
        
        for b in sorted_B:
            y_i = bsearch(b.start, sorted_Y)

            if y_i >= len(sorted_Y): break

            y = sorted_Y[y_i]

            while y_i < len(sorted_Y) and y.start - b.end < 10:
                y = sorted_Y[y_i]
                merge_seqs.append((b.score + y.score, b, y))
                y_i += 1
                    
        for y in sorted_Y:
            b_i = ysearch(y.start, sorted_B) 

            if b_i >= len(sorted_B): break

            b = sorted_B[b_i]

            while b_i < len(sorted_B) and y.start - b.end < 10:
                b = sorted_B[b_i]
                merge_seqs.append((b.score + y.score, b, y))
                b_i += 1
    return merge_seqs

def calc_from_total_overlap(side, b_mass, b_charge, y_mass, y_charge, prec_charge):
    if side: #b overlaps y
        return computational_pipeline.gen_spectra.convert_ion_to_precursor(b_mass, 0, b_charge, prec_charge)
    #y overlaps b
    return computational_pipeline.gen_spectra.convert_ion_to_precursor(y_mass, 1, y_charge, prec_charge)

def calc_from_sequences(start, y_end, pid, max_len, prec_charge):
    db = database_file(max_len, False)
    entries = db.query_sequence(pid, start, y_end)
    if not entries:
        return 0
    entry = entries.pop()
    precursor = computational_pipeline.gen_spectra.convert_ion_to_precursor(entry[0], entry[3], entry[4], prec_charge)
    return precursor
    
def total_overlap(b_pid, y_pid, b_start, y_start, b_end, y_end):
    if b_pid == y_pid:
        if b_start == y_start and b_end <= y_end:
            return True, True
        if b_end == y_end and b_start <= y_start:
            return True, True
    return False, False

def filter_by_precursor(mseqs, obs_prec, tol, precursor_charge, max_len):
    filtered_seqs = []
    for comb_seq in mseqs:
        b_pid, y_pid = comb_seq[1][0], comb_seq[2][0]
        b_start, b_end = comb_seq[1][1], comb_seq[1][2]
        y_start, y_end = comb_seq[2][1], comb_seq[2][2]
        b_charge, y_charge = comb_seq[1][5], comb_seq[2][5]
        b_mass, y_mass = comb_seq[1][4], comb_seq[2][4]
        full, side = total_overlap(b_pid, y_pid, b_start, y_start, b_end, y_end)
        if full:
            combined_precursor = calc_from_total_overlap(side, b_mass, b_charge, y_mass, y_charge, precursor_charge)
        elif b_start <= y_start and b_end <= y_end and y_start < b_end:
            combined_precursor = calc_from_sequences(b_start, y_end, b_pid, max_len, precursor_charge)
        else:
            combined_precursor = computational_pipeline.gen_spectra.calc_precursor_as_disjoint(b_mass, y_mass, b_charge, y_charge, precursor_charge)
        if combined_precursor <= obs_prec + tol:
            filtered_seqs.append(comb_seq)
    return filtered_seqs

def get_overlapping_sequence(b_seq, y_seq, b_start, b_end, y_start):
    if y_start > b_end:
        return b_seq + y_seq
    seq = b_seq[b_start:y_start]
    return seq
    
def overlap(comb_seq):
    b_seq, y_seq = comb_seq[3][4], comb_seq[4][4]
    b_pid, y_pid = comb_seq[3][0], comb_seq[4][0]
    if b_pid == y_pid:
        if b_seq == y_seq:
            return b_seq
        else:
            y_start = comb_seq[4][1]
            b_end = comb_seq[3][2]
            if (y_start - b_end > 0) and (y_start - b_end < 10):
                b_start = comb_seq[3][1]
                return get_overlapping_sequence(b_seq, y_seq, b_start, b_end, y_start)
            
    return b_seq + y_seq

def modified_find_next_mass(cluster, ion, db):
    target_index = cluster[2] + 1 if ion == 'b' else cluster[1] - 1
    target_prot = cluster[0]
    protein = db.proteins[target_prot]
    prot_seq = protein[1]
    to_add = prot_seq[target_index] if (target_index < len(prot_seq) and target_index > 0) else ''
    return to_add

# i don't think this function was finished
# what is it even doing?
def check_for_hybrid_overlap(b_seq, y_seq, ion):
    match = True
    modified_seq = ''
    if ion == 'b':
        for i, char in enumerate(b_seq):
            if char == y_seq[0]:
                k = 0
                for j in range(i, len(b_seq) + 1):
                    if b_seq[j] != y_seq[k]: # k is never changed so it always checks the first character
                        match = False
                        break
        if match:
            print('Match was true for', b_seq)
            modified_seq = b_seq[:i]
    else:
        for i, char in enumerate(y_seq):
            if char == y_seq[0]:
                k = 0
                for j in range(i, len(b_seq) + 1):
                    if b_seq[j] != y_seq[k]: # k is never changed so it always checks the first character
                        match = False
                        break
        if match:
            print('Match was true for', b_seq)
            modified_seq = b_seq[:i]
    return match, modified_seq

def grab_y_matches(indexed_clusters, target_val):
    matches = [key for key in indexed_clusters if key <= target_val] # if key is valid
    return matches
        
def get_hybrid_matches(b_sorted_clusters, y_sorted_clusters, obs_prec, precursor_tol, prec_charge):
    merged_seqs = dict()
    tol = ppm_to_da(obs_prec, precursor_tol)
    for conv_prec in b_sorted_clusters.keys():
        if conv_prec <= obs_prec + tol:
            diff = obs_prec + tol - conv_prec + (prec_charge * PROTON_MASS) + WATER_MASS
            merges = grab_y_matches(y_sorted_clusters, diff)
            merged_seqs[conv_prec] = [x for x in merges]
    return merged_seqs

# this is crazy 
def distribute_merges(merges, b_sorted_clusters, y_sorted_clusters):
    merged_clusters = []
    i, MAX_ITER = 0, 100000000
    for key in merges:
        for y_conv in merges[key]:
            for b in b_sorted_clusters[key]:
                for y in y_sorted_clusters[y_conv]:
                    if b.score + y.score > 4:
                        if i < MAX_ITER:
                            merged_clusters.append((b.score + y.score, b, y))
                            i += 1
                        else:
                            return merged_clusters

    return merged_clusters

def get_search_space(b_sorted_clusters, y_sorted_clusters, prec_charge): #This will eventually be reworked throughout the entire codebase but this is for proof of concept
    b_searches, y_searches = dict(), dict()
    for b in b_sorted_clusters.values():
        for component in b.components:
            mass = component[0]
            charge = component[4]
            prec = computational_pipeline.gen_spectra.convert_ion_to_precursor(mass, 0, charge, prec_charge)
            if prec not in b_searches:
                b_searches[prec] = []
            b_searches[prec].append(component)
                
    for y in y_sorted_clusters.values():
        for component in y.components:
            mass = component[0]
            charge = component[4]
            prec = computational_pipeline.gen_spectra.convert_ion_to_precursor(mass, 1, charge, prec_charge)
            if prec not in y_searches:
                y_searches[prec] = []
            y_searches[prec].append(component)
            
    return b_searches, y_searches