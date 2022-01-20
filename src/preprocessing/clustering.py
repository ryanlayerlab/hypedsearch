import collections
import operator
import os
from utils import ppm_to_da
from gen_spectra import get_precursor

def write_cluster(cluster, filepath, spec_num, ion):
    if len(cluster) == 0 : return None
    O = []
    O.append(len(cluster))
    O.append(cluster[0].pid)
    max_len = 0
    max_hit = None
    for hit in cluster:
        l = hit.end - hit.start + 1
        if l > max_len:
            max_len = l
            max_hit = hit
    O.append(max_hit.seq)
    O.append(max_hit.mz)
    O.append(max_hit.start)
    O.append(max_hit.end)
    for hit in cluster:
        O.append( (hit.start, hit.end, hit.seq, hit.mz) )    # b_hits
    with open(os.path.join(filepath, "spec_" + str(spec_num) + "_" + ion + '_clusters.txt'), 'a') as c:
        c.write( '\t'.join( [str(o) for o in O] ) )
        c.write('\n')

def parse_hits(Hit, file_name):
    hits = []
    for l in open(file_name):
            A = l.rstrip().split('\t')
            pid = int(A[2])
            start = int(A[4].split('-')[0])
            end = int(A[4].split('-')[1])
            seq = A[3]
            mz = float(A[1])

            hits.append( Hit(pid=pid, start=start, end=end, seq=seq, mz=mz) )
    return hits

def create_clusters(ion, filepath, spec_num):
    Hit = collections.namedtuple('Hit', 'pid start end seq mz')
    if ion == 'b':
        file_name = os.path.join(filepath, 'b_hits.txt')
        hits = parse_hits(Hit, file_name)
        sorted_hits = sorted(hits, key=operator.attrgetter('pid', 'start', 'end'))
        last_pid = None
        last_start = None
        cluster = []
        for hit in sorted_hits:
            if last_pid == hit.pid and last_start == hit.start:
                cluster.append(hit)
            else:
                write_cluster(cluster, filepath, spec_num, ion)
                cluster = [hit]
            last_pid = hit.pid
            last_start = hit.start
    else:
        file_name = os.path.join(filepath, 'y_hits.txt')
        hits = parse_hits(Hit, file_name)
        sorted_hits = sorted(hits, key=operator.attrgetter('pid', 'end', 'start'))
        last_pid = None
        last_start = None
        cluster = []
        for hit in sorted_hits:
            if last_pid == hit.pid and last_end == hit.end:
                cluster.append(hit)
            else:
                write_cluster(cluster, filepath, spec_num, ion)
                cluster = [hit]
            last_pid = hit.pid
            last_end = hit.end

def parse_indices(index_set):
    indices = []
    for index in index_set:
        string = str(index)
        A = string.rstrip().split(',')
        start = A[0]
        end = A[1]
        seq = A[2]
        mz = A[3]
        disallowed_characters = " ()\'"
        for character in disallowed_characters:
            start = start.replace(character, "")
            end = end.replace(character, "")
            seq = seq.replace(character, "")
            mz = mz.replace(character, "")
        target_tuple = (int(start), int(end), seq, float(mz))
        indices.append(target_tuple)
    return indices
def sort_clusters_by_post_prob(cluster_filepath, ion):
    cluster = collections.namedtuple('cluster', 'score pid start end seq mz indices')
    if ion == 'b':
        b_cluster_array = []
        with open(cluster_filepath, 'r') as c:
            for line in c:
                A = line.rstrip().split('\t')
                score = int(A[0])
                pid = int(A[1])
                seq = A[2]
                mz = float(A[3])
                start = int(A[4])
                end = int(A[5])
                indices = []
                [indices.append(A[x]) for x in range(6,len(A))]
                indices = parse_indices(indices)

                b_cluster_array.append(cluster(score=score, pid=pid, start=start, end=end, seq=seq, mz=mz, indices=indices) )

        b_sorted_clusters = sorted(b_cluster_array, key=operator.attrgetter('score', 'pid'), reverse = True)
        return b_sorted_clusters
    else:
        y_cluster_array = []
        with open(cluster_filepath, 'r') as c:
            for line in c:
                A = line.rstrip().split('\t')
                score = int(A[0])
                pid = int(A[1])
                seq = A[2]
                mz = float(A[3])
                start = int(A[4])
                end = int(A[5])
                indices = []
                [indices.append(A[x]) for x in range(6,len(A))]
                indices = parse_indices(indices)

                y_cluster_array.append(cluster(score=score, pid=pid, start=start, end=end, seq=seq, mz=mz, indices=indices) )

        y_sorted_clusters = sorted(y_cluster_array, key=operator.attrgetter('score', 'pid'), reverse = True)
        return y_sorted_clusters

def min_info(cluster):
    return (cluster.pid, cluster.start, cluster.end, cluster.score, cluster.seq)

def bsearch(key, Y):
        lo = -1
        hi = len(Y)
        mid = -1
        while (hi - lo > 1):
            mid = int((hi+lo) / 2)
            if Y[mid].start < key:
                lo = mid
            else:
                hi = mid
        return hi

def Ryan_merge(b_sorted_clusters, y_sorted_clusters):
    merge_seqs = []

    B = {}
    for c in b_sorted_clusters:
        if c.pid not in B:
            B[c.pid] = []
        B[c.pid].append(c)

    Y = {}
    for c in y_sorted_clusters:
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
                merge_seqs.append((b.prob * y.prob, b.end - y.start, y.end-b.start,min_info(b), min_info(y)))
                y_i += 1
    return merge_seqs

# def get_top_X(b_clusters, y_clusters, top_num):
#     filtered_b = []
#     filtered_y = []
#     b_len = top_num if len(b_clusters) >= top_num else len(b_clusters)
#     y_len = top_num if len(y_clusters) >= top_num else len(y_clusters)
#     for x in range(0,b_len):
#         filtered_b.append(b_clusters[x])
#     for x in range(0,y_len):
#         filtered_y.append(y_clusters[x])
#     return filtered_b, filtered_y

# def combine(b_cluster, y_cluster):
#     b_start, b_end, y_start, y_end = b_cluster.start, b_cluster.end, y_cluster.start, y_cluster.end
#     if b_cluster.pid == y_cluster.pid:
#         score_add = 2
#         hybrid = False
#         if  (b_end <= y_end) and (b_start <= y_start) and (b_end >= y_start): #overlap
#             overlap = True
#             score_add = 2
#             seq = b_cluster.seq
#             rem_chars = y_start - b_end
#             while (rem_chars >= 0):
#                 seq = seq + y_cluster.seq[len(y_cluster.seq)-1 - rem_chars]
#                 rem_chars = rem_chars - 1
#         else:                                                                #no overlap
#             hybrid = False
#             overlap = False
#             score_add = 2
#             seq = b_cluster.seq + '-' + y_cluster.seq
#     else:                                                                    #hybrid
#         hybrid = True
#         overlap = False
#         score_add = 0
#         seq = b_cluster.seq + '-' + y_cluster.seq
#     return seq, score_add, hybrid, overlap

# def filter_by_validity(b_cluster, y_cluster):
#     valid = True
#     for b in b_cluster.indices:
#         for y in y_cluster.indices:
#             if b[3] == y[3]:
#                 valid = False
#     return valid

# def filter_by_precursor(seq, pc, overlap, obs_prec, precursor_tolerance):
#     new_seq = seq.replace("-", "") if overlap == False else seq
#     tol = ppm_to_da(obs_prec, precursor_tolerance)
#     if get_precursor(new_seq, charge=pc) > obs_prec + tol:
#         return False
#     else:
#         return True

# def filter_by_dist(b, y, x):
#     if y.start - b.end > x:
#         return False
#     else:
#         return True

# def merge_clusters(b_sorted_clusters, y_sorted_clusters, target_precursor, precursor_tolerance):
#     # filtered_b, filtered_y = get_top_X(b_sorted_clusters, y_sorted_clusters, 50)
#     #Start with printing overlapping. Then will incorportate boundary overlaps between last of b and first of y
#     interesting_combos = []
#     for b_cluster in b_sorted_clusters:
#         for y_cluster in y_sorted_clusters:
#             if b_cluster.start <= y_cluster.end:
#                 # b_indices = parse_indices(b_cluster.indices)
#                 # y_indices = parse_indices(y_cluster.indices)
#                 if filter_by_validity(b_cluster, y_cluster):
#                     comb_seq, score_add, hybrid, overlap = combine(b_cluster, y_cluster)
#                     if filter_by_precursor(comb_seq, 2, overlap, target_precursor, precursor_tolerance):
#                         if filter_by_dist(b_cluster, y_cluster, 10):
#                             tup = (comb_seq, b_cluster.score + y_cluster.score + score_add, overlap)
#                             interesting_combos.append(tup)

#     interesting_combos.sort(key=lambda a: a[1], reverse=True)
#     return interesting_combos