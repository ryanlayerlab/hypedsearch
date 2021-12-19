import sys
import os
import collections
import operator

module_path = os.path.abspath(os.path.join('..', '..'))
if module_path not in sys.path:
    sys.path.append(module_path)
module_path = os.path.abspath(os.path.join('..', 'hypedsearch', 'src'))
if module_path not in sys.path:
    sys.path.append(module_path)

Cluster = collections.namedtuple('Cluster', 'score pid seq mass start end ion hits')

import collections
import operator


Cluster = collections.namedtuple('Cluster', 'score pid seq mass start end ion hits')

def get_seq(f, pid, start, end):
    return f[pid][start-1:end]

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

fasta_file = os.path.abspath(os.path.join('../../../data/database/sample_database.fasta'))



solution_array = []

B = {}
for l in b_sorted_clusters:
    if l.pid not in B:
        B[l.pid] = []

    B[l.pid].append(l)

Y = {}
for l in y_sorted_clusters:
    if l.pid not in Y:
        Y[l.pid] = []
    Y[l.pid].append(l)

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
            seq = get_seq(fa, b.pid, b.start, y.end)
            tup = (b.score + y.score, b.end - y.start, y.end - b.start, seq, min_info(b), min_info(y))
            solution_array.append(tup)
            y_i += 1
            
solution_array = sorted(solution_array, key = lambda x: x[0], reverse=True)
with open('Solutions.txt', 'w') as s:
    [s.write(x + '\n') for x in solution_array]