import sys
import os
import collections
import operator

module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
module_path = os.path.abspath(os.path.join('..', 'hypedsearch', 'src'))
if module_path not in sys.path:
    sys.path.append(module_path)

Cluster = collections.namedtuple('Cluster', 'score pid seq mass start end ion hits')

def load_fasta(fasta_file):
    f = []
    for l in open(fasta_file):
        if l[0] == '>':
            continue
        f.append(l.rstrip())
    return f

def get_seq(f, pid, start, end):
    return f[pid][start-1:end]

def min_info(cluster):
    return (cluster.pid, cluster.start, cluster.end, cluster.score, cluster.seq, cluster.ion)

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

b_file = os.path.abspath(os.path.join("src/testing_framework/b_sorted_clusters.txt"))
y_file = os.path.abspath(os.path.join("src/testing_framework/y_sorted_clusters.txt"))
fasta_file = os.path.abspath(os.path.join("data/database/sample_database.fasta"))

fa = load_fasta(fasta_file)

B = {}
for l in open(b_file):
    A = l.rstrip().split()
    c = Cluster(score=int(A[2]),
                pid=int(A[3]),
                seq=A[4],
                mass=float(A[5]),
                start=int(A[6]),
                end=int(A[7]),
                ion='b',
                hits=A[8:])
    if c.pid not in B:
        B[c.pid] = []

    B[c.pid].append(c)

Y = {}
for l in open(y_file):
    A = l.rstrip().split()
    c = Cluster(score=int(A[2]),
                pid=int(A[3]),
                seq=A[4],
                mass=float(A[5]),
                start=int(A[6]),
                end=int(A[7]),
                ion='y',
                hits=A[8:])
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
            seq = get_seq(fa, b.pid, b.start, y.end)
            print(b.score + y.score, b.end - y.start, y.end - b.start, seq, min_info(b), min_info(y))
            y_i += 1