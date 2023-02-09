import os
import matplotlib.pyplot as plt
import database

open_path = "/home/naco3124/jaime_hypedsearch/hypedsearch/data/output"

output_files = []
for (root, _, filenames) in os.walk(open_path):
    for fname in filenames:
        output_files.append(os.path.join(root, fname))

hybrid_list, hybrid_score_list, nat_list, nat_scores = [], [], [], []
for file in output_files:
    with open(file, 'r') as f:
        for line in f:
            A = line.split('\t')
            if A[1] == "Hybrid":
                hybrid_list.append((A[2], A[5], A[6]))
                hybrid_score_list.append(A[3])
            else:
                nat_list.append((A[2], A[5], A[6]))
                nat_scores.append(A[3])
                
prot_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/database/BMEM_searches_database.fasta'
proteins = database.build(prot_path)
                
with open("/home/naco3124/jaime_hypedsearch/hypedsearch/data/database/Adjusted_"+os.path.basename(prot_path), 'w') as d:
    for protein in proteins.proteins:
        d.write('>' + protein[0] + '\n')
        prot_seq = protein[1]
        seq_len = len(prot_seq)
        curr_write = 0
        while seq_len >= 70:
            d.write(prot_seq[curr_write:curr_write+70] + '\n')
            seq_len = seq_len-70
            curr_write = curr_write + 70
        d.write(prot_seq[curr_write:seq_len]+'\n')
    for hybrid in hybrid_list:
        parent1, parent2 = hybrid[1], hybrid[2]
        d.write(">sp|Hypedsearch Hybrid|Hybrid: "+str(parent1)+" - "+str(parent2)+"|"+hybrid[0]+"\n")
        d.write(hybrid[0].replace("-", "") + "\n\n")
        
        
hyb_index = []
for i in range(0,len(hybrid_score_list)):
    hyb_index.append(i)
plt1, ax1 = plt.subplots()
ax1.scatter(hyb_index, hybrid_score_list)
plt.xlabel("Index")
plt.ylabel("scores")
plt.title("Hybrid scores")
plt.savefig("Hybrid_scores")

nat_index = []
for i in range(0,len(nat_list)):
    nat_index.append(i)
plt2, ax2 = plt.subplots()
ax2.scatter(nat_index, nat_scores)
plt.xlabel("Index")
plt.ylabel("scores")
plt.title("Natural scores")
plt.savefig("Natural_scores")