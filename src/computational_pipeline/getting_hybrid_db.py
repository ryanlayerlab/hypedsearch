import os
import matplotlib.pyplot as plt
import computational_pipeline.database_generator as database_generator
from postprocessing.postprocessing_utils import make_db_mapping_by_key

open_path = "/home/naco3124/jaime_hypedsearch/hypedsearch/data/truth_table/Results"

output_files = []
for (root, _, filenames) in os.walk(open_path):
    for fname in filenames:
        output_files.append(os.path.join(root, fname))

hybrid_list, hybrid_score_list, nat_list, nat_scores = set(), [], [], []
for file in output_files:
    with open(file, 'r') as f:
        for line in f:
            A = line[1:len(line)].split(', ') #This code needs to be reworked if only grabbing hypedsearch hybrids
            verdict = A[1][:len(A[1])].split(": ")
            if verdict[1] == "'Hybrid'":
                sequence_entry = A[2][:len(A[2])].split(": ")
                sequence = str(sequence_entry[1][1:len(sequence_entry[1])-1])
                extended_sequence_entry = A[9][:len(A[9])].split(": ")
                extended_sequence = str(extended_sequence_entry[1][1:len(extended_sequence_entry[1])-3])
                left_parent_entry = A[5][:len(A[5])].split(": ")
                left_parent = str(left_parent_entry[1][1:len(left_parent_entry[1])-1])
                right_parent_entry = A[6][:len(A[6])].split(": ")
                right_parent = str(right_parent_entry[1][1:len(right_parent_entry[1])-1])
                hybrid_list.add((extended_sequence, left_parent, right_parent, sequence))
                # hybrid_score_list.append(A[3])
            else:
                nat_list.append((A[2], A[5], A[6], A[2]))
                # nat_scores.append(A[3])
                
def further_extend(hybrid_list, proteins):
    new_hybrids = set()
    db_mapping = make_db_mapping_by_key(proteins)
    for hybrid in hybrid_list:
        left_parent, right_parent = hybrid[1], hybrid[2]
        left_part, right_part = hybrid[0].split("-")
        full_left_parent, full_right_parent = db_mapping[left_parent],db_mapping[right_parent]
        left_parent_seq, right_parent_seq = full_left_parent[1], full_right_parent[1]
        left_start = left_parent_seq.find(left_part)
        left_extensions = left_parent_seq[max(left_start-20, 0):left_start]
        left_extended = left_extensions + left_part
        right_end = right_parent_seq.find(right_part) + len(right_part)
        right_extensions = right_parent_seq[right_end:min(right_end+25,len(right_parent_seq))]
        right_extended = right_part + right_extensions
        
        new_hybrids.add((left_extended+"-"+right_extended, left_parent, right_parent, hybrid[3]))
    
    return new_hybrids
                
prot_path = '/home/naco3124/jaime_hypedsearch/hypedsearch/data/database/UniProt_mouse.fasta'
proteins = database_generator.build_database(prot_path)

hybrid_list = further_extend(hybrid_list, proteins)
                
with open("/home/naco3124/jaime_hypedsearch/hypedsearch/data/database/Hybrid_"+os.path.basename(prot_path), 'w') as d:
    for hybrid in hybrid_list:
        parent1, parent2 = hybrid[1], hybrid[2]
        d.write(">sp|Hypedsearch Hybrid|Hybrid: " +parent1+" - "+parent2+"|" + hybrid[3]+"\n")
        replaced_hybrid = hybrid[0].replace("-", "")
        if len(replaced_hybrid) > 70:
            target_string = replaced_hybrid[:70] + "\n" + replaced_hybrid[70:]
        else:
            target_string = hybrid[0]
        d.write(target_string.replace("-", "") + "\n\n")
 
# hyb_index = []
# for i in range(0,len(hybrid_score_list)):
#     hyb_index.append(i)
# plt1, ax1 = plt.subplots()
# ax1.scatter(hyb_index, hybrid_score_list)
# plt.xlabel("Index")
# plt.ylabel("scores")
# plt.title("Hybrid scores")
# plt.savefig("Hybrid_scores")

# nat_index = []
# for i in range(0,len(nat_list)):
#     nat_index.append(i)
# plt2, ax2 = plt.subplots()
# ax2.scatter(nat_index, nat_scores)
# plt.xlabel("Index")
# plt.ylabel("scores")
# plt.title("Natural scores")
# plt.savefig("Natural_scores")