import matplotlib.pyplot as plt
import os
import sys
module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)
module_path = os.path.abspath(os.path.join('..', 'hypedsearch', 'src'))
if module_path not in sys.path:
    sys.path.append(module_path)

def check_for_second_score(rank, GoB, second_to_good, second_to_bad, score):
    if rank == 1:
        if GoB == True:
            second_to_good.append(score)
        else:
            second_to_bad.append(score)
    else:
        return None
def check_if_top_score(rank, assessment, good_top_scores, bad_top_scores):
    if (rank == 0):
        if assessment == 'True':
            good_top_scores.append(score)
            return True
        else:
            bad_top_scores.append(score)
            return False
def plot_agreement(post_array, score_array):
    fig1, ax1 = plt.subplots()
    ax1.scatter(x = post_array, y = score_array)
    plt.ylabel('Score by cluster scoring')
    plt.xlabel('Posterior probability')
    plt.title("Posterior probability vs cluster scoring")
    plt.savefig("Post_prob_vs_cluster_scoring.jpg")

def plot_rank_by_post_prob(b_bad_indexes, b_bad_probs, b_good_indexes, b_good_probs, y_bad_indexes, y_bad_probs, y_good_indexes, y_good_probs):
    fig2, ax2 = plt.subplots()
    ax2.scatter(b_bad_indexes, b_bad_probs, color = 'r', label = 'Bad hits')
    ax2.scatter(b_good_indexes, b_good_probs, color = 'g', label = 'Good hits')
    plt.title('Rank of b hits when sorted by posterior probability')
    plt.xlabel('Rank of hits')
    plt.ylabel('Posterior Probability')
    plt.legend()
    plt.savefig("Rank of b hits and Posterior Probability")
    fig3, ax3 = plt.subplots()
    ax3.scatter(y_bad_indexes, y_bad_probs, color = 'r', label = 'Bad hits')
    ax3.scatter(y_good_indexes, y_good_probs, color = 'g', label = 'Good hits')
    plt.title('Rank of y hits when sorted by posterior probability')
    plt.xlabel('Rank of hits')
    plt.ylabel('Posterior Probability')
    plt.legend()
    plt.savefig("Rank of y hits and Posterior Probability")

def plot_rank_by_score(b_bad_indexes, b_bad_scores, b_good_indexes, b_good_scores, y_bad_indexes, y_bad_scores, y_good_indexes, y_good_scores):
    fig2, ax2 = plt.subplots()
    ax2.scatter(b_bad_indexes, b_bad_scores, color = 'r', label = 'Bad hits', s=1)
    ax2.scatter(b_good_indexes, b_good_scores, color = 'g', label = 'Good hits', s=20)
    plt.title('Rank of b hits when sorted by score')
    plt.xlabel('Rank of hits')
    plt.ylabel('Score')
    plt.legend()
    plt.savefig("Rank of b hits and Score")
    fig3, ax3 = plt.subplots()
    ax3.scatter(y_bad_indexes, y_bad_scores, color = 'r', label = 'Bad hits', s = 20)
    ax3.scatter(y_good_indexes, y_good_scores, color = 'g', label = 'Good hits', s = 20)
    plt.title('Rank of y hits when sorted by score')
    plt.xlabel('Rank of hits')
    plt.ylabel('Score')
    plt.legend()
    plt.savefig("Rank of y hits and Score")

def calc_difference(score1, score2):
    return score1 - score2

def get_diff_array(good_top_scores, bad_top_scores, second_to_good, second_to_bad):
    good_diff_array = []
    bad_diff_array = []
    for i, score in enumerate(good_top_scores):
        good_diff_array.append(calc_difference(score, second_to_good[i])) 
    for i, score in enumerate(bad_top_scores):
        bad_diff_array.append(calc_difference(score, second_to_bad[i])) 

    return good_diff_array, bad_diff_array

def plot_diff_array(good_top_scores, bad_top_scores, second_to_good, second_to_bad):
    good_diff_array, bad_diff_array = get_diff_array(good_top_scores, bad_top_scores, second_to_good, second_to_bad)
    print(good_top_scores[0])
    fig1, ax1 = plt.subplots()
    ax1.scatter(good_top_scores, good_diff_array, color = 'g', label = 'Good Top Scores', s = 20)
    ax1.scatter(bad_top_scores, bad_diff_array, color = 'r', label = 'Bad Top Scores', s = 20)
    plt.title('Top_hits vs diff between top hit and second hit')
    plt.xlabel('Top_scores')
    plt.ylabel('Diff between top hit and second hit')
    plt.legend()
    plt.savefig("Top vs second best")
    # good_diff_array, bad_diff_array = get_diff_array(y_top_scores, y_second_scores)
    # fig2, ax2 = plt.subplots()
    # ax1.scatter(y_good_top_scores, good_diff_array, color = 'g', label = 'Good Top Scores', s = 20)
    # ax1.scatter(y_bad_top_scores, bad_diff_array, color = 'r', label = 'Bad Top Scores', s = 20)
    # plt.title('y_top_hits vs diff between top hit and second hit')
    # plt.xlabel('Top_scores')
    # plt.ylabel('Diff between top hit and second hit')
    # plt.legend()
    # plt.savefig("y top vs second best")

def plot_size_vs_weight(good_size_array, bad_size_array, good_weight_array, bad_weight_array, prot_size_array, prot_weight_array):
    fig1, ax1 = plt.subplots()
    ax1.scatter(good_size_array, good_weight_array, color = 'g')
    ax1.scatter(bad_size_array, bad_weight_array, color = 'r')
    plt.title('Size vs Weight')
    plt.xlabel('Size')
    plt.ylabel('Weight')
    plt.legend()
    plt.savefig("Size vs Weight")
    fig1, ax2 = plt.subplots()
    ax2.scatter(prot_size_array, prot_weight_array)
    plt.title('Prot Size vs Weight')
    plt.xlabel('Prot Size')
    plt.ylabel('Prot Weight')
    plt.legend()
    plt.savefig("Prot Size vs Weight")

print("Collecting data...")
b_good_indexes = []
b_bad_indexes = []
b_good_cluster_prob = []
b_bad_cluster_prob = []
y_good_indexes = []
y_bad_indexes = []
y_good_cluster_prob = []
y_bad_cluster_prob = []
posterior_array = []
cluster_score_array = []
b_bad_scores = []
b_good_scores = []
y_good_scores = []
y_bad_scores = []

good_top_scores = []
bad_top_scores = []
second_to_good = []
second_to_bad = []
look_for_second = False
prev_spectrum = None
prev_ion = None
GoB = False
# filepath = os.path.join('..', 'hypedsearch', 'src', 'testing_framework', 'data', 'total_data.txt') #For running locally
filepath = os.path.join('data', 'total_data.txt') #For running virtually
with open(filepath, 'r') as d:
    for i, line in enumerate(d):
        A = line.rstrip().split('\t')
        spectrum_num = int(A[0])
        score = int(A[1])
        post_prob = float(A[2])
        seq = A[3]
        rank = int(A[4])
        assessment = A[5]
        ion = A[6]
        if ion == 'b':
            posterior_array.append(post_prob)
            cluster_score_array.append(score)

            if assessment == 'True':
                b_good_scores.append(score)
                b_good_indexes.append(i)
                b_good_cluster_prob.append(post_prob)
                check_for_second_score(rank, GoB, second_to_good, second_to_bad, score)
                GoB = check_if_top_score(rank, assessment, good_top_scores, bad_top_scores)
            else:
                b_bad_scores.append(score)
                b_bad_indexes.append(i)
                b_bad_cluster_prob.append(post_prob)
                check_for_second_score(rank, GoB, second_to_good, second_to_bad, score)
                GoB = check_if_top_score(rank, assessment, good_top_scores, bad_top_scores)
        else:
            posterior_array.append(post_prob)
            cluster_score_array.append(score)

            if assessment == 'True':
                y_good_scores.append(score)
                y_good_indexes.append(i)
                y_good_cluster_prob.append(post_prob)
                check_for_second_score(rank, GoB, second_to_good, second_to_bad, score)
                GoB = check_if_top_score(rank, assessment, good_top_scores, bad_top_scores)
            else:
                y_bad_scores.append(score)
                y_bad_indexes.append(i)
                y_bad_cluster_prob.append(post_prob)
                check_for_second_score(rank, GoB, second_to_good, second_to_bad, score)
                GoB = check_if_top_score(rank, assessment, good_top_scores, bad_top_scores)            


print('Done')
print('Plotting...')

plot_agreement(posterior_array, cluster_score_array)
plot_rank_by_post_prob(b_bad_indexes, b_bad_cluster_prob, b_good_indexes, b_good_cluster_prob, y_bad_indexes, y_bad_cluster_prob, y_good_indexes, y_good_cluster_prob)
plot_rank_by_score(b_bad_indexes, b_bad_scores, b_good_indexes, b_good_scores, y_bad_indexes, y_bad_scores, y_good_indexes, y_good_scores)
plot_diff_array(good_top_scores, bad_top_scores, second_to_good, second_to_bad)
plot_size_vs_weight(good_size_array, bad_size_array, good_weight_array, bad_weight_array, prot_size_array, prot_weight_array)
print('Done')