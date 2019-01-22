'''
    File name: viterbi.py
    Author: Pedro Vasco
    Date created: 12/9/2018
    Date last modified: 12/14/2013
    Python Version: 3.7
'''

import numpy as np
import sys

def vit(S, hmm, ep):
    # changes propabilities to logaritmic scale to avoid overflow
    hmm = np.log2(hmm)
    ep = np.log2(ep)

    m = hmm.shape
    m = m[0]

    s = np.zeros((len(S)), dtype=int)

    # codes sequence of genomic elements into an array of integers
    for i in range(len(S)):
        if S[i] is 'A':
            s[i] = 0
        elif S[i] is 'C':
            s[i] = 1
        elif S[i] is 'G':
            s[i] = 2
        else:
            s[i] = 3

    res = np.zeros((len(S)), dtype=int)
    score = np.zeros((m, len(S)))
    ptr = np.zeros((m, len(S)), dtype=int)
    scorei = np.zeros(m)

    # computes the first scores assuming equal initial probabiities of states
    for i in range(m):
        score[i][0] = np.log2(1 / m) + ep[s[1]][i]

    for i in range(1, len(s)):
        for j in range(m):
            for k in range(m):
                scorei[k] = score[k][i - 1] + hmm[k][j]
            ptr[j][i], = np.where(scorei == max(scorei))
            score[j][i] = ep[s[i]][j] + max(scorei)

    res[len(S) - 1], = np.where(score[:, len(S) - 1] == max(score[:, len(S) - 1]))

    # cicle to determine the most probable sequence based in the traceback technique
    for i in range(len(S) - 1, 1, -1):
        res[i - 1] = ptr[res[i]][i]

    print('Most probable sequence of states:')
    print(res)

if __name__ == '__main__':
    # state transitions matrix (markov model)
    hmm = [[0.5, 0.5], [0.4, 0.6]]
    # emition probabilities matrix with order (A-C-G-T) in rows
    ep = [[0.2, 0.3], [0.3, 0.2], [0.3, 0.2], [0.2, 0.3]]
    # observed sequence
    S = ["G", "G", "C", "A", "C", "T", "G", "A", "A"]

    vit(S, hmm, ep)

    sys.exit(0)
