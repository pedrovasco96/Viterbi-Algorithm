# Viterbi-Algorithm
Viterbi Algorithm for genetic sequences in MATLAB and Python

S: array of characters with each carachters as a nucleotyde
Example: S = ["G", "G", "C", "A", "C", "T", "G", "A", "A"]

hmm: state transition matrix
Example: hmm = [[0.5, 0.5], [0.4, 0.6]]
0.4 is the transition probability from state 2 to state 1

ep: matrix with emition probabilities of each nucleotyde given the state
(A-C-G-T) is the rows order. Each column represents one state
Example: ep = [[0.2, 0.3], [0.3, 0.2], [0.3, 0.2], [0.2, 0.3]]
if state is state 1, the probability that this state emits A is 0.2 and 
that emits G is 0.3
