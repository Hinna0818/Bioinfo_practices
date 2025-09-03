## two sequence alignment globally
import numpy as np

# prepare score_matrix for global 2-seq alignment
def matrix_global(seq1, seq2, rules = [4, 3, 4]):
    if seq1 is None or seq2 is None:
        return 0
    
    len1 = len(seq1) + 1
    len2 = len(seq2) + 1
    match = rules[0]
    mismatch = rules[1]
    space = rules[2]
    matrix = np.zeros((len1, len2), dtype=int)

    ## set penalty for the first row and column
    for i in range(len1):
        if i == 0:
            matrix[i, 0] = 0
        else:
            matrix[i, 0] = -i * space
    
    for j in range(len2):
        if j == 0:
            matrix[0, j] = 0
        else:
            matrix[0, j] = -j * space
    
    ## initialize score_matrix
    for i in range(1, len1):
        for j in range(1, len2):
            diag = left = right = 0
            diag = matrix[i-1, j-1] + (match if seq1[i-1] == seq2[j-1] else - mismatch)
            left = matrix[i, j-1] - space
            right = matrix[i-1, j] - space
            matrix[i, j] = max(diag, left, right)
    
    return matrix

def traceback_global(seq1, seq2, score_matrix, rules = [4, 3, 4]):
    match, mismatch, space = rules
    i, j = len(seq1), len(seq2)
    aligned1, aligned2 = [], []

    while i > 0 or j > 0:
        current = score_matrix[i, j]
        if i > 0 and j > 0:
            diag = score_matrix[i-1, j-1]
            score = match if seq1[i-1] == seq2[j-1] else -mismatch
            if current == diag + score:
                aligned1.append(seq1[i-1])
                aligned2.append(seq2[j-1])
                i -= 1
                j -= 1
                continue
        
        if i > 0:
            up = score_matrix[i-1, j]
            if current == up - space:
                aligned1.append(seq1[i-1])
                aligned2.append('-')
                i -= 1
                continue
        
        if j > 0:
            left = score_matrix[i, j-1]
            if current == left - space:
                aligned1.append('-')
                aligned2.append(seq2[j-1])
                j -= 1
    
    return aligned1[::-1], aligned2[::-1]


if __name__ == "__main__":
    a1 = ['a', 't', 't', 'c', 'c', 'a', 'a', 'g']
    a2 = ['t', 't', 'c', 'g', 'a', 'g', 't']

    res = matrix_global(a1, a2)
    align1, align2 = traceback_global(a1, a2, res)
    print(res)
    print('Alignment 1:', ''.join(align1))
    print('Alignment 2:', ''.join(align2))