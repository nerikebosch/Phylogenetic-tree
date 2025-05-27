import numpy as np

def matrix_building(seq1, seq2, gap_value):
    """
        Initialize the scoring matrix with gap penalties for the Needleman-Wunsch algorithm.

        Args:
            seq1 (str): First sequence.
            seq2 (str): Second sequence.
            gap_value (float): Penalty value for introducing a gap.

        Returns:
            np.ndarray: Initialized scoring matrix with gap penalties.
    """

    m = len(seq1)
    n = len(seq2)
    zero_matrix = np.zeros((m+1,n+1))

    # Fill out first column
    for i in range(0, m + 1):
        zero_matrix[i][0] = gap_value * i

    # Fill out first row
    for j in range(0, n + 1):
        zero_matrix[0][j] = gap_value * j

    return zero_matrix

def algorithm(seq1, seq2, array, gap_value, match_value, mismatch_value):
    """
        Fill in the scoring matrix using the Needleman-Wunsch algorithm.

        Args:
            seq1 (str): First sequence.
            seq2 (str): Second sequence.
            array (np.ndarray): Initialized scoring matrix.
            gap_value (float): Penalty for a gap.
            match_value (float): Reward for a match.
            mismatch_value (float): Penalty for a mismatch.

        Returns:
            np.ndarray: Completed scoring matrix.
    """

    m = len(seq1)
    n = len(seq2)

    # recursion algorithm that is found on https://rna.informatik.uni-freiburg.de/Teaching/index.jsp?toolName=Needleman-Wunsch

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            # for the row - will be penalty
            if seq1[i-1] == seq2[j-1]:
                score = array[i-1][j-1] + match_value
            else:
                score =  array[i-1][j-1] + mismatch_value

            score2 = array[i, j - 1] + gap_value
            score3 = array[i - 1, j] + gap_value
            # getting the max score
            array[i][j] = max(score, score2, score3)

    return array

def reconstruct_alignment(seq1, seq2, path):
    """
        Reconstruct two aligned sequences from a traceback path.

        Args:
            seq1 (str): First sequence.
            seq2 (str): Second sequence.
            path (list): Traceback path from the alignment matrix.

        Returns:
            tuple: (aligned_seq1, aligned_seq2), both containing possible gaps.
    """

    aligned_seq1 = ""
    aligned_seq2 = ""

    for idx in range(1, len(path)):
        i_curr, j_curr = path[idx - 1]
        i_next, j_next = path[idx]

        if i_next == i_curr + 1 and j_next == j_curr + 1:
            # diagonal: match/mismatch
            aligned_seq1 += seq1[i_curr]
            aligned_seq2 += seq2[j_curr]
        elif i_next == i_curr + 1 and j_next == j_curr:
            # up: gap in seq2
            aligned_seq1 += seq1[i_curr]
            aligned_seq2 += "-"
        elif i_next == i_curr and j_next == j_curr + 1:
            # left: gap in seq1
            aligned_seq1 += "-"
            aligned_seq2 += seq2[j_curr]

    return aligned_seq1, aligned_seq2

def traceback(score_matrix, seq1, seq2, gap, match, mismatch):
    """
        Trace back through the scoring matrix to recover the optimal alignment path.

        Args:
            score_matrix (np.ndarray): Completed scoring matrix.
            seq1 (str): First sequence.
            seq2 (str): Second sequence.
            gap (float): Gap penalty.
            match (float): Match reward.
            mismatch (float): Mismatch penalty.

        Returns:
            list: List of (i, j) tuple indices representing the optimal alignment path.
    """

    i, j = len(seq1), len(seq2)
    path = [(i, j)]

    while i > 0 or j > 0:
        current_score = score_matrix[i, j]

        if i > 0 and j > 0:
            diag = score_matrix[i-1, j-1]
            score_diag = match if seq1[i-1] == seq2[j-1] else mismatch
            if current_score == diag + score_diag:
                i -= 1
                j -= 1
                path.append((i, j))
                continue

        if i > 0 and current_score == score_matrix[i-1, j] + gap:
            i -= 1
            path.append((i, j))
            continue

        if j > 0 and current_score == score_matrix[i, j-1] + gap:
            j -= 1
            path.append((i, j))
            continue

    path.reverse()
    return path