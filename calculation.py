import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import streamlit as st
import itertools

# mostly just to ensure that the input is the same
def set_sequences(seq1, seq2):
    """
        Clean and normalize two input sequences.

        Removes whitespace, newlines, and tabs from each sequence,
        and converts them to uppercase. Intended to ensure consistent formatting
        before alignment.

        Args:
            seq1 (str): First sequence.
            seq2 (str): Second sequence.

        Returns:
            tuple: A tuple containing the cleaned versions of (seq1, seq2).
    """

    # Remove spaces, tabs, and newlines
    seq1 = seq1.replace(" ", "").replace("\n", "").replace("\t", "").upper()
    seq2 = seq2.replace(" ", "").replace("\n", "").replace("\t", "").upper()

    return seq1, seq2

def calculate_identity_percentage(msa_sequences):
    """
        Calculate average pairwise identity percentage for an MSA.

        Args:
            msa_sequences (list): List of aligned sequences.

        Returns:
            float: Average identity percentage (0–100), rounded to 2 decimal places.
    """

    num_sequences = len(msa_sequences)
    alignment_length = len(msa_sequences[0])
    total_identity = 0
    pair_count = 0

    for i in range(num_sequences):
        for j in range(i + 1, num_sequences):
            matches = 0
            valid_positions = 0
            for a, b in zip(msa_sequences[i], msa_sequences[j]):
                if a != '-' and b != '-':
                    valid_positions += 1
                    if a == b:
                        matches += 1
            if valid_positions > 0:
                identity = (matches / valid_positions) * 100
                total_identity += identity
                pair_count += 1

    average_identity = total_identity / pair_count if pair_count > 0 else 0
    return round(average_identity, 2)


# Compute MSA score
def calculate_msa_score(msa, match, mismatch, gap):
    """
        Calculate the total MSA score based on pairwise scoring.

        Args:
            msa (list): List of aligned sequences.
            match (float): Score for a match.
            mismatch (float): Penalty for mismatch.
            gap (float): Penalty for a gap.

        Returns:
            float: Total score for the MSA.
    """

    total_score = 0
    num_seqs = len(msa)
    seq_length = len(msa[0])

    for i in range(num_seqs):
        for j in range(i + 1, num_seqs):
            seq1 = msa[i]
            seq2 = msa[j]
            for k in range(seq_length):
                a = seq1[k]
                b = seq2[k]
                if a == '-' or b == '-':
                    total_score += gap
                elif a == b:
                    total_score += match
                else:
                    total_score += mismatch
    return total_score

def count_msa_statistics(msa_sequences):
    """
        Count matches, mismatches, and gaps across all pairs in the MSA.

        Args:
            msa_sequences (list): List of aligned sequences.

        Returns:
            tuple: (matches, mismatches, gaps)
    """

    alignment_length = len(msa_sequences[0])
    num_sequences = len(msa_sequences)

    matches = 0
    mismatches = 0
    gaps = 0

    for i in range(alignment_length):
        column = [seq[i] for seq in msa_sequences]
        for a, b in itertools.combinations(column, 2):
            if a == '-' or b == '-':
                gaps += 1
            elif a == b:
                matches += 1
            else:
                mismatches += 1

    return matches, mismatches, gaps
