import datetime

from calculation import *
from file_utils import *
from needleman_algorithm import *
from msa_algorithm import *
from upgma_algorithm import *
from phytree import *
from graph import *
import streamlit as st
import numpy as np
import itertools


def app_creation():
    """
        Launches the Streamlit web application for MSA and phylogenetic tree generation.

        This function handles:
        - Sequence input via text fields or file upload.
        - Distance matrix input via file upload.
        - Running MSA (Star Algorithm) with Needleman-Wunsch alignment.
        - Calculating alignment statistics and score.
        - Displaying and downloading distance matrix and UPGMA tree in Newick format.
        - Rendering phylogenetic tree graphics and downloadable report.
    """

    st.set_page_config(layout="wide")
    st.title("MSA Star Algorithm")

    # Init session state
    if "sequences" not in st.session_state:
        st.session_state.sequences = [""]
        st.session_state.names_of_sequences = [""]

    if "distance_matrix" not in st.session_state:
        st.session_state.distance_matrix = None

    def add_sequence():
        st.session_state.sequences.append("")
        st.session_state.names_of_sequences.append("")

    def remove_sequence(index):
        if len(st.session_state.sequences) > 1:
            st.session_state.sequences.pop(index)
            st.session_state.names_of_sequences.pop(index)

    new_sequences = []
    new_name_list = []

    st.markdown("### Input Sequences")
    col1, col2 = st.columns(2)
    with col1:
        for i, seq in enumerate(st.session_state.sequences):
            cols = st.columns([3, 3, 1])
            name = cols[0].text_input(f"Name {i + 1}", value=st.session_state.names_of_sequences[i], key=f"name_{i}")
            new_name_list.append(name)
            new_value = cols[1].text_input(f"Sequence {i + 1}", value=seq, key=f"seq_{i}")
            new_sequences.append(new_value)
            if len(st.session_state.sequences) > 1 and cols[2].button("❌", key=f"remove_{i}"):
                remove_sequence(i)
                st.rerun()

        st.session_state.sequences = new_sequences
        st.session_state.names_of_sequences = new_name_list
        if st.button("➕ Add a new sequence"):
            add_sequence()
            st.rerun()

    with col2:
        seq_file = st.file_uploader("Choose a file for all the sequences", type="fasta", key="seq_file")
        if seq_file is not None:
            st.session_state.sequences, st.session_state.names_of_sequences = load_fasta_sequences(seq_file)

        distance_matrix_file = st.file_uploader(
            "Choose a file for the distance matrix. The structure should have the names of the sequences at the top and then from the next line the distance matrix",
            type="txt",
            key="distance_matrix_file"
        )
        if distance_matrix_file is not None:
            st.session_state.distance_matrix, st.session_state.names_of_sequences = read_input_for_distance_matrix(distance_matrix_file)

    # validation for the matrix

    dm = st.session_state.distance_matrix
    names = st.session_state.names_of_sequences

    # Check for square matrix
    if len(dm) != len(dm[0]):
        st.error("Distance matrix must be square (same number of rows and columns).")
        st.stop()

    # Check matrix size matches names
    if len(dm) != len(names):
        st.error("Number of sequence names does not match matrix dimensions.")
        st.stop()

    # Check for symmetry
    for i in range(len(dm)):
        for j in range(len(dm)):
            if dm[i][j] != dm[j][i]:
                st.error("Distance matrix must be symmetric.")
                st.stop()

    st.markdown("### Current Sequences")
    st.write(st.session_state.names_of_sequences)
    st.write(st.session_state.sequences)
    st.divider()

    col1, col2, col3 = st.columns(3)
    with col1:
        gap_value = st.number_input("Gap penalty", value=-2.0, step=1.0, format="%.2f")
    with col2:
        match_value = st.number_input("Match reward", value=1.0, step=1.0, format="%.2f")
    with col3:
        mismatch_value = st.number_input("Mismatch penalty", value=-1.0, step=1.0, format="%.2f")

    if st.button("Run MSA Star"):
        sequences_exist = any(seq.strip() for seq in new_sequences)
        matrix_exists = st.session_state.distance_matrix is not None

        if not sequences_exist and not matrix_exists:
            st.error("You must provide either sequences or a distance matrix.")

        elif sequences_exist:
            # Validate input sequences
            valid_chars = set("ACGTURYKMSWBDHVN-")
            invalid_sequences = []
            empty_sequences = []

            for i, seq in enumerate(new_sequences):
                if not seq.strip():
                    empty_sequences.append(i + 1)

            if empty_sequences:
                st.error(f"Sequence(s) {empty_sequences} are empty. Please provide valid input.")
            elif invalid_sequences:
                bad_info = "\n".join([f"Sequence {i}: {s}" for i, s in invalid_sequences])
                st.error(f"The following sequence(s) contain invalid characters:\n{bad_info}")
            else:
                st.success("Running MSA Star with these sequences:")
                st.write(st.session_state.sequences)

                length = len(new_sequences)
                final_score_matrix = np.zeros((length, length))
                for i, j in itertools.combinations(range(length), 2):
                    seq1_clean, seq2_clean = set_sequences(new_sequences[i], new_sequences[j])
                    mat = matrix_building(seq1_clean, seq2_clean, gap_value)
                    score_mat = algorithm(seq1_clean, seq2_clean, mat, gap_value, match_value, mismatch_value)
                    score = score_mat[-1, -1]
                    final_score_matrix[i, j] = score
                    final_score_matrix[j, i] = score

                row_sums = final_score_matrix.sum(axis=1)
                max_row_index = int(np.argmax(row_sums))

                center_seq = new_sequences[max_row_index]
                aligned_center = center_seq
                aligned_others = []
                order = [i for i in range(length) if i != max_row_index]

                for idx in order:
                    seq = new_sequences[idx]
                    seq1, seq2 = set_sequences(aligned_center, seq)
                    seq2 = project_onto_master(seq1, seq2)
                    mat = matrix_building(seq1, seq2, gap_value)
                    score_mat = algorithm(seq1, seq2, mat, gap_value, match_value, mismatch_value)
                    path = traceback(score_mat, seq1, seq2, gap_value, match_value, mismatch_value)
                    new_c, new_s = reconstruct_alignment(seq1, seq2, path)
                    aligned_center, realigned_s, aligned_others = merge_alignment(aligned_center, new_c, new_s, aligned_others)

                msa = [aligned_center] + aligned_others

                identity_percent = calculate_identity_percentage(msa)
                scoring = calculate_msa_score(msa, match_value, mismatch_value, gap_value)
                matches, mismatches, gaps = count_msa_statistics(msa)

                st.markdown("### MSA Score")
                st.write(f"Total Alignment Score: **{scoring}**")
                st.markdown("### Identity Percentage")
                st.write(f"Average Identity: **{identity_percent}%**")
                st.markdown("### Alignment Statistics")
                st.write(f"Matches: **{matches}**, Mismatches: **{mismatches}**, Gaps: **{gaps}**")

                msa_display = "\n".join([f"s{idx + 1}: {' '.join(row)}" for idx, row in enumerate(msa)])
                st.markdown("### Final Multiple Sequence Alignment")
                st.code(msa_display)

                st.markdown("### MSA Heatmap")
                fig = plot_msa_heatmap(msa)
                st.pyplot(fig)

                if st.session_state.distance_matrix is None:
                    st.session_state.distance_matrix = msa_distance_matrix(msa)

                st.markdown("### Distance Matrix")
                st.code(st.session_state.distance_matrix)

                dictionary = {}
                final_cluster = upgma(
                    st.session_state.distance_matrix,
                    len(msa),
                    dictionary,
                    st.session_state.names_of_sequences
                )
                result = ''.join(print_cluster(dictionary, final_cluster)) + ";"
                st.markdown("### Newick Tree Output")
                st.code(result)

                fig = create_image_of_phytree(result)
                st.pyplot(fig)

                report_text = get_text(
                    input_type="MSA",
                    match_value=match_value,
                    mismatch_value=mismatch_value,
                    gap_value=gap_value,
                    sequence_names=st.session_state.names_of_sequences,
                    sequences=st.session_state.sequences,
                    distance_matrix=st.session_state.distance_matrix,
                    newick_tree=result
                )

                report_filename = f"phylo_report_{datetime.datetime.now().strftime('%d-%m-%Y')}.txt"
                save_to_text_file(report_filename, report_text)

                st.divider()
                st.download_button("Download Results", data=report_text, file_name=report_filename)

        elif matrix_exists:
            st.success("Using uploaded distance matrix to build phylogenetic tree.")

            st.markdown("### Distance Matrix")
            st.code(st.session_state.distance_matrix)

            dictionary = {}
            final_cluster = upgma(
                st.session_state.distance_matrix,
                len(st.session_state.distance_matrix),
                dictionary,
                st.session_state.names_of_sequences
            )
            result = ''.join(print_cluster(dictionary, final_cluster)) + ";"
            st.markdown("### Newick Tree Output")
            st.code(result)

            fig = create_image_of_phytree(result)
            st.pyplot(fig)

            report_text = get_text(
                input_type="Distance Matrix",
                match_value=match_value,
                mismatch_value=mismatch_value,
                gap_value=gap_value,
                sequence_names=st.session_state.names_of_sequences,
                sequences=st.session_state.sequences,
                distance_matrix=st.session_state.distance_matrix,
                newick_tree=result
            )

            report_filename = f"phylo_report_{datetime.datetime.now().strftime('%d-%m-%Y')}.txt"
            save_to_text_file(report_filename, report_text)

            st.divider()
            st.download_button("Download Results", data=report_text, file_name=report_filename)