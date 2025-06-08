from io import StringIO
from Bio import SeqIO

# load sequence from fasta file
def load_fasta_sequences(input_file):
    """
    Loads and parses sequences from a FASTA file.

    Args:
        input_file (UploadedFile): Uploaded file object from Streamlit's file uploader.

    Returns:
        tuple:
            - fasta_sequences (list[str]): List of sequence strings.
            - names_of_sequences (list[str]): Corresponding sequence names.
    """

    fasta_sequences = []
    names_of_sequences = []

    stringio = StringIO(input_file.getvalue().decode('utf-8'))
    for sequence in SeqIO.parse(stringio, "fasta"):
        fasta_sequences.append(str(sequence.seq))
        names_of_sequences.append(str(sequence.name))

    return fasta_sequences, names_of_sequences


def get_text(input_type, match_value, mismatch_value, gap_value,
             sequence_names=None, sequences=None, distance_matrix=None, newick_tree=None):
    """
    Generates a formatted text report for phylogenetic analysis using UPGMA.

    Args:
        input_type (str): The type of input used ("Distance Matrix" or "Raw Sequences").
        match_value (int): Score for a character match.
        mismatch_value (int): Penalty for a mismatch.
        gap_value (int): Penalty for a gap.
        sequence_names (list[str], optional): List of sequence names.
        sequences (list[str], optional): List of sequences (if applicable).
        distance_matrix (list[list[float]], optional): Distance matrix (if applicable).
        newick_tree (str, optional): Newick format tree string.

    Returns:
        str: A formatted text summary including input data, scoring values, and the UPGMA tree.
    """

    text = ""

    text += f"""Scoring Parameters:
    - Match value: {match_value}
    - Mismatch value: {mismatch_value}
    - Gap value: {gap_value}
    
    """

    if input_type == "Distance Matrix":
        text += "Input Type: Distance Matrix\n\n"
        text += "Sequence Names:\n"
        for name in sequence_names:
            text += f"- {name}\n"
        text += "\nDistance Matrix:\n"
        for row in distance_matrix:
            text += "\t".join(map(str, row)) + "\n"
    else:
        text += "Input Type: Raw Sequences\n\n"
        text += "Sequences:\n"
        for name, seq in zip(sequence_names, sequences):
            text += f"{name}: {seq}\n\n"


    text += f"UPGMA Tree (Newick Format):\n{newick_tree}\n"

    return text



def save_to_text_file(filename, text):
    """
        Saves a given text string to a local file.

        Args:
            filename (str): The name or path of the file to be created.
            text (str): The textual content to write into the file.

        Returns:
            str: The same text content that was saved to the file.
    """

    with open(filename, 'w') as output_file:
        output_file.write(text)

    return text