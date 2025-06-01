from io import StringIO
from Bio import SeqIO

# load sequence from fasta file
def load_fasta_sequences(input_file):
    """
        Load and parse sequences from a FASTA file.

        Args:
            input_file (UploadedFile): Uploaded file object from Streamlit's file uploader.

        Returns:
            list[str]: List of sequences as strings.
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
    Generate a formatted report string for phylogenetic tree analysis using UPGMA.

    Returns:
        str: Formatted summary including input data, scoring, and tree.
    """

    text = ""

    # 1. Scoring details
    text += f"""Scoring Parameters:
    - Match value: {match_value}
    - Mismatch value: {mismatch_value}
    - Gap value: {gap_value}
    
    """

    # 2. Input data
    if input_type is "Distance Matrix":
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


    # 4. Tree
    text += f"UPGMA Tree (Newick Format):\n{newick_tree}\n"

    return text



def save_to_text_file(filename, text):
    """
        Save the provided text to a local file.

        Args:
            filename (str): Local file name.
            text (str): Text content to save.

        Returns:
            str: The same text that was saved.
    """

    with open(filename, 'w') as output_file:
        output_file.write(text)

    return text