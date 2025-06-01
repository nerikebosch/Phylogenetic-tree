import numpy as np
from phytree import *


def read_input_for_distance_matrix(uploaded_file):
    """
    Parses a distance matrix from a file-like UploadedFile object.

    The file is expected to have the sequence names in the first line,
    followed by rows of distance values.

    Args:
        uploaded_file (UploadedFile): The uploaded file from Streamlit.

    Returns:
        tuple: A tuple containing:
            - matrix (list[list[float]]): 2D distance matrix.
            - names (list[str]): List of sequence names.
    """

    matrix = []
    names = []

    content = uploaded_file.getvalue().decode('utf-8').splitlines()

    # First line contains the names
    names = content[0].strip().split()

    # Remaining lines contain the matrix rows
    for line in content[1:]:
        parts = line.strip().split()
        row = list(map(float, parts))
        matrix.append(row)

    return matrix, names

def minimum_of_matrix(matrix):
    """
        Finds the indices of the smallest non-diagonal value in a matrix.

        Args:
            matrix (list[list[float]]): A square distance matrix.

        Returns:
            tuple: Indices (i, j) of the minimum non-diagonal value.
    """

    temp = np.array(matrix)
    np.fill_diagonal(temp, np.inf)

    min_value = np.min(temp[np.nonzero(temp)])
    result = np.where(temp == min_value)
    min_index_i, min_index_j = result[0][0], result[1][0]
    return min_index_i, min_index_j


def upgma(matrix, length, dictionary, names_of_sequences):
    """
        Constructs a UPGMA (Unweighted Pair Group Method with Arithmetic Mean) tree.

        Args:
            matrix (list[list[float]]): Initial distance matrix.
            length (int): Number of sequences (rows/columns in the matrix).
            dictionary (dict): Dictionary to hold cluster tree structure.
            names_of_sequences (list[str]): List of initial sequence names.

        Returns:
            str: The name of the final root cluster node.
    """

    leaves = []
    count = 0

    leaves = names_of_sequences.copy()

    number_of_clusters = length

    while (length > 1):

        number_of_clusters = number_of_clusters + 1
        count = count + 1

        min_index_i, min_index_j = minimum_of_matrix(matrix)

        leaves.append("S" + str(number_of_clusters))
        distance = matrix[min_index_i][min_index_j] / float(2)

        size = 0
        if leaves[min_index_i] not in dictionary.keys():
            size = 1
            distance1 = distance
        else:
            size = dictionary[leaves[min_index_i]][4]
            distance1 = distance - max(dictionary[leaves[min_index_i]][0], dictionary[leaves[min_index_i]][2])

        if leaves[min_index_j] not in dictionary.keys():
            size = size + 1
            distance2 = distance
        else:
            size = size + dictionary[leaves[min_index_j]][4]
            distance2 = distance - max(dictionary[leaves[min_index_j]][0], dictionary[leaves[min_index_j]][2])
        dictionary["S" + str(number_of_clusters)] = [distance1, leaves[min_index_i], distance2, leaves[min_index_j], size]

        # Create a new row and column
        matrix = np.insert(matrix, length, values=float(0), axis=0)
        matrix = np.insert(matrix, length, values=float(0), axis=1)

        for i in range(0, length):
            matrix[-1][i] = matrix[i][-1] = (matrix[i][min_index_i] + matrix[i][min_index_j]) / 2

        # Delete the minimum value
        if min_index_i < min_index_j:
            matrix = np.delete(matrix, min_index_i, 0)
            matrix = np.delete(matrix, min_index_i, 1)
            matrix = np.delete(matrix, (min_index_j) - 1, 0)
            matrix = np.delete(matrix, (min_index_j) - 1, 1)
            length = len(matrix)
            del leaves[min_index_j]
            del leaves[min_index_i]

        else:
            matrix = np.delete(matrix, min_index_i, 0)
            matrix = np.delete(matrix, min_index_i, 1)
            matrix = np.delete(matrix, min_index_j, 0)
            matrix = np.delete(matrix, min_index_j, 1)
            length = len(matrix)
            del leaves[min_index_i]
            del leaves[min_index_j]

    return "S" + str(number_of_clusters)

def print_cluster(dictionary, final_cluster):
    """
        Generates a Newick-style tree string from the cluster dictionary.

        Args:
            dictionary (dict): Dictionary representing the UPGMA cluster structure.
            final_cluster (str): The root node of the final cluster.

        Returns:
            list[str]: List of strings that, when joined, form the Newick format tree.
    """

    stack = []
    result = []
    stack.append(final_cluster)
    while stack:

        current = stack.pop()
        if isinstance(current, float):
            if isinstance(current_prev, float):
                result.pop()
                result.append(")")
            result.append(":" + str(current))
            result.append(",")

        elif current in dictionary.keys():

            stack.append(dictionary[current][0])
            stack.append(dictionary[current][1])
            stack.append(dictionary[current][2])
            stack.append(dictionary[current][3])
            result.append("(")
        else:
            result.append(current)
        current_prev = current

    result.pop()
    result.append(")")
    return result