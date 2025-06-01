import numpy as np
from phytree import *


def read_input(uploaded_file):
    """
        Parse a distance matrix from a file-like UploadedFile object.

        Args:
            uploaded_file (UploadedFile): The uploaded file from Streamlit.

        Returns:
            tuple: (2D distance matrix as list of lists, list of sequence names)
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
    temp = np.array(matrix)
    np.fill_diagonal(temp, np.inf)

    min_value = np.min(temp[np.nonzero(temp)])
    result = np.where(temp == min_value)
    min_index_i, min_index_j = result[0][0], result[1][0]
    return min_index_i, min_index_j


def upgma(matrix, length, dictionary, names_of_sequences):
    leaves = []
    count = 0

    leaves = names_of_sequences.copy()

    number_of_clusters = length

    while (length > 1):

        number_of_clusters = number_of_clusters + 1
        count = count + 1

        min_index_i, min_index_j = minimum_of_matrix(matrix, length)

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