import numpy as np
from phytree import *


def read_input():
    matrix = []
    with open('UPGMA_Input.txt') as f:
        matrix = [[float(x) for x in ln.split()] for ln in f]
    matrix = np.asarray(matrix)
    length = len(matrix)
    print("Distance matrix:")
    print(matrix)
    print("")
    print("Number of sequences:", length)
    print("")
    return matrix, length

def minimum_of_matrix(matrix, length):
    min_index_i = 0
    min_index_j = 0
    minimum = float('inf')

    for i in range(length):
        temp = matrix[i]
        min_value = np.min(temp[np.nonzero(temp)])
        j = temp.tolist().index(min_value)

        if min_value < minimum:
            minimum = min_value
            min_index_i = i
            min_index_j = j

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


'''
This is the function used to print the cluster resulting as a result of the upgma method
'''


def print_cluster(dictionary, finalCluster):
    stack = []
    result = []
    stack.append(finalCluster)
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

if __name__ == "__main__":
    matrix, length = read_input()
    dictionary = {}
    finalCluster = upgma(matrix, length, dictionary)
    result = print_cluster(dictionary, finalCluster)
    result=''.join(result)
    result = result+";"
    print(result)

    create_image_of_phytree(result)

