import sys
import numpy as np
import pandas as pd
import mykmeanssp as spk


def main(original_k: int, goal: str, file_name: str):
    data = pd.read_csv(file_name, header=None)
    data_arr = data.values.tolist()
    original_n, original_dim = len(data_arr), len(data_arr[0])

    if goal == 'spk':
        observations = spk.c_code(original_k, data_arr, 0, original_n, original_dim)  # should be a 2d array
        n = len(observations)
        dim = len(observations[0])
        k = len(observations[0])

        z = 0
        np.random.seed(0)
        first_index = np.random.choice(n)
        first_point = observations[first_index]
        initial_centroid_indices = [first_index]
        centroids = [first_point]
        while z < k - 1:
            probabilities = []
            distances = []
            for i, point in enumerate(observations):
                distance = find_min_distance(centroids, point)
                distances.append(distance)
            sum_of_distances = sum(distances)

            for i, distance in enumerate(distances):
                probabilities.append(distance / sum_of_distances) if sum_of_distances != 0 else probabilities.append(0)

            chosen_centroid_index = np.random.choice(n, size=None, p=probabilities)
            chosen_centroid = observations[chosen_centroid_index]
            initial_centroid_indices.append(chosen_centroid_index)
            centroids.append(chosen_centroid)

            z += 1

        print(*initial_centroid_indices, sep=',')

        spk.c_code_kmeans(initial_centroid_indices, observations, k, n, dim)  # C code
        print("back to python")

    elif goal == 'wam':
        spk.c_code(original_k, data_arr, 1, original_n, original_dim)
    elif goal == 'ddg':
        spk.c_code(original_k, data_arr, 2, original_n, original_dim)
    elif goal == 'lnorm':
        spk.c_code(original_k, data_arr, 3, original_n, original_dim)
    elif goal == 'jacobi':
        spk.c_code(original_k, data_arr, 4, original_n, original_dim)

    print("passed elif's")


def find_min_distance(centroids, vector):
    min_distance = float('inf')
    for i in range(len(centroids)):
        current_distance = calc_distance(centroids[i], vector)
        if min_distance > current_distance:
            min_distance = current_distance
    return min_distance


def calc_distance(point_1, point_2):
    dis = 0
    for i in range(len(point_1)):
        dis += (point_1[i] - point_2[i]) ** 2
    return dis


if __name__ == "__main__":
    goal = sys.argv[2]
    if goal != 'wam' and goal !='ddg' and goal!='lnorm' and goal!='jacobi' and goal!='spk':
        print("Invalid Input!")
        exit(1)
    main(int(sys.argv[1]), sys.argv[2], sys.argv[3])
    print("passed main")
