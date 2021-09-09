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
        print("back to python")
        n = len(observations)
        dim = len(observations[0])
        k = len(observations[0])

        distances = np.full(n, np.inf)  # distances[i] = distance from observations[i] to closest centroid
        centroids = np.full(k, -1)  # 1d array of size k, holds indices of centroids

        np.random.seed(0)
        centroid_index = np.random.choice(n)
        z = 0
        while z < k:
            centroids[z] = centroid_index
            distances = update_distances(distances, centroids, observations, z, n)
            centroid_index = choose_next_centroid(distances, observations)
            z += 1
        initial_centroids = centroids.astype(int).tolist()
        original_indices = [int(initial_centroids[i]) for i in range(k)]  # indices from left column
        print(*original_indices, sep=',')

        spk.c_code_kmeans(initial_centroids, observations.tolist(), k, n, dim)  # C code

    elif goal == 'wam':
        spk.c_code(original_k, data_arr, 1, original_n, original_dim)
    elif goal == 'ddg':
        spk.c_code(original_k, data_arr, 2, original_n, original_dim)
    elif goal == 'lnorm':
        spk.c_code(original_k, data_arr, 3, original_n, original_dim)
    elif goal == 'jacobi':
        spk.c_code(original_k, data_arr, 4, original_n, original_dim)


def update_distances(distances, centroids, observations, z, n):
    for i in range(n):
        distances[i] = find_min_distance(centroids, observations, i, z)
    return distances


def find_min_distance(centroids, observations, obs_index, z):
    min_distance = np.inf
    for i in range(z+1):
        current_distance = calculate_distance(observations[centroids[i]], observations[obs_index])
        if current_distance<min_distance:
            min_distance=current_distance
    return min_distance


def calculate_distance(centroid, observation):
    distance = 0
    for i in range(len(observation)):
        distance+= (observation[i]-centroid[i])**2
    return distance


def choose_next_centroid(distances, observations):
    sum = np.sum(distances)
    probabilities = [(i/sum) for i in distances]
    res = np.random.choice(len(observations), p=probabilities)
    return res


if __name__ == "__main__":
    goal = sys.argv[2]
    if goal != 'wam' and goal !='ddg' and goal!='lnorm' and goal!='jacobi' and goal!='spk':
        print("Invalid Input!")
        exit(1)
    main(int(sys.argv[1]), sys.argv[2], sys.argv[3])
