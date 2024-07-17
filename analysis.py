import symnmf 
import numpy as np
import sys
import kmeans
from sklearn.metrics import silhouette_score

iterations = 300
epsilon = 1e-4

def derive_clustering_from_centroids(centroids, vectors):
    vector_to_cluster = []
    cluster_to_vectors = [[] for i in range(len(centroids))]
    for vector in vectors:
        index = kmeans.find_closest_cluster(vector, centroids)
        vector_to_cluster.append(index)
        cluster_to_vectors[index].append(vector)
    return vector_to_cluster, cluster_to_vectors
    

def max_index(arr):
    max_val = arr[0]
    index = 0
    for i in range(1,len(arr)):
        if arr[i] > max_val:
            max_val = arr[i]
            index = i
    return index


def derive_clustering_from_symnmf(answer_matrix, vectors):
    vector_to_cluster = []
    cluster_to_vectors = [[] for i in range(len(answer_matrix[0]))]
    for i in range(len(vectors)):
        index = max_index(answer_matrix[i])
        vector_to_cluster.append(index)
        cluster_to_vectors[index].append(vectors[i])
    return vector_to_cluster, cluster_to_vectors


def main(arguments):
    if (len(arguments) != 3):
        print("An Error Has Occurred")
        return -1
    np.random.seed(0)
    k = int(arguments[1])
    file_name = arguments[2]
    vectors = symnmf.initialize_vectors(file_name)
    n = len(vectors)
    d = len(vectors[0])
    
    symnmf_answer = symnmf.calculate_symnmf(vectors, k, False)
    symnmf_vector_to_cluster, symnmf_cluster_to_vectors = derive_clustering_from_symnmf(symnmf_answer, vectors)
    kmeans_answer = kmeans.calc_kmeans(file_name, k, n, d, iterations, epsilon)
    kmeans_vector_to_cluster, kmeans_cluster_to_vectors = derive_clustering_from_centroids(kmeans_answer, vectors)

    X = np.array(vectors)
    labels = np.array(symnmf_vector_to_cluster)
    symnmf_score = silhouette_score(X, labels)

    X = np.array(vectors)
    labels = np.array(kmeans_vector_to_cluster)
    kmeans_score = silhouette_score(X, labels)

    to_print1 = "nmf: "  + ("%.4f" % symnmf_score)
    to_print2 = "kmeans: "  + ("%.4f" % kmeans_score)

    print(to_print1)
    print(to_print2)

    

if __name__ == "__main__":
    main(sys.argv)
