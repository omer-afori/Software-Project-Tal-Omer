import symnmf 
import numpy as np
import sys
import kmeans


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


def average_distance_vector_cluster(vector, cluster):
    return sum([kmeans.distance(vector, v) for v in cluster]) / len(cluster)


def silhouette_coefficients(vectors, i, vector_to_cluster, cluster_to_vectors):
    vector = vectors[i]
    my_cluster_index = vector_to_cluster[i]
    my_cluster = cluster_to_vectors[my_cluster_index]
    a = average_distance_vector_cluster(vector, my_cluster)
    b_list =[]
    for i in range(len(cluster_to_vectors)):
        if i != my_cluster_index:
            b_list.append(average_distance_vector_cluster(vector, cluster_to_vectors[i]))
    b = min(b_list)
    return (b-a)/max(a,b)


def silhouette_score(vectors, vector_to_cluster, cluster_to_vectors):
    return sum([silhouette_coefficients(vectors, i, vector_to_cluster, cluster_to_vectors) for i in range(len(vectors))]) / len(vectors)


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
    #print("symnmf done")
    symnmf_vector_to_cluster, symnmf_cluster_to_vectors = derive_clustering_from_symnmf(symnmf_answer, vectors)
    #print("clusterting for symnmf done")
    kmeans_answer = kmeans.calc_kmeans(file_name, k, n, d, iterations, epsilon)
    #print("kmeans done")
    kmeans_vector_to_cluster, kmeans_cluster_to_vectors = derive_clustering_from_centroids(kmeans_answer, vectors)
    #print("clustering for kmeans done")

    symnmf_score = silhouette_score(vectors, symnmf_vector_to_cluster, symnmf_cluster_to_vectors)
    kmeans_score = silhouette_score(vectors, kmeans_vector_to_cluster, kmeans_cluster_to_vectors)

    print("nmf:", symnmf_score)
    print("kmeans:", kmeans_score)
    

if __name__ == "__main__":
    main(sys.argv)
