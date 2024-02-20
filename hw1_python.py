import math
import sys
import os


def get_vector(file):
    vector = file.readline().split(',')
    for coordinate in vector:
        if not is_float(coordinate):
            print("An Error Has Occurred")
            return -1
    return [float(coordinate) for coordinate in vector]


def is_float(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


def initialize(input_file, k, d):
    centroids = []
    with open(input_file) as f:
        for i in range(k):
            vector = get_vector(f)
            if vector == -1:
                exit()
            if len(vector) != d:
                print("An Error Has Occurred")
                exit()
            centroids.append(vector)
    return centroids


def average(cluster, d):
    sums = [0] * d
    for vector in cluster:
        for i in range(d):
            sums[i] += vector[i]
    return [coord/len(cluster) for coord in sums]


def distance(u, v):
    return math.sqrt(sum([math.pow(u[i]-v[i], 2) for i in range(len(u))]))


def find_closest_cluster(vector, centroids):
    index = 0
    min_d = distance(vector, centroids[0])
    for i in range(1, len(centroids)):
        d = distance(vector, centroids[i])
        if d < min_d:
            min_d = d
            index = i
    return index


def update_centroids(file, centroids, n, d):
    with open(file) as f:
        clusters = [[] for j in range(len(centroids))]
        for i in range(n):
            vector = get_vector(f)
            if vector == -1 or len(vector) != d:
                print("An Error Has Occurred")
                exit()
            index = find_closest_cluster(vector, centroids)
            clusters[index].append(vector)
    return [average(cluster, d) for cluster in clusters]


def convergence(centroids, new_centroids):
    for i in range(len(centroids)):
        if distance(centroids[i], new_centroids[i]) >= 0.001:
            return False
    return True


def print_vector(vector):
    vector_string = ""
    for coordinate in vector:
        vector_string += ("%.4f" % coordinate) + ","
    print(vector_string[:-1])


def main():
    arguments = sys.argv
    if len(arguments) == 5:
        k, n, d, iterations, input_file = arguments[1], arguments[2], arguments[3], 200, arguments[4]
    elif len(arguments) == 6:
        k, n, d, iterations, input_file = arguments[1], arguments[2], arguments[3], arguments[4], arguments[5]
    else:
        print("An Error Has Occurred")
        exit()

    if not (input_file.endswith('.txt') and os.path.isfile(input_file)):
        print("An Error Has Occurred")
        exit()

    if not n.isdigit():
        print("Invalid number of points!")
        exit()
    n = int(n)

    if not k.isdigit() or int(k) < 2 or int(k) >= n:
        print("Invalid number of clusters!")
        exit()
    k = int(k)

    if not iterations == 200 and (not iterations.isdigit() or int(iterations) <= 1 or int(iterations) >= 1000):
        print("Invalid maximum iteration!")
        exit()
    iterations = int(iterations)

    if not d.isdigit():
        print("Invalid dimension of point!")
        exit()
    d = int(d)

    centroids = initialize(input_file, k, d)

    for i in range(iterations):
        new_centroids = update_centroids(input_file, centroids, n, d)
        if convergence(centroids, new_centroids):
            break
        centroids = new_centroids

    for centroid in centroids:
        print_vector(centroid)


if __name__ == "__main__":
    main()
