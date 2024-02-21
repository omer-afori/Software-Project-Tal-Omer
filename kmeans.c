#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define BUFF_SIZE 200
#define error_msg ("An Error Has Occurred\n")
#define true (1)
#define false (0)
#define EPSILON 0.001

typedef struct node
{
	double* vec;
	struct node* next;
}node_t;

node_t* head = 0;
int k = 0;
int dimension = 0;
int vector_count = 1;


void free_list_with_content(node_t* starting_node){
	node_t* current_node = starting_node;
	if (!current_node->next){
		free(current_node->vec);
		free(current_node);
		return;
	}
	free_list_with_content(current_node->next);
	free(current_node->vec);
	free(current_node);
}


double distance(double* vec1, double* vec2){
	/*assumes both vectors are of the same length, len*/
	double sum = 0;
	int i = 0;
	for (i = 0 ; i < dimension; i++){
		sum += pow((vec1[i]-vec2[i]),2);
	}
	return sqrt(sum);
}


int count_char(char* string, char ch){
	int count = 0;
	size_t i = 0;
	for(i = 0; i < strlen(string); i++){
		if (string[i] == ch){
			count++;
		}
	}
	return count;
}

/*void print_vec(double* vec){
	int i;
	for (i = 0; i < dimension-1; i++){
		printf("%f,",vec[i]);
	}
	printf("%f\n", vec[dimension-1]);
}*/

void print_vec(double* vec){
	int i;
	for (i = 0; i < dimension-1; i++){
		printf("%.4f,",vec[i]);
	}
	printf("%.4f\n", vec[dimension-1]);
}

double* get_vector(){
	char *lineptr;
	size_t size;
	int linelen;
	int index = 0;
	double* vector;
	int counter = 0;
	int i = 0;
	linelen = getline(&lineptr, &size, stdin);
	if(linelen < 0){
		free(lineptr);
		printf("in get vector\n");
		printf(error_msg);
		exit(1);
	}
	for (i = 0; i < linelen; i++){
		if (lineptr[i] == ',' || lineptr[i] == '\n'){
			lineptr[i] = '\0';
		}
	}
	
	vector = malloc(sizeof(double)*dimension);
	if (!vector){
		/*FREE*/
		printf("in get vector\n");
		printf(error_msg);
		exit(1);
	}
	for(i = 0; i < linelen; i++){
		if (lineptr[i] == '\0'){
			vector[counter] = atof(lineptr + index*sizeof(char));
			index = i+1;
			/*printf("%f ", vector[counter]);*/
			counter++;
		}
	}
	/*
	printf("\nk,d,n are: %d %d %d\n", k, dimension, vector_count);
	for(i = 0; i < 3; i++){
		printf("%f ", vector[i]);
	}
	printf("\n");
	print_vec(vector);*/
	return vector;
}

void set_globals(){
	int i = 0;
	double* cur_vec;
	node_t* current_node;
	node_t* prev_node ;
	head = (node_t*)malloc(sizeof(node_t));
	current_node = head;
	if (!head){
		printf(error_msg);
		exit(1);
	}
	for(i = 0; i < vector_count; i++){
		/*at each point, assumes current_node is allocated already*/
		cur_vec = get_vector();
		current_node->vec = cur_vec;
		current_node->next = (node_t*)malloc(sizeof(node_t));
		if (!current_node->next){
			printf(error_msg);
			exit(1);
		}
		prev_node = current_node;
		current_node = current_node->next;
	}
	free(prev_node->next); /*at the end, current node is empty*/
	prev_node->next = 0;
	return;
}
double* average_cluster(node_t* cluster){
	double* centroid = malloc(sizeof(double)*dimension);
	int i = 0;
	int j = 0;
	int cluster_len = 0;
	node_t* current_node = cluster;
	if (!centroid){
		printf("in average_cluster\n");
		printf(error_msg);
		/*FREE*/
		exit(1);
	}
	for (i = 0; i < dimension; i++){
		centroid[i] = 0.0;
	}
	while(current_node) {
		cluster_len += 1;
		for (j = 0; j < dimension; j++){
			centroid[j] += current_node->vec[j];
		}
		current_node = current_node->next;
	}
	for (i = 0; i < dimension; i++){
		centroid[i] = centroid[i]/cluster_len;
	}
	return centroid;
}
int find_closest_cluster(double* vec, double** centroids){
	double min_d = distance(vec, centroids[0]);
	double d;
	int min_index = 0;
	int i = 0;
	for (i = 1; i < k; i++){
		d = distance(vec, centroids[i]);
		if (d < min_d){
			min_d = d;
			min_index = i;
		}
	} 
	return min_index;
}
double** update_centroids(double** centroids){
	node_t* current_node = head;
	double** new_centroids;
	double* cur_vec;
	int i;
	int j;
	int index;
	int* lengths;
	lengths = malloc(sizeof(int)*k);
	new_centroids = malloc(sizeof(double*)*k);
	if (!new_centroids || !lengths){
		printf(error_msg);
		exit(1);
	}
	for(i = 0; i < k; i++){
		new_centroids[i] = malloc(sizeof(double)*dimension);
		if (!new_centroids){
			printf(error_msg);
			exit(1);
		}
		for (j = 0; j < dimension; j++){
			new_centroids[i][j] = 0;
		}
		lengths[i] = 0;
	}

	while(current_node){
		cur_vec = current_node->vec;
		index = find_closest_cluster(cur_vec, centroids);
		lengths[index] += 1;
		for(i = 0; i < dimension; i++){
			new_centroids[index][i] += cur_vec[i];
		}
		current_node = current_node->next;
	}
	/*DEBUG*/
	/*
	for(i = 0; i < k; i++){
		for(j = 0; j < dimension; j++){
			printf("%f ", new_centroids[i][j]);
		}
		printf("\n");
	}
	*/
	for(i = 0; i < k; i++){
		for(j = 0; j < dimension; j++){
			new_centroids[i][j] = new_centroids[i][j]/lengths[i];
		}
	}
	free(lengths);

	return new_centroids;
}

/*
double** update_centroids(double** centroids){
	node_t** clusters;
	node_t** last_nodes;
	double* cur_vec;
	int index;
	double** new_centroids;
	int i = 0;
	node_t* current_node = head;

	last_nodes = malloc(sizeof(node_t*)*k);
	clusters = malloc(sizeof(node_t*)*k);
	if (!clusters || !last_nodes){
		printf("update_centroids - 1\n");
		printf(error_msg);
		exit(1);
	}
	new_centroids = malloc(sizeof(double*)*k);
	if (!new_centroids){
		printf("update_centroids - 4\n");
		printf(error_msg);
		exit(1);
	}
	for (i = 0; i < k; i++){
		clusters[i] = (node_t*)malloc(sizeof(node_t));
		if(!clusters[i]){
			printf("update_centroids - 3\n");
			printf(error_msg);
			exit(1);
		}
	}
	while(current_node){
		cur_vec = current_node->vec;
		index = find_closest_cluster(vec, centroids);



		if (!last_nodes[index]){
			last_nodes[index] = (node_t*)malloc(sizeof(node_t));
		}
		else{
			last_nodes[index]->next = (node_t*)malloc(sizeof(node_t));
			last_nodes[index] = last_nodes[index]->next;
		}
		if(!last_nodes[index]){
			printf(error_msg);
			exit(1);
		}
		last_nodes[index]->vec = cur_vec;
		last_nodes[index]->next = 0;
	}

	
	for(i = 0; i < k; i++){
		new_centroids[i] = average_cluster(clusters[i]);
	}
	for(i = 0; i < k; i++){
		free_list_without_content(clusters[i]);
	}
	free(clusters);
	free(last_nodes);
	return new_centroids;
}
*/

double* copy_vec(double* vec){
	double* ans_vec;
	int i;
	ans_vec = malloc(sizeof(double)*dimension);
	if (!ans_vec){
		printf(error_msg);
		exit(1);
	}
	for(i = 0; i < dimension; i++){
		ans_vec[i] = vec[i];
	}
	return ans_vec;
}

double** initialize_centroids(){
	double** centroids = malloc(sizeof(double*)*k);
	node_t* current_node = head;
	double* vec;
	int i = 0;
	if (!centroids){
		/*FREE*/
		printf("init_centroids\n");
		printf(error_msg);
		exit(1);
	}
	for (i = 0; i < k; i++){
		if (!current_node){
			printf(error_msg);
			exit(1);
		}
		vec = copy_vec(current_node->vec);
		current_node = current_node->next;
		centroids[i] = vec;
	}
	return centroids;
}


int convergence(double** centroids, double** new_centroids){
	int i;
	for (i = 0; i < k; i++){
		if (distance(centroids[i], new_centroids[i]) >= EPSILON){
			return false;
		}
	}
	return true;
}

void free_centroids(double** centroids){
	int j;
	for (j = 0; j < k; j++){
		free(centroids[j]);
	}
	free(centroids);
}

void free_globals(){
	free_list_with_content(head);
}

int main(int argc, char* argv[]){
	int iter;
	double** centroids;
	double** new_centroids;
	int i = 0;
	if (argc < 4 || argc > 5){
		printf(error_msg);
		exit(1);
	}
	if (argc == 4){
		iter = 200;
		k = atoi(argv[1]);
		if (k <= 0){
			printf("Invalid number of clusters!\n");
			exit(1);
		}
		vector_count = atoi(argv[2]);
		if (vector_count <= 0){
			printf("Invalid number of points!\n");
			exit(1);
		}
		dimension = atoi(argv[3]);
		if (dimension <= 0){
			printf("Invalid dimension of point!");
		}
	}

	if (argc == 5){
		k = atoi(argv[1]);
		if (k <= 0){
			printf("Invalid number of clusters!\n");
			exit(1);
		}
		vector_count = atoi(argv[2]);
		if (vector_count <= 0){
			printf("Invalid number of points!\n");
			exit(1);
		}
		dimension = atoi(argv[3]);
		if (dimension <= 0){
			printf("Invalid dimension of point!");
		}
		iter = atoi(argv[4]);
		if (iter <= 0){
			printf("Invalid maximum iteration!\n");
			exit(1);
		}
	}
	set_globals();
	if (k < 2 || k >= vector_count){
		printf("Invalid number of clusters!\n");
		exit(1);
	}
	if (iter < 2 || iter > 999){
		printf("Invalid maximum iteration!\n");
		exit(1);
	}
	centroids = initialize_centroids();

	for (i = 0; i < iter; i++){
		/*printf("round num: %d\n", i);*/

		new_centroids = update_centroids(centroids);
		/*
		for(j = 0; j < k; j++){
			print_vec(new_centroids[j]);
		}
		printf("\n");*/
		if (convergence(centroids, new_centroids)){
			free_centroids(new_centroids);
			break;
		}
		free_centroids(centroids);
		centroids = new_centroids;
	}

	for(i = 0; i < k; i++){
		print_vec(centroids[i]);
	}
	free_centroids(centroids);
	free_globals();
	return 0;
}


