/*
Priya Kudva
CIS 630 Project 1, Part 2
*/

#include <iostream>
#include <mpi.h>
#include <string>
#include <cstdlib>
#include <sstream>
#include <ctime>
#include <bits/stdc++.h>

#define SIZE 99999999

using namespace std;

int main(int argc, char *argv[]){

	if(argc != 5) {
		cerr << "Compile with make, Execute with make run" << endl;
		exit(1);
	}

	clock_t total_time = clock();

	// Initialize the MPI environment
	MPI_Init(&argc, &argv);

	// Find out rank, size
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); //proc id
	int num_of_tasks;
	MPI_Comm_size(MPI_COMM_WORLD, &num_of_tasks); // num of processes

	char partition_file[100];
	char graph_file[100];
	FILE* partition_fp;
	FILE* graph_fp;

	int fileSize;
	int* degree;

	int* neighbors; //neighbor info
	int* p_array; // partition info


	strcpy(graph_file, argv[1]);
	strcpy(partition_file, argv[2]);

	int rounds = atoi(argv[3]);
	int partitions = atoi(argv[4]);

	if(partitions != num_of_tasks){
		perror("The number of processes is not the same as the partition number read");
		return 1;
	}
	char* buffer = (char*) malloc(sizeof(char) * 1024 * 1024);

	degree = (int*) malloc(SIZE * sizeof(int));


	p_array = (int*) malloc(SIZE * sizeof(int));

	partition_fp = fopen(partition_file, "r");

	if(partition_fp == NULL){
		perror("Could not open partition file");
		exit(1);
	}

	//find size of file
	fseek(partition_fp, 0L, SEEK_END);
	fileSize = ftell(partition_fp);
	fseek(partition_fp, 0L, SEEK_SET);

	char* nodeID = (char*) malloc(sizeof(char) * 1024);
	char* nodeDegree = (char*) malloc(sizeof(char) * 1024);
	char partitionID;

	int tab_counter = 0;
	int tab_counter1 = 0;
	int tab_counter2 = 0;
	
	int node_length = 0;
	int degree_length = 0;
	int int_node, int_degree, int_partition;
	int buffer_size;
	int array_size = SIZE;

	clock_t start_time = clock();
	clock_t tRead = clock();



	/*------- Read partition file -------*/

	int maxID = 0;
	while(!feof(partition_fp)){
		buffer_size = fread(buffer, 1, 1024*1024, partition_fp);

		for(int i = 0; i < buffer_size; i++){
			if(buffer[i] == '\n'){
				nodeID[node_length] = '\0';
				int_node = atoi(nodeID);
				if(maxID < int_node){
			  		maxID = int_node;
				}

				if(maxID > array_size){ //double the size to fit everything
					degree = (int*) realloc(degree, (maxID+1) * sizeof(int) * 2);
					p_array = (int*) realloc(p_array, (maxID + 1) * sizeof(int) * 2);
					array_size = (maxID + 1) * 2;
				}

				nodeDegree[degree_length] = '\0';
				int_degree = atoi(nodeDegree);

				int_partition = 0 + (partitionID - '0');
				degree[int_node-1] = int_degree; 
				p_array[int_node - 1] = int_partition;

				node_length = 0;
				degree_length = 0;
				tab_counter = 0;
			}
			else if(buffer[i] == ' ' || buffer[i] == '\t'){
				tab_counter++;
			}   
			
			else if(tab_counter == 0){       
				nodeID[node_length] = buffer[i];
				node_length++;
			}     
			else if(tab_counter == 1){
				nodeDegree[degree_length] = buffer[i];
				degree_length++;
			}
			else if(tab_counter == 2){
				partitionID = buffer[i];
			} 
		}
	}
	
	fclose(partition_fp);

	if(maxID < array_size){ //increase size of array if not big enough, not sure if needed?
		degree = (int*) realloc(degree, (maxID + 1) * sizeof(int));
		p_array = (int*) realloc(p_array, (maxID + 1) * sizeof(int));
	}


	/*------- Read graph file -------*/

	graph_fp = fopen(graph_file,"r");

	if(graph_fp == NULL){
		perror("Could not open graph file");
		exit(1);
	}

	//find size of file
	fseek(graph_fp, 0L, SEEK_END);
	fileSize = ftell(graph_fp);
	fseek(graph_fp, 0L, SEEK_SET);

	neighbors = (int*) malloc(sizeof(int) * SIZE);
	int neighbor_size = SIZE;

	int neighbor_index = 0;
	int lines = 0;

	char* node1 = (char*) malloc(sizeof(char) * 1024);
	char* node2 = (char*) malloc(sizeof(char) * 1024);
	int node1_index = 0;
	int node2_index = 0;
	int n1;
	int n2;

	buffer_size = 0;
	tab_counter = 0;
	start_time = clock();

	while(!feof(graph_fp)){
		buffer_size = fread(buffer, 1, 1024*1024, graph_fp);
		for(int i = 0; i < buffer_size; i++){
			if(buffer[i] == '\n'){

				tab_counter = 0;
				node1[node1_index] = '\0'; 
				node2[node2_index] = '\0';
				n1 = atoi(node1);
				n2 = atoi(node2);

				neighbors[neighbor_index] = n1; //assigning nodes to the neighbor array
				neighbors[neighbor_index + 1] = n2;

				neighbor_index = neighbor_index + 2;
				lines++;
		    	
		    	if(lines > neighbor_size){ //double the size to fit everything
					neighbors = (int*) realloc(neighbors, (lines + 1) * sizeof(int) * 2);
					neighbor_size = (lines + 1) * 2;
				}

				node1_index = 0;
				node2_index = 0;
			}
			else if(buffer[i] == ' ' || buffer[i] == '\t'){
				tab_counter = 1;
			}
			else if(tab_counter == 0){
				node1[node1_index] = buffer[i];
				node1_index++;
			}
			else if(tab_counter == 1){
				node2[node2_index] = buffer[i];
				node2_index++;
			}
		}
	}
	fclose(graph_fp);

	printf("time to read input files, partition %d = %.6fsec\n", rank, (double)(clock() - tRead)/CLOCKS_PER_SEC);


	/*------- Calculate Credits per Round -------*/
	double* credit;
	credit = (double *) malloc((maxID) * sizeof(double));
	for(int i = 0; i <= maxID; i++){
		credit[i] = 1;
	}

	double round_time = 0;
	int partition_counter = 0;
	for(int n = 0; n < rounds; n++){
		start_time = clock();
		partition_counter = 1;
		double* send_buffer;
		send_buffer = (double*) malloc((maxID) * sizeof(double));
		for(int i = 0; i <= maxID; i++){
			send_buffer[i] = 0;
		}

		double* receive_buffer;
		receive_buffer = (double*) malloc((maxID) * sizeof(double));
		for(int i = 0; i <= maxID; i++){
			receive_buffer[i] = 0;
		}

		for(int i = 0; i < lines * 2; i = i + 2){
			int f = neighbors[i];
			int j = i + 1;
			int s = neighbors[j];

			if(p_array[f - 1] == rank){
				send_buffer[f - 1] += credit[s - 1] / degree[s - 1];
				send_buffer[s - 1] += credit[f - 1] / degree[f - 1];
			}
		}

		MPI_Allreduce(send_buffer, receive_buffer, maxID, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		credit = receive_buffer;
		ofstream output;
		stringstream fn;
		fn << "round" << n + 1 << ".txt";
		output.open(fn.str());

		start_time = clock();

		for(int i = 0; i <= maxID; i++){
			output << i + 1 << "\t" << degree[i] << "\t" << credit[i] << "\t" << "\n";
		}
		round_time = (double)(clock() - start_time)/CLOCKS_PER_SEC;
   		printf("time for round %d, partition %d = %.6f seconds\n", n + 1, rank, round_time);

   		double total_round_time;
   		MPI_Reduce(&round_time, &total_round_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

   		int total_partition_counter;
   		MPI_Reduce(&partition_counter, &total_partition_counter, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   		if(total_partition_counter == partitions){
   			printf("total time for round %d: %.6f seconds\n", n + 1, total_round_time);
   		}
	}	
	
	MPI_Finalize();
	return 0;
}
