//
//  main.c
//  mknapsack
//
//  Created by Bai on 14/03/2020.
//  Copyright © 2019 UNNC. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

/* global parameters */
int RAND_SEED[] = {1,20,30,40,50,60,70,80,90,100,110, 120, 130, 140, 150, 160, 170, 180, 190, 200};
int NUM_OF_RUNS = 1;
int MAX_TIME = 30;  //max amount of time permited (in sec)
int num_of_problems;    //number of problems


/* parameters for evlutionary algorithms*/
static int POP_SIZE = 50;   //please modify these parameters according to your problem
int MAX_NUM_OF_GEN = 10000; //max number of generations
float CROSSOVER_RATE = 0.5;
float MUTATION_RATE = 0.01;
int MATING_POOL_SIZE = 50;
int TOURM_SIZE = 2;
clock_t time_start, time_fin;
int K= 2; // k-opt is used
int VNS_SWAP_NUM = 10000;

struct solution_struct best_sln;  //global best solution

//return a random number between 0 and 1
float rand_01()
{
    float number;
    number = (float) rand();
    number = number/RAND_MAX;
    return number;
}

//return a random nunber ranging from min to max (inclusive)
int rand_int(int min, int max)
{
    int div = max-min+1;
    int val =rand() % div + min;
    return val;
}


struct item_struct{
    int dim; //no. of dimensions
    int* size; //volume of item in all dimensions
    int p;
};

struct problem_struct{
    int n; //number of items
    int dim; //number of dimensions
    struct item_struct* items;
    int* capacities;  //knapsack capacities
    int optimal;       //optimal solution
};

void free_problem(struct problem_struct* prob)
{
    if(prob!=NULL)
    {
        if(prob->capacities !=NULL) free(prob->capacities);
        if(prob->items!=NULL)
        {
            for(int j=0; j<prob->n; j++)
            {
                if(prob->items[j].size != NULL)
                    free(prob->items[j].size);
            }
            free(prob->items);
        }
        free(prob);
    }
}

void init_problem(int n, int dim, struct problem_struct** my_prob)
{
    struct problem_struct* new_prob = malloc(sizeof(struct problem_struct));
    new_prob->n=n; new_prob->dim=dim;
    new_prob->items=malloc(sizeof(struct item_struct)*n);
    for(int j=0; j<n; j++)
        new_prob->items[j].size= malloc(sizeof(int)*dim);
    new_prob->capacities = malloc(sizeof(int)*dim);
    *my_prob = new_prob;
}


//example to create problem instances, actual date should come from file
struct problem_struct** load_problems(char* data_file)
{
    int i,j,k;
    //int num_of_probs;
    FILE* pfile = fopen(data_file, "r");
    if(pfile==NULL)
        {printf("Data file %s does not exist. Please check!\n", data_file); exit(2); }
    fscanf (pfile, "%d", &num_of_problems);
 
    struct problem_struct** my_problems = malloc(sizeof(struct problem_struct*)*num_of_problems);
    for(k=0; k<num_of_problems; k++)
    {
        int n, dim, obj_opt;
        fscanf (pfile, "%d", &n);       // number of items
        fscanf (pfile, "%d", &dim);     // nmuber of dimensions
        fscanf (pfile, "%d", &obj_opt); // optimal solution value (zero if unavailable)

        init_problem(n, dim, &my_problems[k]);  //allocate data memory

        my_problems[k]->optimal = obj_opt;

        for(j=0; j<n; j++)
        {
            my_problems[k]->items[j].dim=dim;
            fscanf(pfile, "%d", &(my_problems[k]->items[j].p)); //read profit data
        }
        for(i=0; i<dim; i++)
        {
            for(j=0; j<n; j++)
            {
                fscanf(pfile, "%d", &(my_problems[k]->items[j].size[i])); //read size data
            }
        }
        for(i=0; i<dim; i++){
            fscanf(pfile, "%d", &(my_problems[k]->capacities[i]));
        }
    }
    fclose(pfile); //close file
    return my_problems;
}

struct solution_struct{
    struct problem_struct* prob; //maintain a shallow copy of the problem data
    float objective;    //the price of the current solution
    int feasibility; //indicate the feasiblity of the solution
    int* x; //chromosome vector
    int* cap_left; //capacity left in all dimensions
};

void free_solution(struct solution_struct* sln)
{
    if(sln!=NULL)
    {
        free(sln->x);
        free(sln->cap_left);
        sln->objective=0;
        sln->prob=NULL;
        sln->feasibility=false;
    }
}

//copy a solution from another solution
bool copy_solution(struct solution_struct* dest_sln, struct solution_struct* source_sln)
{
    
    if(source_sln ==NULL) return false;
    if(dest_sln==NULL)
    {   
        dest_sln = malloc(sizeof(struct solution_struct));        
    }
    else{       
        free(dest_sln->cap_left);
        free(dest_sln->x);
    }

    int n = source_sln->prob->n;
    int m =source_sln->prob->dim;
    dest_sln->x = malloc(sizeof(int)*n);
    dest_sln->cap_left=malloc(sizeof(int)*m);
    for(int i=0; i<m; i++)
        dest_sln->cap_left[i]= source_sln->cap_left[i];
    for(int j=0; j<n; j++)
        dest_sln->x[j] = source_sln->x[j];
    dest_sln->prob= source_sln->prob;
    dest_sln->feasibility=source_sln->feasibility;
    dest_sln->objective=source_sln->objective;
    return true;
}



void evaluate_solution(struct solution_struct* sln)
{
    //evaluate the feasibility and objective of the solution
    sln->objective =0; sln->feasibility = 1;
    struct item_struct* items_p = sln->prob->items;
    
    for(int i=0; i< items_p->dim; i++)
    {
        sln->cap_left[i]=sln->prob->capacities[i];
        for(int j=0; j<sln->prob->n; j++)
        {
            sln->cap_left[i] -= items_p[j].size[i]*sln->x[j];
            if(sln->cap_left[i]<0) {
                sln->feasibility = -1*i; //exceeding capacity
                return;
            }
        }
    }

    if(sln->feasibility>0)
    {
        for(int j=0; j<sln->prob->n; j++)
        {
            sln->objective += sln->x[j] * items_p[j].p;
        }
    }
}

//output a given solution to a file
void output_solution(struct solution_struct* sln, char* out_file)
{
    if(out_file !=NULL){
        FILE* pfile = fopen(out_file, "a"); //append solution data
        fprintf(pfile, "%i\n", (int)sln->objective);
        for(int i=0; i<sln->prob->n; i++)
        {
            fprintf(pfile, "%i ", sln->x[i]);
        }
        fprintf(pfile, "\n");
        fclose(pfile);
    }
    else
        printf("sln.feas=%d, sln.obj=%f\n", sln->feasibility, sln->objective);
}


//check the  feasiblity and obj values of solutons from solution_file.
//return 0 is all correct or the index of the first infeasible problem [1, num_of_problems].
int check_solutions(struct problem_struct** my_problems, char* solution_file)
{
    FILE * pfile= fopen(solution_file, "r");
    if(pfile==NULL)
    {
        printf("Solution file %s does not exist. Please check!\n", solution_file);
        exit(2);
    }
    float val_obj;
    int val;
    fscanf (pfile, "%i", &val);
    if(val != num_of_problems)
    {
        printf("The stated number of solutions does not match the number of problems.\n");
        exit(3);
    }
    struct solution_struct temp_sln;
    int count=0, k=0;
    int n, dim;
    while(fscanf (pfile, "%f", &val_obj)!=EOF && k<num_of_problems)
    {
        n= my_problems[k]->n;  dim= my_problems[k]->dim;
        temp_sln.x = malloc(sizeof(int)*n);
        temp_sln.cap_left=malloc(sizeof(int)*dim);
        temp_sln.prob = my_problems[k];
        while(fscanf (pfile, "%i", &val)!=EOF)
        {
            if(val<0 || val>1) {fclose(pfile);  return k+1;} //illeagal values
            temp_sln.x[count] = val;
            count++;
            if(count==n)
            {
                evaluate_solution(&temp_sln);
                if(!temp_sln.feasibility || fabs(temp_sln.objective - val_obj)>0.01)
                {
                    fclose(pfile);
                    return k+1;  //infeasible soltuion or wrong obj
                }
                else{
                    break;
                }
            }
        }
        count=0; k++;
        
        free(temp_sln.x); free(temp_sln.cap_left);
    }
    fclose(pfile);
    return 0;
}

//intialise the population with random solutions
void init_population(struct problem_struct* prob, struct solution_struct* pop){   
    for(int i = 0; i < POP_SIZE; i++){
        pop[i].prob = prob;
        pop[i].x = malloc(sizeof(int)*prob->n);
        pop[i].cap_left = malloc(sizeof(int)*prob->dim);        
        //initialize all items as unpacked
        for(int j = 0; j < prob->n; j++){
            pop[i].x[j] = 0;
        }
        //initialize capacities in all dimensions
        for(int k = 0; k < prob->dim; k++){
            pop[i].cap_left[k]=prob->capacities[k];
        }

        //randomly initialize x        
        while(1){
            int index = rand_int(0, prob->n-1);
            pop[i].x[index] = 1;
            bool space_left = true;
            for(int k = 0; k < prob->dim; k++){
                pop[i].cap_left[k] -= prob->items[index].size[k];
                if(pop[i].cap_left[k] < 0){
                    space_left = false;
                }
            }
            //unpack the last item causing overload
            if(!space_left){   
                pop[i].x[index] = 0;
                for(int k = 0; k < prob->dim; k++){
                    pop[i].cap_left[k] += prob->items[index].size[k];
                }
                break;
            }
        }
        evaluate_solution(&pop[i]);
    }
}

int arr_max(float array[]){
    float max = array[0];
	int size = TOURM_SIZE;
	int i = 0;
	for (i = 0; i < size; i++){
		if (max < array[i])
			max = array[i];
	}
    return max;
}

//tournament selection
void selection(struct solution_struct* mating_pool, struct solution_struct* source_pop){
    for(int i = 0; i < MATING_POOL_SIZE; i++){
        mating_pool[i].x = malloc(sizeof(int)*source_pop->prob->n);
        mating_pool[i].cap_left = malloc(sizeof(int)*source_pop->prob->dim);  

        int candidate;
        for(int j = 0; j < TOURM_SIZE; j++){                        
            int candidates_index[TOURM_SIZE];      //storing candidates objective
            float candidates_obj[TOURM_SIZE];
            for(int j = 0; j < TOURM_SIZE; j++){
                int random_index = rand_int(0, POP_SIZE-1);
                candidates_index[j] = random_index;
                candidates_obj[j] = source_pop[random_index].objective;
            }
            int parent = arr_max(candidates_obj);
            for(int j = 0; j < TOURM_SIZE; j++){
                if(parent == candidates_obj[j]){
                    candidate = candidates_index[j];
                    break;
                }
            }  
        }

        copy_solution(&mating_pool[i], &source_pop[candidate]);      
    }
}

//Uniform Crossover
void cross_over(struct solution_struct* pop)
{
    int chorom_length = pop->prob->n;
    for(int i = 0; i < MATING_POOL_SIZE/2; i++){
        for(int j = 0; j < chorom_length; j++){
            float crossover_rate = rand_01();
            int temp_x;
            if(crossover_rate < CROSSOVER_RATE){
            	// crossover between head solution and tale solution
                temp_x = pop[MATING_POOL_SIZE-1-i].x[j];
                pop[MATING_POOL_SIZE-1-i].x[j] = pop[i].x[j];
                pop[i].x[j] = temp_x;
            }
        }
    }
}

void mutation(struct solution_struct* pop)
{
    int chorom_length = pop->prob->n;
    for(int i = 0; i < MATING_POOL_SIZE; i++){        
        for(int j = 0; j < chorom_length; j++){
            float mutation_rate = rand_01();
            if(mutation_rate < MUTATION_RATE){
                if(pop[i].x[j] == 0){
                    pop[i].x[j] = 1;
                }
                else{
                    pop[i].x[j] = 0;
                }
            }
        }
        evaluate_solution(&pop[i]);
    }
}

//modify the solutions that violate the capacity constraints
//drop by considering the price
void feasibility_repair(struct solution_struct* pop)
{
    int chorom_length = pop->prob->n;
    for(int i = 0; i < MATING_POOL_SIZE; i++) {
        while(pop[i].feasibility != 1){
            //initialize min_price and min_index
            int min_price, min_index;
            for(int j = 0; j < chorom_length; j++){
                if(pop[i].x[j] == 1){
                    min_price = pop[i].prob->items[0].p;
                    min_index = j;
                    break;
                }
            }
            for(int j = 0; j < chorom_length; j++){                
                if(pop[i].prob->items[j].p < min_price && pop[i].x[j] == 1){
                    min_price = pop[i].prob->items[j].p;
                    min_index = j;
                }
            }
            pop[i].x[min_index] = 0;
            evaluate_solution(&pop[i]); 
        }
    }
}

bool can_swap(struct solution_struct* sln, int out, int in)
{
    for(int d =0; d<sln->prob->dim; d++)
    {
        if(sln->cap_left[d]+sln->prob->items[out].size[d] < sln->prob->items[in].size[d])
            return false;
    }
    return true;
}

bool can_move(int nb_indx, int* move, struct solution_struct* curt_sln ){
    bool ret=true;
    if(nb_indx==1){
        ret=can_swap(curt_sln, move[0], move[1]);
    }
    else if(nb_indx==2){//3-item swap
        int i= move[0], j= move[1], k= move[2];
        if(i<0 || j<0 || k<0) return false;
        //1-2 swap
        for(int d=0; d<curt_sln->prob->dim; d++){
            if(curt_sln->cap_left[d] + curt_sln->prob->items[i].size[d] <
                curt_sln->prob->items[j].size[d] + curt_sln->prob->items[k].size[d])
                return false;
        }        
    }
    else ret=false;
    return ret;
}

bool apply_move(int nb_indx, int* move, struct solution_struct* sln ){
    bool ret=true;
    if(nb_indx==1){
        for(int d=0; d<sln->prob->dim; d++){
            sln->cap_left[d] = sln->cap_left[d] + sln->prob->items[move[0]].size[d]-
                sln->prob->items[move[1]].size[d];
        }
        sln->objective += sln->prob->items[move[1]].p-sln->prob->items[move[0]].p;
        sln->x[move[0]]=0; sln->x[move[1]]=1;
    }
    else if(nb_indx==2){//3-item swap
        int i= move[0], j= move[1], k= move[2];
        if(i<0 || j<0 || k<0) return false;
        //1-2 swap
        for(int d=0; d<sln->prob->dim; d++){
            sln->cap_left[d] = sln->cap_left[d]+sln->prob->items[i].size[d] -
                sln->prob->items[j].size[d] - sln->prob->items[k].size[d];
        }
        sln->objective += sln->prob->items[j].p+sln->prob->items[k].p-sln->prob->items[i].p;
        sln->x[i]=0; sln->x[j]=1; sln->x[k]=1;

        
    }
    else ret=false;
    return ret;
}

//two neibourhood: 1-1 swap and 1-2 swap
struct solution_struct* best_descent_vns(int nb_indx, struct solution_struct* curt_sln){
    struct solution_struct* best_neighb = malloc(sizeof(struct solution_struct));
    best_neighb->cap_left = malloc(sizeof(int)*curt_sln->prob->dim);
    best_neighb->x = malloc(sizeof(int)*curt_sln->prob->n);
    copy_solution(best_neighb, curt_sln);
    int n=curt_sln->prob->n;
    int curt_move[] ={-1,-1,-1}, best_move []={-1,-1,-1}, delta=0, best_delta=0;  //storing best neighbourhood moves
    //used for spliting the packed and unpacked items into two lists
    int packed[n], unpacked[n];
    int packed_index = 0;
    int unpacked_index = 0;

    switch (nb_indx)
    {
        case 1:
        	//pair swap
            //divide the items into 2 groups to reduce the time complexity
            packed_index = 0; unpacked_index = 0;
            for(int i = 0; i < n; i++){
                if(curt_sln->x[i] == 1){
                    packed[packed_index] = i;
                    packed_index++;
                }
                else{
                    unpacked[unpacked_index] = i;
                    unpacked_index++;
                }
            }

            for(int i=0; i<packed_index; i++){
                for(int j=0; j<unpacked_index; j++){
                    curt_move[0]= packed[i]; curt_move[1]= unpacked[j]; curt_move[2]=-1;
                    if(can_move(nb_indx, &curt_move[0], best_neighb)){
                        delta = curt_sln->prob->items[unpacked[j]].p -curt_sln->prob->items[packed[i]].p;
                        if(delta > best_delta){
                            best_delta = delta; best_move[0] = packed[i]; best_move[1] = unpacked[j]; best_move[2]=-1;
                        }
                    }
                }
            }
            if(best_delta>0) { apply_move(nb_indx, &best_move[0], best_neighb);}
            break;
        case 2:
            // 1-2 swap
            //divide the items into 2 groups to reduce the time complexity
            packed_index = 0; unpacked_index = 0;
            for(int i = 0; i < n; i++){
                if(curt_sln->x[i] == 1){
                    packed[packed_index] = i;
                    packed_index++;
                }
                else{
                    unpacked[unpacked_index] = i;
                    unpacked_index++;
                }
            }
            
            //select items to 1-2 swap
            int iteration = 0;
            while(iteration < VNS_SWAP_NUM){
                int i = rand_int(0, packed_index-1);
                int j = rand_int(0, unpacked_index-1);
                int k = rand_int(0, unpacked_index-1);
                while(j == k){
                    k = rand_int(0, unpacked_index-1);
                }
                curt_move[0]=packed[i]; curt_move[1]=unpacked[j]; curt_move[2]=unpacked[k];
                if(can_move(nb_indx, &curt_move[0], best_neighb)){
                    delta = curt_sln->prob->items[unpacked[k]].p +curt_sln->prob->items[unpacked[j]].p-curt_sln->prob->items[packed[i]].p;
                    if(delta > best_delta){
                        best_delta = delta; best_move[0] = packed[i]; best_move[1] = unpacked[j]; best_move[2] = unpacked[k];
                    }
                }
                iteration++;
            }           
            if(best_delta>0) { apply_move(nb_indx, &best_move[0], best_neighb);}
            break;
            
        default:
            printf("Neighbourhood index is out of the bounds, nothing is done!\n");
    }
    return best_neighb;
}

//best descent VNS
void varaible_neighbourhood_search(struct solution_struct* pop){
    for(int i = 0; i < MATING_POOL_SIZE; i++){
        float ratio = rand_01();
        if(ratio >= 1){
            continue;
        }
        int nb_indx = 0; //neighbourhood index
        time_fin=clock();
        double time_spent = (double)(time_fin-time_start)/CLOCKS_PER_SEC;
        struct solution_struct* curt_sln = &pop[i];

        while(time_spent < MAX_TIME && nb_indx<K){      
            struct solution_struct* neighb_s = best_descent_vns(nb_indx+1, curt_sln); //best solution in neighbourhood nb_indx
            if(neighb_s->objective > curt_sln->objective){
                copy_solution(curt_sln, neighb_s);
                nb_indx=1;
            }
            else nb_indx++;
            free_solution(neighb_s);free(neighb_s); 
        }
        pop[i] = *curt_sln;
    } 
}

//free population
void free_population(struct solution_struct* pop, int size) {
    for(int p=0; p<size; p++) {   
        if(pop[p].x != NULL && pop[p].cap_left != NULL) {
        free(pop[p].x);
        free(pop[p].cap_left);
        }
    }
}

int cmpfunc_sln (const void * a, const void * b) {
    const struct solution_struct* sln1 = a;
    const struct solution_struct* sln2 = b;
    if(sln1->objective > sln2 ->objective) return -1;
    if(sln1->objective < sln2 ->objective) return 1;
    return 0;
}

void replacement(struct solution_struct* curt_pop, struct solution_struct* new_pop){
    struct solution_struct mix_pop[POP_SIZE+MATING_POOL_SIZE];
    for(int i = 0; i < POP_SIZE+MATING_POOL_SIZE; i++){
        mix_pop[i].x = malloc(sizeof(int)*curt_pop->prob->n);
        mix_pop[i].cap_left = malloc(sizeof(int)*curt_pop->prob->dim); 
    }
    //add mating_pool individuals to mix_pop
    for(int i = 0; i < MATING_POOL_SIZE; i++) {
        copy_solution(&mix_pop[i], &curt_pop[i]);
    }
    //add parent_pop individuals to mix_pop
    for(int i = 0; i < POP_SIZE; i++) {
        copy_solution(&mix_pop[MATING_POOL_SIZE+i], &new_pop[i]);
    }
    //qsort the list
    qsort(mix_pop, POP_SIZE+MATING_POOL_SIZE, sizeof(struct solution_struct), cmpfunc_sln);
    for(int i = 0; i < POP_SIZE; i++){  
        copy_solution(&new_pop[i], &mix_pop[i]);
    }
    //free
    free_population(mix_pop, POP_SIZE+MATING_POOL_SIZE);
    free_population(curt_pop, MATING_POOL_SIZE);

}

//update global best solution from sln
void update_best_solution(struct solution_struct* sln)
{
    if(best_sln.objective < sln->objective)
    copy_solution(&best_sln, sln);
}

//memetic algorithm
int memeticAlgorithm(struct problem_struct* prob)
{
    time_start = clock();
    double time_spent=0;
    int iter =0;        //number of generation
    struct solution_struct parent_pop[POP_SIZE];
    struct solution_struct mating_pool[MATING_POOL_SIZE];
    init_population(prob, parent_pop);

    while(iter<MAX_NUM_OF_GEN && time_spent < MAX_TIME)
    {
        selection(mating_pool, parent_pop);
        cross_over(mating_pool);
        mutation(mating_pool);
        feasibility_repair(mating_pool);
        varaible_neighbourhood_search(mating_pool); 
        replacement(mating_pool, parent_pop);
        iter++;
        time_fin=clock();
        time_spent = (double)(time_fin-time_start)/CLOCKS_PER_SEC;
        //printf("gap: %f  worst: %f   iter: %d  time: %f\n", (parent_pop->prob->optimal-parent_pop[0].objective)/parent_pop->prob->optimal, (parent_pop->prob->optimal-parent_pop[POP_SIZE-1].objective)/parent_pop->prob->optimal, iter, time_spent);
        //printf("best: %f  worst: %f   iter: %d  time: %f\n", parent_pop[0].objective, parent_pop[POP_SIZE-1].objective, iter, time_spent);   
    }
    update_best_solution(parent_pop);
    //printf("optimal: %d\n", parent_pop->prob->optimal);
    //printf("gap: %f  worst: %f   iter: %d  time: %f\n", (parent_pop->prob->optimal-parent_pop[0].objective)/parent_pop->prob->optimal, (parent_pop->prob->optimal-parent_pop[POP_SIZE-1].objective)/parent_pop->prob->optimal, iter, time_spent);    
    printf("best: %f  worst: %f   iter: %d  time: %f\n", parent_pop[0].objective, parent_pop[POP_SIZE-1].objective, iter, time_spent);   
    return 0;
}


int main(int argc, const char * argv[]) {
    
    printf("Starting the run!\n");
    char data_file[50]={"somefile"}, out_file[50]={}, solution_file[50]={};  //max 50 problem instances per run
    if(argc<3)
    {
        printf("Insufficient arguments. Please use the following options:\n   -s data_file (compulsory)\n   -o out_file (default my_solutions.txt)\n   -c solution_file_to_check\n   -t max_time (in sec)\n");
        return 1;
    }
    else if(argc>9)
    {
        printf("Too many arguments.\n");
        return 2;
    }
    else
    {
        for(int i=1; i<argc; i=i+2)
        {
            if(strcmp(argv[i],"-s")==0)
                strcpy(data_file, argv[i+1]);
            else if(strcmp(argv[i],"-o")==0)
                strcpy(out_file, argv[i+1]);
            else if(strcmp(argv[i],"-c")==0)
                strcpy(solution_file, argv[i+1]);
            else if(strcmp(argv[i],"-t")==0)
                MAX_TIME = atoi(argv[i+1]);
        }
        //printf("data_file= %s, output_file= %s, sln_file=%s, max_time=%d", data_file, out_file, solution_file, MAX_TIME);
    }
    struct problem_struct** my_problems = load_problems(data_file);
    
    if(strlen(solution_file)<=0)
    {
        if(strcmp(out_file,"")==0) strcpy(out_file, "my_solutions.txt"); //default output
        FILE* pfile = fopen(out_file, "w"); //open a new file
        fprintf(pfile, "%d\n", num_of_problems); fclose(pfile);
        for(int k=0; k<num_of_problems; k++)
        {
            best_sln.objective=0; best_sln.feasibility=0;
            for(int run=0; run<NUM_OF_RUNS; run++)
            {
                srand(RAND_SEED[run]);
                
                memeticAlgorithm(my_problems[k]); // call MA method
            }
            output_solution(&best_sln,out_file);
        }
    }
    for(int k=0; k<num_of_problems; k++)
    {
       free_problem(my_problems[k]); //free problem data memory
    }
    free(my_problems); //free problems array
    if(best_sln.x!=NULL && best_sln.cap_left!=NULL){ free(best_sln.cap_left); free(best_sln.x);} //free global
    return 0;
}
