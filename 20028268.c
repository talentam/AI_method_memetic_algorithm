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
static int POP_SIZE = 100;   //please modify these parameters according to your problem
int MAX_NUM_OF_GEN = 10000; //max number of generations
float CROSSOVER_RATE = 0.8;
float MUTATION_RATE = 0.05;
int MATING_POOL_SIZE = 50;
int TOURM_SIZE = 4;

/* declare parameters for simulated annealing here */


/* declare parameteres for tabu search here*/


/* declare parameters for variable neighbourhood search here*/


struct solution_struct best_sln;  //global best solution

//return a random number between 0 and 1
float rand_01()
{
    float number;
    number = (float) rand();
    number = number/RAND_MAX;
    //printf("rand01=%f\n", number);
    return number;
}

//return a random nunber ranging from min to max (inclusive)
int rand_int(int min, int max)
{
    int div = max-min+1;
    int val =rand() % div + min;
    //printf("rand_range= %d \n", val);
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
        for(j=0; j<n; j++)
        {
            my_problems[k]->items[j].dim=dim;
            fscanf(pfile, "%d", &(my_problems[k]->items[j].p)); //read profit data
            //printf("item[j].p=%d\n",my_problems[k]->items[j].p);
        }
        for(i=0; i<dim; i++)
        {
            for(j=0; j<n; j++)
            {
                fscanf(pfile, "%d", &(my_problems[k]->items[j].size[i])); //read size data
                //printf("my_problems[%i]->items[%i].size[%i]=%d\n",k,j,i,my_problems[k]->items[j].size[i]);
            }
        }
        for(i=0; i<dim; i++){
            fscanf(pfile, "%d", &(my_problems[k]->capacities[i]));
            //printf("capacities[i]=%d\n",my_problems[k]->capacities[i] );
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
        /*for(int j=0; j<sln->prob->n; j++)
            fprintf(pfile, "%i ", sln->prob->items[j].p);
        fprintf(pfile, "\n");*/
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
        //val_obj = val;
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
                    //printf("feasb=%i, obj= %f, val=%i\n",temp_sln.feasibility, temp_sln.objective, val_obj);
                    //output_solution(&temp_sln, "my_debug.txt");
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
///////////////////////////////////////////////////////////////////////////


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
        //printf("feasibility%d %d \n", i, pop[i].feasibility);
        // printf("item_num: %d ", prob->n);
        // printf("objective: %f ", pop[i].objective);
        // for(int j=0; j<prob->n; j++) {
        //     printf("%d ", pop[i].x[j]);
        // }
        //printf("%d\n", i);
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

// void selection(int mating_pool[], struct solution_struct* source_pop){
//     //Tournament Selection
//     for(int i = 0; i < MATING_POOL_SIZE; i++){  //50 times
//         int candidate;
//         int candidates_index[TOURM_SIZE];      //storing candidates objective
//         float candidates_obj[TOURM_SIZE];
//         for(int j = 0; j < TOURM_SIZE; j++){
//             int random_index = rand_int(0, source_pop->prob->n-1);
//             candidates_index[j] = random_index;
//             candidates_obj[j] = source_pop[random_index].objective;
//         }
//         int parent = arr_max(candidates_obj);
//         for(int j = 0; j < TOURM_SIZE; j++){
//             if(parent == candidates_obj[j]){
//                 candidate = candidates_index[j];
//                 break;
//             }
//         }
//         mating_pool[i] = candidate;
//     }
    
// }

void selection(struct solution_struct* mating_pool, struct solution_struct* source_pop){
    for(int i = 0; i < MATING_POOL_SIZE; i++){
        ////////////////////////////////
        //防止copy_solution的free导致segmentation fault
        mating_pool[i].x = malloc(sizeof(int)*source_pop->prob->n);
        mating_pool[i].cap_left = malloc(sizeof(int)*source_pop->prob->dim);  


        int candidate;
        for(int j = 0; j < TOURM_SIZE; j++){                        
            int candidates_index[TOURM_SIZE];      //storing candidates objective
            float candidates_obj[TOURM_SIZE];
            for(int j = 0; j < TOURM_SIZE; j++){
                int random_index = rand_int(0, source_pop->prob->n-1);
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
        //printf("%d ", candidate);
        //select which one to be the pool
        copy_solution(&mating_pool[i], &source_pop[candidate]);
        //printf("!!!!!!\n");
        //printf("feasibility%d %d \n", i, mating_pool[i].feasibility);
        //test for whether overcross
        // if(i == 0 || i == 1){
        //     for(int a = 0; a < MATING_POOL_SIZE; a++){
        //         printf("%d ", mating_pool[i].x[a]);
        //     }  
        //     printf("\n");
        // }

        
        //printf("%d\n", i);        
    }
}

void cross_over(struct solution_struct* pop)
{
    //todo
    //1 pop[2i]
    //2 pop[2i+1]
    //for now only crossover between adjacent solution

    for(int i = 0; i < MATING_POOL_SIZE/2; i++){    //i  0-24
        float crossover_rate = rand_01();
        int chorom_length = pop->prob->n;
        //test for whether overcross
        // if(i == 0){
        //     for(int a = 0; a < chorom_length; a++){
        //         printf("%d ", pop[i].x[a]);
        //     }
        //     printf("\n");
        //     for(int a = 0; a < chorom_length; a++){
        //        printf("%d ", pop[i+1].x[a]);
        //     }
        //     printf("\n");
        // }

        if(crossover_rate < CROSSOVER_RATE){
            int split_point = rand_int(1, chorom_length);
            int temp_x[split_point];
            //swap the fist segment
            for(int j = 0; j < split_point; j++){
                temp_x[j] = pop[2*i].x[j];          //store the x value of 1st chorom
                pop[2*i].x[j] = pop[2*i+1].x[j];
                pop[2*i+1].x[j] = temp_x[j];
            }
        }
        evaluate_solution(&pop[2*i]);
        evaluate_solution(&pop[2*i+1]);
        //printf("feasibility%d %d \n", 2*i, pop[2*i].feasibility);
        //printf("feasibility%d %d \n", 2*i+1, pop[2*i+1].feasibility);
        //test for whether overcross
        // if(i == 0){
        //     for(int a = 0; a < chorom_length; a++){
        //         printf("%d ", pop[i].x[a]);
        //     }
        //     printf("\n");
        //     for(int a = 0; a < chorom_length; a++){
        //        printf("%d ", pop[i+1].x[a]);
        //     }
        //     printf("\n");
        // }
    }
}

void mutation(struct solution_struct* pop)
{
    //todo
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
        //printf("feasibility %d \n", pop[i].feasibility);
    }
}
////////////////////////////////////////////////////////////////////////////

//modify the solutions that violate the capacity constraints
void feasibility_repair(struct solution_struct* pop)
{
    //todo
    //sort by price
    int chorom_length = pop->prob->n;
    //printf("chorom_length: %d\n", chorom_length);
    for(int i = 0; i < MATING_POOL_SIZE; i++) {
        //sort first
        //then remove the least one
        //evaluate()
        while(pop[i].feasibility != 1){
            //initialize min_price and min_index
            int min_price, min_index;
            for(int j = 0; j < chorom_length; j++){
                if(pop[i].x[j] == 1){
                    min_price = pop[i].prob->items[0].p;
                    min_index = j;
                    //printf("min_price %d \n", min_price);
                    break;
                }
            }
            //printf("min_price %d \n", min_price);
            //for(int j = chorom_length; j > 0; j--){
            for(int j = 0; j < chorom_length; j++){                
                if(pop[i].prob->items[j].p < min_price && pop[i].x[j] == 1){
                    min_price = pop[i].prob->items[j].p;
                    min_index = j;
                }
            }
            pop[i].x[min_index] = 0;
            evaluate_solution(&pop[i]);
            
        }
        //printf("feasibility %d \n", pop[i].feasibility);
    }
    // for(int a = 0; a < MATING_POOL_SIZE; a++){
    //     printf("feasibility %d \n", pop[a].feasibility);
    // }
    // printf("MATING_POOL_SIZE %d \n", MATING_POOL_SIZE);
}

//local search
void local_search_first_descent(struct solution_struct* pop)
{
    //todo
    //pair swap
    int chorom_length = pop->prob->n;
    //printf("chorom_length: %d\n", chorom_length);
    for(int i = 0; i < MATING_POOL_SIZE; i++){
        //evaluate_solution(&pop[i]);
        //printf("feasibility %d \n", pop[i].feasibility);
        int improvement = 1;
        while(improvement == 1){
            improvement = 0;
            //printf("goin!\n");
            for(int j = 0; j < chorom_length-1; j++){
                for(int k = j+1; k < chorom_length; k++){
                    //printf("j %d k %d \n", j, k);
                    if(pop[i].x[j] + pop[i].x[k] == 1){             //test whether 1/0  or   0/1
                        //swap
                        //printf("feasibility %d \n", pop[i].feasibility);
                        int temp_obj = pop[i].objective;
                        int temp_x = pop[i].x[k];
                        pop[i].x[k] = pop[i].x[j];
                        pop[i].x[j] = temp_x;
                        evaluate_solution(&pop[i]);
                        //printf("feasibility %d \n", pop[i].feasibility);
                        if(pop[i].feasibility <= 0 || pop[i].objective <= temp_obj){      //bad swap
                            //if bad swap, swap it back
                            int temp_x = pop[i].x[k];
                            pop[i].x[k] = pop[i].x[j];
                            pop[i].x[j] = temp_x;
                            evaluate_solution(&pop[i]);
                            //printf("feasibility %d \n", pop[i].feasibility);
                            improvement = 0;
                        }
                        else if(pop[i].feasibility > 0 || pop[i].objective > temp_obj){       //good swap
                            //printf("feasibility %d \n", pop[i].feasibility);
                            improvement = 1;
                            break;
                        }
                    }  
                }
                if(improvement == 1){
                    break;
                }  
            }
            // for(int a = 0; a < chorom_length; a++){
            //    printf("%d ", pop[1].x[a]);
            // }
            // printf("\n");
        }
        //printf("go!\n");
    }
}

void free_population(struct solution_struct* pop, int size) {
    for(int p=0; p<size; p++) {   
        if(pop[p].x != NULL && pop[p].cap_left != NULL) {
        free(pop[p].x);
        free(pop[p].cap_left);
        }
    }
}

void replacement(struct solution_struct* curt_pop, struct solution_struct* new_pop){
    //todo
    //  curt_pop: mating_pool
    //  new_pop: parent_pop
    struct solution_struct mix_pop[POP_SIZE+MATING_POOL_SIZE];
    //struct solution_struct temp_pop;
    //防止copy_solution的free导致segmentation fault?????????????????????????????????????
    for(int i = 0; i < POP_SIZE+MATING_POOL_SIZE; i++){
        mix_pop[i].x = malloc(sizeof(int)*curt_pop->prob->n);
        mix_pop[i].cap_left = malloc(sizeof(int)*curt_pop->prob->dim); 
    }
    //temp_pop.x = malloc(sizeof(int)*curt_pop->prob->n);
    //temp_pop.cap_left = malloc(sizeof(int)*curt_pop->prob->n);
    //add mating_pool individuals to mix_pop
    for(int i = 0; i < MATING_POOL_SIZE; i++) {
        copy_solution(&mix_pop[i], &curt_pop[i]);
    }
    //add parent_pop individuals to mix_pop
    for(int i = 0; i < POP_SIZE; i++) {
        copy_solution(&mix_pop[MATING_POOL_SIZE+i], &new_pop[i]);
    }
    //select the individuals with higher objective      100 times
    for(int i = 0; i < POP_SIZE; i++){  
        int max_obj = mix_pop[0].objective;
        int max_index = 0;
        for(int j = 0; j < POP_SIZE+MATING_POOL_SIZE; j++){
            if(max_obj <= mix_pop[j].objective){
                max_obj = mix_pop[j].objective;
                max_index = j;
                copy_solution(&new_pop[i], &mix_pop[max_index]);
                break;
            }
        }
        //copy_solution(&temp_pop, &mix_pop[max_index]);
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
    
    clock_t time_start, time_fin;
    time_start = clock();
    double time_spent=0;
    int iter =0;        //number of generation
    /////////////////
    struct solution_struct parent_pop[POP_SIZE];
    struct solution_struct mating_pool[MATING_POOL_SIZE];
    //int mating_pool[MATING_POOL_SIZE];       //store the index of mating pool
    //struct solution_struct offspring_pop[POP_SIZE];
    init_population(prob, parent_pop);
    //init_population(prob, offspring_pop);
    /////////////////

    while(iter<MAX_NUM_OF_GEN && time_spent < MAX_TIME)
    {
        //////////////////////
        selection(mating_pool, parent_pop);
        //printf("1");
        cross_over(mating_pool);
        mutation(mating_pool);
        feasibility_repair(mating_pool);
        local_search_first_descent(mating_pool);
        printf("1\n");
        replacement(mating_pool, parent_pop);
        ///////////////////////
        iter++;
        time_fin=clock();
        time_spent = (double)(time_fin-time_start)/CLOCKS_PER_SEC;
    }
    //printf("1\n");
    update_best_solution(parent_pop);
    
    //output_solution(&best_sln[0], "my_debug.txt");
    
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