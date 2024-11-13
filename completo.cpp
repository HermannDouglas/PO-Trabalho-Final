#include <stdio.h> // print, fgets function
#include <string.h> //strlen, memcpy function
#include <math.h> // sqrt, pow, cos, exp, M_PI
#include <time.h> // time function
#include <stdlib.h>

void tweak(double* sol);
void initialize(double* sol);
void hill_climbing(double* sol); // Precisa retornar um double
void simul_annealing(double* sol); // Precisa retornar um double

double ackley(double* sol);
double rosenbrock(double* sol);
double rastrigin(double* sol);
double griewank(double* sol);

double fitness(double* sol){
  return griewank(sol);
}

//min - max values per function
//ackley [-32.768, 32.768]
//griewank [-600, 600]
//rastrigin [-5.12, 5.12]
//rosenbrock [-5, 10]

int dimension=10;
double min=-600;
double max=600;

int main(){
    srand(time(NULL));
double sol[dimension];
    initialize(sol);

    //  Hill Climbing
      //implementar conforme slide

    //  Simulated Annealing
      //implementar conforme slide

    printf("Print solution:\t");
    for(int i=0;i<dimension;i++)
        printf("%f ", sol[i]);
    printf("\nFitness: %f\n",fitness(sol));
    printf("\n------------------\n");
    hill_climbing(sol);
    simul_annealing(sol);
    // # teste de tweak
    // tweak(sol);
    // for(int i=0;i<dimension;i++)
    // printf("%f ", sol[i]);
    // printf("\nFitness: %f",fitness(sol));

  return 0;
}

void hill_climbing(double *sol) {
    int i = 0;
    double *vet_copy = (double *)malloc(dimension * sizeof(double));
    memcpy(vet_copy, sol, dimension * sizeof(double));
    tweak(vet_copy);

    while (i < 5) {
        if (fitness(vet_copy) < fitness(sol)) {
            memcpy(sol, vet_copy, dimension * sizeof(double));
            tweak(vet_copy);
            i++;
        } else {
            tweak(vet_copy);
        }
    }

    free(vet_copy);

    printf("\nNovo vetor sol:\n");
    for (i = 0; i < dimension; i++) {
        printf("%f ", sol[i]);
    }
    printf("\n");
}

void simul_annealing(double *sol){
  double t = 0.95;
  double *vet_best = (double *)malloc(dimension * sizeof(double));
  double *vet_copy = (double *)malloc(dimension * sizeof(double));

  memcpy(vet_best, sol, dimension * sizeof(double));

  while(t > 0.0){
    memcpy(vet_copy, sol, dimension * sizeof(double));
    tweak(vet_copy);
    if(fitness(vet_copy) > fitness(sol) || ((double)rand()/RAND_MAX) < (exp(fitness(vet_copy) - fitness(sol))/t)){ // gerar um número aleatorio entre 0 a 1.
      memcpy(sol, vet_copy, dimension * sizeof(double));
      // printf("fitness: %f\n", fitness(sol));
    }
    t = t- 0.09;

    if(fitness(sol) > fitness(vet_best)){
      // printf("fitness best: %f", fitness(vet_best));
      memcpy(vet_best, sol, dimension * sizeof(double));
    }

  }
  // printf("O best final: %f",fitness(vet_best));
  printf("\nO best final:");
  for (int i = 0; i < dimension; i++) {
    printf("%f ", vet_best[i]);
  }
  free(sol);
  free(vet_copy); 
  free(vet_best);
}

void initialize(double *sol){
  for(int i=0; i<dimension; i++)
    sol[i] = min + ((double)rand()/RAND_MAX) * (max - min);
}

void tweak(double *sol){
  //implementar conforme slide
  double p = 1.0;
  double r = 0.5;

  for(int i=0; i<dimension; i++){
    if(p >= ((double)rand()/RAND_MAX)){
      double n;
      do {
        n = -r + ((double)rand()/RAND_MAX) * (2*r);
      } while (!(min <= sol[i]+n && sol[i]+n <= max));
      sol[i] += n;
    }
  }

}

//adapted from jmetal
double griewank(double *sol){
 double sum = 0.0;
    double mult = 1.0;
    double d = 4000.0;
    for (int var = 0; var < dimension; var++){
      sum += sol[var]*sol[var];
      mult *= cos(sol[var]/sqrt(var+1));
    }

    return (1.0/d) * sum - mult + 1.0;
}

//adapted from jmetal
double rastrigin(double* sol){
    double result = 0.0;
    double a = 10.0;
    double w = 2*M_PI;

    for (int i = 0; i < dimension; i++) {
      result += sol[i]*sol[i] - a*cos(w*sol[i]);
    }
    result += a*dimension;

    return result;
}

//adapted from jmetal
double rosenbrock(double* sol){
    double sum = 0.0;
    for (int i = 0; i < dimension - 1; i++) {
      sum += 100.0 * (sol[i+1]-sol[i]*sol[i])*(sol[i+1]-sol[i]*sol[i]) +(sol[i]-1)*(sol[i]-1) ;
    }
    return sum;
}

//adapted from math function
double ackley(double* sol){
double a=20;
double b=0.2;
double c=2*M_PI;
double sum1 = 0;
double sum2 = 0;
for (int i = 0; i < dimension; i++) {
sum1 += pow(sol[i], 2);
sum2 += cos(c*sol[i]);
}
return -a*exp(-b*sqrt(sum1/dimension)) - exp(sum2/dimension) + a + exp(1.0);
}