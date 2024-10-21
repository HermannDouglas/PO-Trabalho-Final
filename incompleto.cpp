#include <stdio.h> // print, fgets function
#include <string.h> //strlen, memcpy function
#include <math.h> // sqrt, pow, cos, exp, M_PI
#include <time.h> // time function

void tweak(double* sol);
void initialize(double* sol);

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
    printf("Fitness: %f\n",fitness(sol));

	return 0;
}
void initialize(double *sol){
  //implementar conforme slide
}
void tweak(double *sol){
  //implementar conforme slide
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

    return 1.0/d * sum - mult + 1.0;
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
