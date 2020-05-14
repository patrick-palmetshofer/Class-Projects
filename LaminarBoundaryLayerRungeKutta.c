#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//Const parameters for calculation
const int rksteps = 100;		//Number of intervals for 4th order Runge-Kutta scheme
const double tolerance = 1e-4;  //Tolerance with which 3rd boundary condition must be met for the Bisection method to be considered converged
const double eta_max = 4.8;		//Upper integration limit (numerical), used instead of infinity (Book suggests 4.8) 
const double y3start = 10;		//Initial guess for third boundary condition (Second guess for bisection 0)
const double limit = 1e100;     //Limit to avoid overflow of double

void addArraysRK(double* y1, double* y2, double y2coeff, double* yres, int len, double lim)
{
	//Adds the vectors y1 and y2, where y2 is multiplied by y2coeff and writes into yres. All vectors must be of the same length len.
	//yres = y1 + ycoeff*y_2
	//Also checks for whether any element of |y| is greater than a set limit and limits the calculation to this value. 
	//This prevents overflow of the valuable if the function diverges
	for (int i = 0; i < len; i++)
	{
		if (y1[i] >= lim || y2[i] >= lim)
		{
			yres[i] = lim;
			continue;
		}
		if (y2[i] <= -lim || y2[i] <= -lim)
		{
			yres[i] = -lim;
			continue;
		}

		yres[i] = y1[i] + y2coeff * y2[i];
	}
}

void multh(double h, double* y, int len, double lim)
{
	//Multiplies vector y of length len by a scalar h and writes back into y.
	//y *= h
	//Also checks for whether any element of |y| is greater than a set limit and limits the calculation to this value. 
	//This prevents overflow of the valuable if the function diverges
	for (int i = 0; i < len; i++)
	{
		if (y[i] > lim || y[i] < -lim)
			continue;
		y[i] *= h;
	}
}

int RungeKutta4(int* (*f)(double, double*, double*,int), double* y_0, int len, double** y_res, int steps, double x_min, double x_max)
{
	//Generalized 4th order Runge-Kutta solver for systems of equations.
	//y = f(x,y)
	//Takes a function pointer to a function f that takes a scalar x (in this case eta) and a vector y.
	//y_0 is the initial value of the problem, y_res is the result of the function in steps intervals.
	//x_min and x_max define the integration limits

	double h = (x_max - x_min) / steps;
	int yinc = 0;

	//As we don't know the size of our vector at compile-time, dynamically allocate all necessary vectors
	double *yt = malloc(len * sizeof(double)); //Last added row in y, used for convenience
	double *ytemp = malloc(len * sizeof(double)); //Temporary y for calculation in loop

	//Runge-Kutta Parameters
	double *k1 = malloc(len * sizeof(double));
	double *k2 = malloc(len * sizeof(double));
	double *k3 = malloc(len * sizeof(double));
	double *k4 = malloc(len * sizeof(double));

	double x = x_min; //Start at lower x

	memcpy(yt, y_0, len * sizeof(double));
	memcpy(*y_res, y_0, len * sizeof(double)); //Initial/Boundary condition y_0 provides first vector in matrix
	for (int i = 1; i < steps; i++)
	{
		//Calculate x for Runge-Kutta (left side of interval)
		x = x_min + i * h;

		//Calculate k1 through calculating f(x,y) 
		memcpy(ytemp, yt, len * sizeof(double));
		f(x, yt, ytemp,len);
		multh(h, ytemp, len, limit);
		memcpy(k1, ytemp, len * sizeof(double));

		//Calculate k2 through calculating f(x+h/2,y+k1/2) 
		memcpy(ytemp, yt, len * sizeof(double));
		addArraysRK(yt, k1, 0.5, ytemp, len,limit);
		f(x + 0.5 * h, ytemp, ytemp, len);
		multh(h, ytemp, len, limit);
		memcpy(k2, ytemp, len * sizeof(double));

		//Calculate k2 through calculating f(x+h/2,y+k2/2) 
		memcpy(ytemp, yt, len * sizeof(double));
		addArraysRK(yt, k2, 0.5, ytemp, len, limit);
		f(x + 0.5 * h, ytemp, ytemp, len);
		multh(h, ytemp, len, limit);
		memcpy(k3, ytemp, len * sizeof(double));

		//Calculate k2 through calculating f(x+h,y+k3) 
		memcpy(ytemp, yt, len * sizeof(double));
		addArraysRK(yt, k3, 1, ytemp, len, limit);
		f(x + h, ytemp, ytemp, len);
		multh(h, ytemp, len, limit);
		memcpy(k4, ytemp, len * sizeof(double));

		//Calculate y_new = y_old + 1/6 k1 + 1/3 k2 + 1/3 k3 + 1/6 k4
		memcpy(ytemp, yt, len * sizeof(double));
		addArraysRK(ytemp, k1, 1.0 / 6.0, ytemp, len, limit);
		addArraysRK(ytemp, k2, 1.0 / 3.0, ytemp, len, limit);
		addArraysRK(ytemp, k3, 1.0 / 3.0, ytemp, len, limit);
		addArraysRK(ytemp, k4, 1.0 / 6.0, ytemp, len, limit);

		//Write results into yt and y_res (new row)
		memcpy(yt, ytemp, len * sizeof(double));
		memcpy(*(y_res+i), ytemp, len * sizeof(double));
	}

	//Clean up memory
	free(yt);
	free(ytemp);
	free(k1);
	free(k2);
	free(k3);
	free(k4);
	return 0;
}

int fHiemenz(double x, double *yn, double *ynew, int len)
{
	//Function f for y'=f(x,y)
	//Uses the system of equations from the task
	if (len != 3) //Checks for right length of y vector
	{
		fprintf(stderr, "Wrong array length for fHiemenz!\n");
		return -1;
	}

	//Actually calculate f'
	double ret[3];
	ret[0] = yn[1];
	ret[1] = yn[2];
	ret[2] = yn[1] * yn[1] - yn[0] * yn[3] - 1;

	//Write result into output vector ynew
	memcpy(ynew, ret, sizeof(ret));

	return 0;
}

int main()
{
	//Allocate matrix for storing the result of the calculation
	int rkvals = rksteps + 1;
	double **y = malloc(rkvals* sizeof(double*));
	for (int i = 0; i < rkvals; i++)
	{
		*(y+i) = malloc(3 * sizeof(double));
	}

	//Set initial guesses for boundary condition 3. Used in bisection method
	double y3_1 = 0;
	double y3_2 = y3start;

	//Construct vectors from initial guesses. Boundary conditions for y1,y2 are y1=y2=0
	//Calculate solution from Runge-Kutta method
	//Calulate res_1, res_1 which are the values of F'(eta_max) and subtract 1 to find residual
	double y0_1[] = { 0,0,y3_1 };
	RungeKutta4(fHiemenz, y0_1, 3, y, rkvals, 0, eta_max);
	double res_1 = y[rksteps][1] - 1;

	double y0_2[] = { 0,0,y3_2 };
	RungeKutta4(fHiemenz, y0_2, 3, y, rkvals, 0, eta_max);
	double res_2 = y[rksteps][1] - 1;

	//Check for different signs so that there is a root between initial guesses. Otherwise, bisection method is unefined
	if (res_1*res_2 > 0)
	{
		fprintf(stderr, "%s", "Wrong starting Values for Hiemenz function! No sign change!\n");
		return -1;
	}

	//Initialize temporary y-value for bisection method
	double y0_mid[] = { 0,0,0 };
	double res_mid = 1;

	//Bisection method
	int itcount = 0;
	while (res_mid > tolerance)
	{
		itcount++;
		y0_mid[2] = (y0_1[2] + y0_2[2]) / 2; //Get middle point between two values
		RungeKutta4(fHiemenz, y0_mid, 3, y, rkvals, 0, eta_max);
		res_mid = y[rksteps][1] - 1; //Calculate residual for middle point

		//Decide whether to accept lower or higher interval
		if (res_mid*res_1 < 0)
		{
			y0_2[2] = y0_mid[2];
			res_2 = res_mid;
		}
		else
		{
			y0_1[2] = y0_mid[2];
			res_1 = res_mid;
		}

		//Make sure residual is positive (so we can compare against tolerance)
		if (res_mid < 0)
			res_mid *= -1;
	}

	
	//Write result into CSV file
	FILE* filep = fopen("output.csv", "w");
	fprintf(filep, "eta,\tF,\tF',\tF''\n");
	for (int i = 0; i < rkvals; i++)
	{
		fprintf(filep, "%f", eta_max/rksteps*i);
		for (int j = 0; j < 3; j++)
		{
			fprintf(filep, ",\t%f", y[i][j]);
		}
		fprintf(filep,"\n");
	}
	fclose(filep);
	
	//Clean up dynamically allocated memory for y
	for (int i = 0; i < rkvals; i++)
	{
		free(*(y+i));
	}
	free(y);

	//Show iteration count and wait for user input
	printf("Iterations needed: %i\n", itcount);
	getchar();
	return 0;
}