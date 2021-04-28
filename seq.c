#include <omp.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>

#define f(a,b) for(int a =0; a<b; a++)

void seq_crout(double **A, double **L, double **U, int n) {
	int i, j, k;
	double sum = 0;

	for (i = 0; i < n; i++) {
		U[i][i] = 1;
	}

	for (j = 0; j < n; j++) {
		for (i = j; i < n; i++) {
			sum = 0;
			for (k = 0; k < j; k++) {
				sum = sum + L[i][k] * U[k][j];	
			}
			L[i][j] = A[i][j] - sum;
		}

		for (i = j; i < n; i++) {
			sum = 0;
			for(k = 0; k < j; k++) {
				sum = sum + L[j][k] * U[k][i];
			}
			if (L[j][j] == 0) {				
				printf("det(L) close to 0!\n Can't divide by 0...\n");
				exit(EXIT_FAILURE);
			}
			U[j][i] = (A[j][i] - sum) / L[j][j];
		}
	}
}

void crout_1(double **A, double **L, double **U, int n, int t) {
	int i, j, k;
	double sum = 0;

	
	for (i = 0; i < n; i++) {
		U[i][i] = 1;
	}
	
	for (j = 0; j < n; j++) {
		#pragma omp parallel for shared(j,L,A,U,n) private(i,k,sum)
		for (i = j; i < n; i++) {
			sum = 0;
			for (k = 0; k < j; k++) {
				sum = sum + L[i][k] * U[k][j];	
			}
			L[i][j] = A[i][j] - sum;
		}
		#pragma omp barrier
		
		#pragma omp parallel for shared(j,L,A,U,n) private(i,k,sum)
		for (i = j; i < n; i++) {
			sum = 0;
			for(k = 0; k < j; k++) {
				sum = sum + L[j][k] * U[k][i];
			}
			if (L[j][j] == 0) {				
				printf("det(L) close to 0!\n Can't divide by 0...\n");
				exit(EXIT_FAILURE);
			}
			U[j][i] = (A[j][i] - sum) / L[j][j];
		}
	}	
}

void crout_2(double **A, double **L, double **U, int n, int t) {
	// int i, j, k;
	double sum = 0;

	for (int i = 0; i < n; i++) {
		U[i][i] = 1;
	}

	for (int j = 0; j < n; j++) {

		double sum = 0;
		for(int k=0; k<j; k++)
		{
			sum = sum + L[j][k]*U[k][j];
		}
		L[j][j] = A[j][j] - sum;

		#pragma omp sections
		{
			#pragma omp section
			for (int i = j+1; i < n; i++) {
				double sum = 0;
				for (int k = 0; k < j; k++) {
					sum = sum + L[i][k] * U[k][j];	
				}
				L[i][j] = A[i][j] - sum;
			}

			#pragma omp section
			for (int i = j; i < n; i++) {
				double sum = 0;
				for(int k = 0; k < j; k++) {
					sum = sum + L[j][k] * U[k][i];
				}
				if (L[j][j] == 0) {
					printf("det(L) close to 0!\n Can't divide by 0...\n");
					exit(EXIT_FAILURE);
				}
				U[j][i] = (A[j][i] - sum) / L[j][j];
			}
		}
	}
}

void crout_3(double **A, double **L, double **U, int n, int t) {
	double sum = 0;

	for (int i = 0; i < n; i++) {
		U[i][i] = 1;
	}

	for (int j = 0; j < n; j++) {

		double sum = 0;
		for(int k=0; k<j; k++)
		{
			sum = sum + L[j][k]*U[k][j];
		}
		L[j][j] = A[j][j] - sum;

		#pragma omp sections
		{
			#pragma omp section
			#pragma omp parallel for
			for (int i = j+1; i < n; i++) {
				double sum = 0;
				for (int k = 0; k < j; k++) {
					sum = sum + L[i][k] * U[k][j];	
				}
				L[i][j] = A[i][j] - sum;
			}

			#pragma omp section
			#pragma omp parallel for
			for (int i = j; i < n; i++) {
				double sum = 0;
				for(int k = 0; k < j; k++) {
					sum = sum + L[j][k] * U[k][i];
				}
				if (L[j][j] == 0) {
					printf("det(L) close to 0!\n Can't divide by 0...\n");
					exit(EXIT_FAILURE);
				}
				U[j][i] = (A[j][i] - sum) / L[j][j];
			}
		}
	}
}



int main(int argc, char* argv[]){
    
    
    int n=atoi(argv[1]),t=atoi(argv[3]),s=atoi(argv[4]);
    
    double** mat=malloc(n*sizeof(double*));
	f(i,n)
		mat[i]=malloc(n*sizeof(double));
	
	FILE *file;
	file=fopen(argv[2], "r");
  	f(i,n){
  		f(j,n){
  			if (!fscanf(file, "%lf", &mat[i][j])) 
           			break;
     // 			printf("%lf\n",mat[i][j]); //Use lf format specifier, \n is for new line
  		}
  	}
	fclose(file);
	
	double** L = malloc(n*sizeof(double*));
	f(i,n)
		L[i]=malloc(n*sizeof(double));
	
	double** U=malloc(n*sizeof(double*));
	f(i,n)
		U[i]=malloc(n*sizeof(double));

    	double start, end;
  
    	start = omp_get_wtime();
	
	printf("s %d\n",s);
	
  	if (s==0){
		seq_crout(mat, L, U, n);
	}
	else if(s==1){
		crout_1(mat,L,U,n,t);
	}
	else if(s==2){
		crout_2(mat,L,U,n,t);
	}
	else if(s==3){
		crout_3(mat,L,U,n,t);
	}

	end = omp_get_wtime();
 
    	// Calculating total time taken by the program.
    double time_taken = (double)(end - start);
    printf("Time taken by program is : %lf sec\n",time_taken);
    
	
	char str[20]="";
	strcat(str,"output_L_");
	strcat(str,argv[4]);
	strcat(str,"_");
	strcat(str,argv[3]);
	strcat(str,".txt");
	FILE *file1 = fopen(str,"w");
	
    str[7] = 'U';
	FILE *file2 = fopen(str,"w");
   
   	f(i,n){
  		f(j,n){
  			fprintf(file1,"%lf ",L[i][j]);
  		}
  		if(i<n-1)
  			fprintf(file1,"\n");
  	}
  	
  	fclose(file1);
  	printf("L written successfully\n");
  	
  	f(i,n){
  		f(j,n){
  			fprintf(file2,"%lf ",U[i][j]);
  		}
  		if(i<n-1)
  			fprintf(file2,"\n");
  	}
  	
  	fclose(file2);
  	printf("U written successfully\n");
  	
  	free(mat);
  	free(L);
  	free(U);  	
   
    return 0;
}
		
