#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <malloc.h>
#include <mpi.h>

#define f(a,b) for(int a =0; a<b; a++)

int main(int argc, char* argv[]){
    
    int n=atoi(argv[1]),t=atoi(argv[3]),s=atoi(argv[4]);
	double** mat=malloc(n*sizeof(double*));
    f(i,n)
		mat[i]=malloc(n*sizeof(double));
	
	int rank;
	int size;
	int i, j, k;
	double sum = 0;
	int i_start, i_end;
	MPI_Status status;
	MPI_Request request;
	FILE *file;
	
	int subsection=0;
	
	double** L = malloc(n*sizeof(double*));
	f(i,n)
		L[i]=malloc(n*sizeof(double));
	
	double** U=malloc(n*sizeof(double*));
	f(i,n)
		U[i]=malloc(n*sizeof(double));

	double start,end;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	
	if(rank == 0){
		start = MPI_Wtime();
		file=fopen(argv[2], "r");
  		f(i,n){
  			f(j,n){
  				if (!fscanf(file, "%lf", &mat[i][j])) 
           			break;
//           		printf("%d %d %lf\n",i,j,mat[i][j]); //Use lf format specifier, \n is for new line
			}
  		}
  		printf("A read succesfully\n");
		fclose(file);
	}

    f(g,n)
		MPI_Bcast(&mat[g][0],n,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	/*
	f(g,n){
		f(h,n){
			printf("%lf ",mat[g][h]);
		}
		printf(" rank %d\n",rank);
	}  
    */
	
	for (int i = 0; i < n; i++) {
		U[i][i] = 1;
	}


	for (j = 0; j < n; j++) {
	    subsection = (n-j)/size;
		i_start = j + rank*subsection;
		if(rank==size-1)
			i_end = n;
		else
			i_end = i_start + subsection;
			
		for(i = j; i<n; i++){
			if(size>1){
				MPI_Irecv(&L[i][j],1, MPI_DOUBLE, MPI_ANY_SOURCE,(i+1)*100, MPI_COMM_WORLD, &request);
		//		printf("process %d recieved L%d%d \n",rank,i,j);
			}
		}
			
		for (i = i_start; i < i_end; i++) {
			sum = 0;
			for (k = 0; k < j; k++) {
				sum = sum + L[i][k] * U[k][j];	
			}
			L[i][j] = mat[i][j] - sum;
		
			if(size>1){
				for(int p=0; p<size; p++){
					MPI_Send( &L[i][j], 1 , MPI_DOUBLE,p, (i+1)*100, MPI_COMM_WORLD); 
			//		printf("process %d sent L%d%d to %d %lf\n",rank,i,j,p,L[i][j]);
				}	
			}
		}
		//printf("End of send: process %d sent from %d to %d with j %d\n",rank, i_start, i_end, j);
		
		for(i=j;i<n;i++){
			//printf( "process %d recieved L%d%d %lf\n",rank,i,j,L[i][j]);
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
		
		
		for(i = j; i<n; i++){
			if(size>1){
				MPI_Irecv(&U[j][i],1, MPI_DOUBLE, MPI_ANY_SOURCE,(i+1)*100, MPI_COMM_WORLD, &request);
				//printf("process %d recieved U%d%d \n",rank,j,i);
			}
		}
		
		for (i = i_start; i < i_end; i++) {
			sum = 0;
			for(k = 0; k < j; k++) {
				sum = sum + L[j][k] * U[k][i];
			}
			/*
			if (L[j][j] == 0) {				
				exit(0);
			}
			*/
			U[j][i] = (mat[j][i] - sum) / L[j][j];
			if(size>1){
				for(int p=0; p<size; p++){
					MPI_Send( &U[j][i], 1 , MPI_DOUBLE,p,(i+1)*100, MPI_COMM_WORLD); 
					//printf("process %d sent U%d%d to %d %lf\n",rank,j,i,p,U[j][i]);
				}	
			}
		}
		
		//printf("Send over : process %d j %d\n",rank,j);
		
		for(i=j;i<n;i++)
			//printf("process %d recieved U%d%d %lf\n",rank,j,i,U[j][i]);

		
		
		
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
	// Writing in files
	if(rank ==0){
	
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
	  	printf("L written succesfully\n");
	  	
	  	
	  	f(i,n){
	  		f(j,n){
	  			fprintf(file2,"%lf ",U[i][j]);
	  		}
	  		if(i<n-1)
	  			fprintf(file2,"\n");
	  	}
	  	
	  	fclose(file2);
	  	printf("U written succesfully\n");
	  	
	  	end = MPI_Wtime();
	  	
	  	printf("time taken : %lf\n",end-start);
	}
	free(mat);
  	free(L);
  	free(U);  
  	MPI_Finalize();	
   
     return 0;
}
		
