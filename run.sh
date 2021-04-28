if [ $4 == "4" ] ; 
then mpiexec -n $3 mpi $1 $2 $3 $4; 
else ./seq $1 $2 $3 $4;
fi

