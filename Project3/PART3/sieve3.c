/*
 *   Sieve of Eratosthenes
 *
 *   Programmed by Michael J. Quinn
 *
 *   Last modification: 7 September 2001
 */

#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN(a,b)  ((a)<(b)?(a):(b))

int main (int argc, char *argv[])
{
   long long   count;        /* Local prime count */
   double elapsed_time; /* Parallel execution time */
   long long    first;        /* Index of first multiple */
   long long    global_count; /* Global prime count */
   long long    high_value;   /* Highest value on this proc */
   long long    i;
   int    id;           /* Process ID number */
   long long    index = 0;        /* Index of current prime */
   long long    low_value;    /* Lowest value on this proc */
   char  *marked;       /* Portion of 2,...,'n' */
   long long    n;            /* Sieving from 2, ..., 'n' */
   int    p;            /* Number of processes */
   long long    proc0_size;   /* Size of proc 0's subarray */
   long long    prime;        /* Current prime */
   long long    size;         /* Elements in 'marked' */
   
   long long    primeArray_size;    /* local size  */
   char *primeArray_marked;          /* Store Local Primes */
   long long    primeArray_firstElement;     /* local First */

   long long cache_size; /* Cache Size */
   long long last;
    

   MPI_Init (&argc, &argv);

   /* Start the timer */

      

   MPI_Comm_rank (MPI_COMM_WORLD, &id);
   MPI_Comm_size (MPI_COMM_WORLD, &p);
   MPI_Barrier(MPI_COMM_WORLD);
   elapsed_time = -MPI_Wtime();

   if (argc != 2) {
      if (!id) printf ("Command line: %s <m>\n", argv[0]);
      MPI_Finalize();
      exit (1);
   }

   n = atoll(argv[1]);


   /* Figure out this process's share of the array, as
      well as the integers represented by the first and
      last array elements */

   /*
    *  3+ to remove 1,2,3
    */

   low_value = 3 +  2 * (id*((n/2)-1)/p);
   high_value = 1 + 2 * ((id+1)*((n/2)-1)/p);
   size = ((high_value - low_value)/2) + 1;

   primeArray_size = (long long) sqrt((double) n);


   /* Bail out if all the primes used for sieving are
      not all held by process 0 */

   proc0_size = ((n/2)-1)/p;

   if (( (2 * proc0_size) + 3) < (long long) sqrt((double) n)/2) {
      if (!id) printf ("Too many processes\n");
      MPI_Finalize();
      exit (1);
   }

   /* Allocate this process's share of the array. */

   marked = (char *) malloc (size);
	
   primeArray_marked = (char *) malloc (primeArray_size);

   if (marked == NULL || primeArray_marked == NULL) {
      printf ("Cannot allocate enough memory\n");
      MPI_Finalize();
      exit (1);
   }

   for (i = 0; i < size; i++) marked[i] = 0;

   for (i = 0; i < primeArray_size; i++) primeArray_marked[i] = 0;   


   cache_size = 1<<19; // pow (2,19) = 2 ^ 9 k = 512 k
   last = 0;

   do {

   //if (!id) index = 0;
   prime = 3;
   index = 0;

   	do {
      		if (prime * prime > low_value)
         		first = (prime * prime - low_value)/2;
      		else {
         		if (!(low_value % prime)) first = 0;
         		else  {
            			first = prime - (low_value % prime);
            			if ( (low_value+first)%2 == 0  ) {
               				first = first + prime;
            		        }
                                first = first / 2 ;
                        }
                }


      		primeArray_firstElement = (prime * prime -3)/2; 
      

      		for (i = first + last ; i < MIN(last + cache_size,size); i += prime) marked[i] = 1;

      		for (i = primeArray_firstElement; i<primeArray_size; i += prime) primeArray_marked[i] = 1;  


      		while (primeArray_marked[++index]);
         		prime = 2 * index + 3;

   	} while (prime * prime <= high_value);
	
	low_value = low_value + cache_size * 2;
        last = last + cache_size; 
   } while (last<size);

   count = 0;
   global_count = 0;

   for (i = 0; i < size; i++)
      if (!marked[i]) count++;

   MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM,
                          0, MPI_COMM_WORLD);

   /* Stop the timer */

   elapsed_time += MPI_Wtime();


   /* Print the results */

   if (!id) {
      global_count++; // 2 is even number and also a prime, so lets add it :D
      printf ("There are %d primes less than or equal to %lld\n",
              global_count, n);
      printf ("SIEVE (%d) %10.6f\n", p, elapsed_time);
   }
   MPI_Finalize ();
   return 0;
}
