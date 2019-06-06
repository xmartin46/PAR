#include "heat.h"
#include <omp.h>
/*
 * Function to copy one matrix into another
 */

void copy_mat (double *u, double *v, unsigned sizex, unsigned sizey)
{
	//#pragma omp parallel for collapse(2)
	for (int i=1; i<=sizex-2; i++)
		for (int j=1; j<=sizey-2; j++) 
			v[ i*sizey+j ] = u[ i*sizey+j ];
}

/*
 * Blocked Jacobi solver: one iteration step
 */
double relax_jacobi (double *u, double *utmp, unsigned sizex, unsigned sizey)
{
    /*double diff, sum=0.0;
  
    int howmany=4;
    for (int blockid = 0; blockid < howmany; ++blockid) {
      int i_start = lowerb(blockid, howmany, sizex);
      int i_end = upperb(blockid, howmany, sizex);
      for (int i=max(1, i_start); i<= min(sizex-2, i_end); i++) {
        for (int j=1; j<= sizey-2; j++) {
	     utmp[i*sizey+j]= 0.25 * ( u[ i*sizey     + (j-1) ]+  // left
	                               u[ i*sizey     + (j+1) ]+  // right
				       u[ (i-1)*sizey + j     ]+  // top
				       u[ (i+1)*sizey + j     ]); // bottom
	     diff = utmp[i*sizey+j] - u[i*sizey + j];
	     sum += diff * diff; 
	 }
      }
    }

    return sum;*/
    
    
	double diff, sum=0.0;

	int howmany = omp_get_num_threads();
	#pragma omp parallel reduction(+: sum) private(diff)
	{
		int blockid = omp_get_thread_num();
		int i_start = lowerb(blockid, howmany, sizex);
		int i_end = upperb(blockid, howmany, sizex);
		for (int i=max(1, i_start); i<= min(sizex-2, i_end); i++) {
			for (int j=1; j<= sizey-2; j++) {
				 utmp[i*sizey+j]= 0.25 * ( u[ i*sizey     + (j-1) ]+  // left
										   u[ i*sizey     + (j+1) ]+  // right
							   u[ (i-1)*sizey + j     ]+  // top
							   u[ (i+1)*sizey + j     ]); // bottom
				 diff = utmp[i*sizey+j] - u[i*sizey + j];
				 sum += diff * diff; 
			}
		}
	}

    return sum;
}

/*
 * Blocked Gauss-Seidel solver: one iteration step
 */
double relax_gauss (double *u, unsigned sizex, unsigned sizey)
{
    double unew, diff, sum=0.0;

    int howmany = omp_get_num_threads();
    
    #pragma omp parallel for ordered(2) private(unew,diff) reduction(+: sum)
    for (int blockid_i = 0; blockid_i < howmany; ++blockid_i) {
		for (int blockid_j = 0; blockid_j < howmany; ++blockid_j) {
			int i_start = lowerb(blockid_i, howmany, sizex);
			int i_end = upperb(blockid_i, howmany, sizex);
			int j_start = lowerb(blockid_j, howmany, sizey);
			int j_end = upperb(blockid_j, howmany, sizey);
			
			#pragma omp ordered depend(sink: blockid_i-1, blockid_j) // No cal i,j-1 perqué està al mateix thread 
			for (int i=max(1, i_start); i<= min(sizex-2, i_end); i++) {
				for (int j=max(1, j_start); j<= min(sizey-2, j_end); j++) {
					
					unew= 0.25 * ( u[ i*sizey	+ (j-1) ]+  // left
						u[ i*sizey	+ (j+1) ]+  // right
						u[ (i-1)*sizey	+ j     ]+  // top
						u[ (i+1)*sizey	+ j     ]); // bottom
					diff = unew - u[i*sizey+ j];
					sum += diff * diff; 
					u[i*sizey+j]=unew;
					
				}
			}
			#pragma omp ordered depend(source)
		}
    }

    return sum;
}
