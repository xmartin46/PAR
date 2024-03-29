#include "heat.h"
#include <omp.h>
/*
 * Function to copy one matrix into another
 */

void copy_mat (double *u, double *v, unsigned sizex, unsigned sizey)
{
	#pragma omp parallel for collapse(2)
	for (int i=1; i<=sizex-2; i++)
		for (int j=1; j<=sizey-2; j++) 
			v[ i*sizey+j ] = u[ i*sizey+j ];
}

/*
 * Blocked Jacobi solver: one iteration step
 */
double relax_jacobi (double *u, double *utmp, unsigned sizex, unsigned sizey)
{
	double diff, sum=0.0;

	#pragma omp parallel reduction(+: sum) private(diff)
	{
		int howmany = omp_get_num_threads();
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
