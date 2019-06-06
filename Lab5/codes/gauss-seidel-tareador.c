#include "heat.h"
#include <stdio.h>
#include <tareador.h>

double relax_gauss (double *u, unsigned sizex, unsigned sizey)
{
    double unew, diff, sum=0.0;

    int howmany=1;
    for (int blockid = 0; blockid < howmany; ++blockid) {
      int i_start = lowerb(blockid, howmany, sizex);
      int i_end = upperb(blockid, howmany, sizex);
      for (int i=max(1, i_start); i<= min(sizex-2, i_end); i++) {
        for (int j=1; j<= sizey-2; j++) {
			tareador_start_task("gauss_seidel_innermost_task");
				
				unew= 0.25 * ( u[ i*sizey	+ (j-1) ]+  // left
					   u[ i*sizey	+ (j+1) ]+  // right
					   u[ (i-1)*sizey	+ j     ]+  // top
					   u[ (i+1)*sizey	+ j     ]); // bottom
				diff = unew - u[i*sizey+ j];
				
				sum += diff * diff; 
				
				u[i*sizey+j]=unew;
			tareador_end_task("gauss_seidel_innermost_task");
        }
      }
    }

    return sum;
}
