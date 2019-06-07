double relax_gauss (double *u, unsigned sizex, unsigned sizey)
{
    double unew, diff, sum=0.0;

    
    #pragma omp parallel private(unew, diff) reduction(+: sum)
    {
		int howmany = omp_get_num_threads();
		
		#pragma omp for ordered(2) 
		for (int blockid_i = 0; blockid_i < howmany; ++blockid_i) {
			for (int blockid_j = 0; blockid_j < howmany; ++blockid_j) {
				
				int i_start = lowerb(blockid_i, howmany, sizex);
				int i_end = upperb(blockid_i, howmany, sizex);
				int j_start = lowerb(blockid_j, howmany, sizey);
				int j_end = upperb(blockid_j, howmany, sizey);
				
				
				#pragma omp ordered depend(sink: blockid_i-1,blockid_j)
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
	}

    return sum;
}
