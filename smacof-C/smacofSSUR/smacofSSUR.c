#include "../include/smacof.h"

int n = 0;

void outputCSV( double *embedding ) {

	// open the file
	FILE *fp = NULL;
	if( (fp = fopen( "smacof_result.csv", "w" )) == NULL ) {

		printf("ERROR: Can't open points output file %s\n", "smacof_result.csv" );
		exit( 0 );
	}

	// output header
	fprintf(fp,"X,Y\nDOUBLE,DOUBLE\n"); 

	// output data to file
	for( int i = 0; i < n; i++ ) {
		for( int j = 0; j < 2; j++ ) {
			fprintf(fp, "%f",  embedding[2*i+j] );
			if( j < 1) 
				fprintf(fp, "," );
		}
		fprintf(fp,"\n");
	}

	// close the file
	fclose( fp );
}

int main(void) {

	
	
	char line[65536];	// line of input buffer
	char item[512];		// single number string
	
	FILE *fp = fopen( "distmat_inv.csv", "r" );
	if( fp == NULL ) {
		printf( "ERROR cannot open %s\n", "distmat.csv" );
		exit( 0 );
	}

	// get dataset statistics
	int line_num = 0;
	
	while( fgets( line, 65536, fp ) != NULL ) {

		// count the number of points (for every line)
		line_num++;
	}

	fclose( fp );
	float n_float = 0.5 + sqrt(0.25+2*line_num);
	n = (int) n_float;
	
	double* distmat = new double[line_num];
	
	fp = fopen( "distmat_inv.csv", "r" );
	if( fp == NULL ) {
		printf( "ERROR cannot open %s\n", "distmat.csv" );
		exit( 0 );
	}
	
	int k = 0;
	
	while( fgets( line, 65536, fp ) != NULL ) {

		int done = 0;
		int i = 0;
		int j = 0;
		while( !done ) {

				// parse character data

  				if( line[i] == '\n' || line[i] == '\0' ) {

					item[j] = '\0';
					distmat[k++] = (double) atof( item );
					done++;
				}
				else if( line[i] != ' ' ) {

					item[j++] = line[i];
				}
				i++;
		}
		
	}
  
  //distmat = delta
  double* x = new double[2*n];
  for (int g = 0; g < 2*n; g++) {
    x[g] = 0.1*g;
  }  
  int p = 2, itmax = 10000;
  bool verbose = true, speedup = true;
  double eps = 1e-14;
  (void)smacofSSUR(distmat, x, &n, &p, &speedup, &itmax, &eps, &verbose);
  printf("\n\n");
  //(void)smacofPrintMatrix(4, 2, 15, 10, x);
  //double *dist = (double *)calloc((size_t) (n * (n - 1) / 2), sizeof(double));
  //(void)smacofDistances(x, n, p, dist);
  //(void)smacofGradientU(dist, distmat, n, p, x);
  //(void)smacofPrintMatrix(4, 2, 15, 10, x);
  outputCSV(x);
  free(distmat);
  free(x);
  return EXIT_SUCCESS;
}


void smacofSSUR(double *delta, double *x, const int *nobjects,
                const int *ndimensions, const bool *speedup, const int *itmax,
                const double *eps, const bool *verbose) {
  int n = *nobjects, p = *ndimensions, acc = *speedup;
  int m = n * (n - 1) / 2, np = n * p, itel = 1;
  double sold = 0.0, snew = 0.0;
  // normalize delta
  (void)smacofNormDeltaU(delta, m);
  // compute initial distances
  double *dist = (double *)calloc((size_t)m, sizeof(double));
  (void)smacofDistances(x, n, p, dist);
  // scale initial configuration
  (void)smacofScaleXU(x, dist, delta, m, np);
  // compute initial stress
  sold = smacofLossU(dist, delta, m);
  while (true > false) {
    (void)smacofGuttmanU(dist, delta, n, p, acc, x);
    (void)smacofDistances(x, n, p, dist);
    snew = smacofLossU(dist, delta, m);
    //if (verbose) {
      //printf("itel = %4d sold = %15.10f snew = %15.10f\n", itel, sold, snew);
    //}
    if ((itel == *itmax  || ((sold - snew) < *eps))) {
      break;
    }
    itel += 1;
    sold = snew;
  }
  free(dist);
  return;
}

void smacofSSUR_dist(double *delta, double *x, const int *nobjects,
                const int *ndimensions, const bool *speedup, const int *itmax,
                const double *eps, const bool *verbose) {
  int n = *nobjects, p = *ndimensions, acc = *speedup;
  int m = n * (n - 1) / 2, np = n * p, itel = 1;
  double sold = 0.0, snew = 0.0;
  // normalize delta
  (void)smacofNormDeltaU(delta, m);
  // compute initial distances
  double *dist = (double *)calloc((size_t)m, sizeof(double));
  (void)smacofDistances(x, n, p, dist);
  // scale initial configuration
  (void)smacofScaleXU(x, dist, delta, m, np);
  // compute initial stress
  sold = smacofLossU(dist, delta, m);
  while (true > false) {
    (void)smacofGuttmanU(dist, delta, n, p, acc, x);
    (void)smacofDistances(x, n, p, dist);
    snew = smacofLossU(dist, delta, m);
    if (verbose) {
      printf("itel = %4d sold = %15.10f snew = %15.10f\n", itel, sold, snew);
    }
    if ((itel == *itmax) || ((sold - snew) < *eps)) {
      break;
    }
    itel += 1;
    sold = snew;
  }
  free(dist);
  return;
}
