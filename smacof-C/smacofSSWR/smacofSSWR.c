#include "../include/smacof.h"
#include <time.h>

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
  double* x_inter = new double[2*n];
  for (int g = 0; g < 2*n; g++) {
    x_inter[g] = 0.1*g;
  }
  double* x = new double[2*n];
  for (int g = 0; g < 2*n; g++) {
    x[g] = 0.1*g;
  }  
  double* w = new double[n*(n-1)/2];
  for (int g = 0; g < (n*(n-1)/2); g++) {
    w[g] = 1.0;
  }
  double min_err = 10000000;
  int num_tries = 10;
  int p = 2, itmax = 10000;
  bool verbose = false, speedup = false;
  double eps = 1e-15;
  for (int i = 0; i <= num_tries; i++) {
    srand48((long int)time(NULL)); 
    for (int g = 0; g < 2*n; g++) {        
      x_inter[g] = drand48();
    }
    double recent_err = smacofSSWR(distmat, w, x_inter, &n, &p, &speedup, &itmax, &eps, &verbose);
    if (recent_err < min_err) {
      min_err = recent_err;
      for (int g = 0; g < 2*n; g++) {
        x[g] = x_inter[g];
      }
    }
  }
  //(void)smacofSSWI(distmat, w, x, &n, &p, &speedup, &itmax, &eps, &verbose);
  printf("\n\n");
  //(void)smacofPrintMatrix(4, 2, 15, 10, x);
  //double *dist = (double *)calloc((size_t) (n * (n - 1) / 2), sizeof(double));
  //(void)smacofDistances(x, n, p, dist);
  //(void)smacofGradientU(dist, distmat, n, p, x);
  //(void)smacofPrintMatrix(4, 2, 15, 10, x);
  outputCSV(x);
  free(distmat);
  free(x);
  free(w);
  free(x_inter);
  return EXIT_SUCCESS;
}

// Single Matrix, Symmetric, Weighted Stress, Ratio Transform

double smacofSSWR(double *delta, double *w, double *x, const int *nobjects,
                const int *ndimensions, const bool *speedup, const int *itmax,
                const double *eps, const bool *verbose) {
  int n = *nobjects, p = *ndimensions, acc = *speedup;
  int m = n * (n - 1) / 2, np = n * p, itel = 1;
  double sold = 0.0, snew = 0.0;
  // normalize w
  (void)smacofNormW(w, m);
  // normalize delta
  (void)smacofNormDeltaW(w, delta, m);
  // Compute the MP inverse of V and test for irreducibility
  double *vinv = (double *)calloc((size_t)m, sizeof(double));
  (void)smacofInvertVW(w, vinv, n, m);
  // Compute initial distances
  double *dist = (double *)calloc((size_t)m, sizeof(double));
  (void)smacofDistances(x, n, p, dist);
  // Scale initial configuration and distances
  (void)smacofScaleXW(x, dist, w, delta, m, np);
  // Compute initial stress
  sold = smacofLossW(dist, w, delta, m);
  while (true) {
    (void)smacofGuttmanW(dist, w, delta, vinv, n, p, acc, x);
    (void)smacofDistances(x, n, p, dist);
    snew = smacofLossW(dist, w, delta, m);
    //if (verbose) {
      //printf("itel = %4d sold = %15.10f snew = %15.10f\n", itel, sold, snew);
    //}
    if ((itel == *itmax) || ((sold - snew) < *eps)) {
      break;
    }
    itel += 1;
    sold = snew;
  }
  free(vinv);
  free(dist);
  return snew;
}

void smacofSSWRdist(double *delta, double *w, double *x, const int *nobjects,
                const int *ndimensions, const bool *speedup, const int *itmax,
                const double *eps, const bool *verbose) {
  int n = *nobjects, p = *ndimensions, acc = *speedup;
  int m = n * (n - 1) / 2, np = n * p, itel = 1;
  double sold = 0.0, snew = 0.0;
  double *dist = (double *)calloc((size_t)m, sizeof(double));
  double *vinv = (double *)calloc((size_t)m, sizeof(double));
  (void)smacofScaleXW(x, dist, w, delta, m, np);
  // Compute initial stress
  sold = smacofLossW(dist, w, delta, m);
  while (true) {
    (void)smacofGuttmanW(dist, w, delta, vinv, n, p, acc, x);
    (void)smacofDistances(x, n, p, dist);
    snew = smacofLossW(dist, w, delta, m);
    if (verbose) {
      printf("itel = %4d sold = %15.10f snew = %15.10f\n", itel, sold, snew);
    }
    if ((itel == *itmax) || ((sold - snew) < *eps)) {
      break;
    }
    itel += 1;
    sold = snew;
  }
  free(vinv);
  return;
}

