#include "../include/smacof.h"

int n = 0;


void smacofSSUO(double *delta, double *x, const int *nobjects,
                const int *ndimensions, const bool *speedup, const int *itmax,
                const double *eps, const bool *verbose) {
  int n = *nobjects, p = *ndimensions, acc = *speedup;
  int m = n * (n - 1) / 2, np = n * p, itel = 1;
  double sold = 0.0, smid = 0.0, snew = 0.0;
  // normalize delta
  (void)smacofNormDeltaU(delta, m);
  // compute initial distances
  double *dist = (double *)calloc((size_t)m, sizeof(double));
  (void)smacofDistances(x, n, p, dist);
  // scale initial configuration
  (void)smacofScaleXU(x, dist, delta, m, np);
  // compute initial stress
  sold = smacofLossU(dist, delta, m);
  // make the ordinal structure
  //(void) smacofBlockSort(delta, NULL, m, 1,NULL);
  while (true > false) {
    (void)smacofGuttmanU(dist, delta, n, p, acc, x);
    (void)smacofDistances(x, n, p, dist);
    smid = smacofLossU(dist, delta, m);
    (void)smacofMissingU(dist, delta, m);
    snew = smacofLossU(dist, delta, m);
    if (verbose) {
      printf("itel = %4d sold = %15.10f smid = %15.10f snew = %15.10f\n", itel,
             sold, smid, snew);
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

void smacofNormDeltaU(double *delta, const int m) {
  double s = 0.0, r = 0.0;
  for (int k = 0; k < m; k++) {
    s += SQUARE(delta[k]);
  }
  r = sqrt(((double)m) / s);
  for (int k = 0; k < m; k++) {
    delta[k] *= r;
  }
  return;
}

void smacofScaleXU(double *x, double *dist, const double *delta, const int m,
                   const int np) {
  double sd1 = 0.0, sd2 = 0.0;
  for (int k = 0; k < m; k++) {
    sd1 += dist[k] * delta[k];
    sd2 += SQUARE(dist[k]);
  }
  double lbd = sd1 / sd2;
  for (int k = 0; k < m; k++) {
    dist[k] *= lbd;
  }
  for (int k = 0; k < np; k++) {
    x[k] *= lbd;
  }
}

double smacofLossU(const double *dist, const double *delta, const int m) {
  double stress = 0.0;
  for (int k = 0; k < m; k++) {
    stress += SQUARE(delta[k] - dist[k]);
  }
  return stress;
}

void smacofGuttmanU(const double *dist, const double *delta, const int n,
                    const int p, const bool speedup, double *x) {
  int k;
  for (int s = 1; s <= p; s++) {
    double *y = (double*) calloc((size_t)n, sizeof(double));
    for (int i = 1; i <= n; i++) {
      for (int j = 1; j <= n; j++) {
        if (j == i) {
          continue;
        }
        if (i > j) {
          k = SINDEX(i, j, n);
        }
        if (j > i) {
          k = SINDEX(j, i, n);
        }
        if (dist[k] < 1e-15) {
          continue;
        }
        y[VINDEX(i)] +=
            delta[k] * (x[MINDEX(i, s, n)] - x[MINDEX(j, s, n)]) / dist[k];
      }
      if (speedup) {
        y[VINDEX(i)] = 2 * y[VINDEX(i)] / ((double)n) - x[MINDEX(i, s, n)];
      } else {
        y[VINDEX(i)] /= ((double)n);
      }
    }
    for (int i = 1; i <= n; i++) {
      x[MINDEX(i, s, n)] = y[VINDEX(i)];
    }
    free(y);
  }
}

void smacofGradientU(double *dist, double *delta, const int n, const int p,
                     double *x) {
  for (int s = 1; s <= p; s++) {
    double *y = (double *)calloc((size_t)n, sizeof(double));
    for (int i = 1; i <= n; i++) {
      for (int j = 1; j <= n; j++) {
        int k = 0;
        if (j == i) {
          continue;
        }
        if (i > j) {
          k = SINDEX(i, j, n);
        }
        if (j > i) {
          k = SINDEX(j, i, n);
        }
        if (dist[k] < 1e-15) {
          continue;
        }
        y[VINDEX(i)] += (1 - delta[k] / dist[k]) *
                        (x[MINDEX(i, s, n)] - x[MINDEX(j, s, n)]);
      }
    }
    for (int i = 1; i <= n; i++) {
      x[MINDEX(i, s, n)] = y[VINDEX(i)];
    }
    free(y);
  }
  return;
}

void smacofMissingU(const double *dist, double *delta, int m) {
  double acum = 0.0, bcum = 0.0, ccum = 0.0, dcum = 0.0;
  double s = 0.0, deltamin = INFINITY;
  for (int k = 0; k < m; k++) {
    acum += delta[k];
    bcum += dist[k];
    deltamin = MIN(deltamin, delta[k]);
  }
  acum /= (double)m;
  bcum /= (double)m;
  for (int k = 0; k < m; k++) {
    ccum += (delta[k] - acum) * (dist[k] - bcum);
    dcum += SQUARE(delta[k] - acum);
  }
  double a0 = ccum / dcum;
  double b0 = -(a0 * acum - bcum);
  if ((a0 > 1e-15) && ((a0 * deltamin + b0) > 1e-15)) {
    s = 0.0;
    for (int k = 0; k < m; k++) {
      delta[k] = a0 * delta[k] + b0;
      s += SQUARE(delta[k]);
    }
  } else {
    ccum = 0.0;
    dcum = 0.0;
    for (int k = 0; k < m; k++) {
      ccum += (delta[k] - deltamin) * dist[k];
      dcum += SQUARE(delta[k] - deltamin);
    }
    double a1 = ccum / dcum;
    double b1 = -a1 * deltamin;
    if (a1 > 1e-15) {
      s = 0.0;
      for (int k = 0; k < m; k++) {
        delta[k] = a0 * delta[k] + b0;
        s += SQUARE(delta[k]);
      }
    } else {
      return;
    }
  }
  double r = sqrt(((double)m) / s);
  for (int k = 0; k < m; k++) {
    delta[k] *= r;
  }
  return;
}

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
  (void)smacofSSUO(distmat, x, &n, &p, &speedup, &itmax, &eps, &verbose);
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
