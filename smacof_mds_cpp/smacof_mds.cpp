#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <Eigen/Dense>
#include <iostream>

//full distance-matrix needed!

int n = 0;

double smacof_single(double* dissimilarities, double* x, double* x_inter, int n_samples, bool init=false, bool metric=true, int n_components=2, int max_iter=1000, bool verbose=0, double eps=1e-5, int random_state=0, bool normalized_stress = false) {
   

  int m = n_samples*(n_samples-1)/2;
  srand(time(NULL));
  if (!init) {
    for (int it = 0; it < 2*n_samples; it++) {
      x_inter[it] = drand48();
      //printf("x_inter[%d] = %f",it,x_inter[it]);
    }
  }
  else {
    x_inter[0] = 0.0;
    x_inter[1] = 0.0;
    x_inter[2] = 1.0;
    x_inter[3] = -1.0;
    x_inter[4] = -1.0;
    x_inter[5] = 1.0;
    x_inter[6] = 1.0;
    x_inter[7] = 1.0;
    x_inter[8] = -1.0;
    x_inter[9] = -1.0;
  }
  double old_stress = 1e15;
  double* disparities = new double[n_samples*n_samples];
  double* dis = new double[n_samples*n_samples];
  double stress = 0;
  //isotonic regression
  for (int it = 0; it < max_iter; it++) {
    for (int f = 0; f < n_samples; f++) {
      for (int g = 0; g < n_samples; g++) {
        dis[n_samples*f+g] = sqrt(pow(x_inter[2*f]-x_inter[2*g],2)+pow(x_inter[2*f+1]-x_inter[2*g+1],2));
        //printf("dis[%d] = %f\n",(n_samples*f+g), dis[n_samples*f+g]);
        
      }
    }
    if (metric) {
      for (int i = 0; i < n_samples*n_samples; i++) {
        disparities[i] = dissimilarities[i];
        //printf("disparities[%d] = %f\n",(i), disparities[i]);
      }
    }
    /* non-metric case
    else {
      disparities_flat = ir.fit_transform(sim_flat_w aka dissimilarities, dis)
      disparities = dis.copy
      disparities[sim_flat != 0] = disparities_flat
      disparities *= sqrt(n_samples * n_samples-1)/2 / (disparities**2).sum() / 2)
    }
    */
    double normalizer;
    for (int i = 0; i < n_samples*n_samples; i++) {
      //printf("difference[%d] = %f\n",(i),dis[i] - disparities[i]);
      stress += pow(dis[i] - disparities[i], 2);
    }
    stress *= 0.5;
    //printf("stress = %f\n", stress);
    if (normalized_stress) { 
      for (int i = 0; i < n_samples*n_samples; i++) {
        normalizer += pow(disparities[i],2);
      }
      normalizer *= 0.5;
      stress = sqrt(stress / normalizer);
      //printf("stress = %f\n", stress);
    }
    Eigen::MatrixXd B(n_samples,n_samples);
    Eigen::MatrixXd ratio(n_samples,n_samples);
    Eigen::MatrixXd X_m(2, n_samples);
    
    
    for (int h = 0; h < 2*n_samples; h += 2) {
      X_m(0, h/2) = x_inter[h];
      X_m(1, h/2) = x_inter[h+1];
    }
    for (int i = 0; i < n_samples*n_samples; i++) {
      if (dis[i] == 0) {dis[i] = 1e-5; 
      }
    }
    for (int i = 0; i < n_samples; i++) {
      for (int j = 0; j < n_samples; j++) {
        ratio(i, j) = disparities[i*n_samples+j]/dis[i*n_samples+j];
        B(i,j) = -ratio(i,j);
      } 
    }
    //std::cout << ratio << std::endl;
    //std::cout << B << std::endl;
    for (int i = 0; i < n_samples; i++) {
      double sum = 0;
      for (int r = 0; r < n_samples; r++) {
        sum += ratio(i, r);
      }
      B(i,i) += sum;
    }
    //std::cout << B << std::endl;
    Eigen::MatrixXd X_c(2, n_samples);
    X_c=X_m;
    X_m = X_c*B;
    X_m = 1.0/n_samples*X_m;
    //std::cout << X_m << std::endl;
    Eigen::MatrixXd X_sq(2, n_samples);
    X_sq = (X_m.array())*(X_m.array());
    //std::cout << X_sq << std::endl;
    Eigen::VectorXd interm(n_samples);
    for (int i = 0; i < n_samples; i++) {
      double sum = 0;
      for (int r = 0; r < 2; r++) {
        sum += (X_sq)(r, i);
      }
      interm(i) = sqrt(sum);
    }
    //std::cout << interm << std::endl;
    double discrepancy = interm.sum();
    //printf("discrepancy = %f\n",discrepancy);
    //printf("old_stress - stress: %f\n", old_stress - stress / discrepancy);
    if ((old_stress - stress / discrepancy) < eps)  {
      break;
    }
    old_stress = stress/discrepancy;
    //printf("old_stress %f\n", old_stress);
    for (int e = 0; e < 2*n_samples; e+=2) {
      x_inter[e] = X_m(0, e/2);
      x_inter[e+1] = X_m(1, e/2);
    }
    //printf("discrepancy %f\n", discrepancy);
    //printf("stress %f\n", stress);
  }
  /*
  for (int a = 0; a < 2*n_samples; a++) {
    printf("%f \n",x_inter[a]);
  }
  */
  delete disparities;
  delete dis;
  return stress;
}
void smacof(double* dissimilarities, double* x, double* x_inter, int n_iterations, int n_samples, bool init=false, bool metric=true, int n_components=2, int max_iter=1000, bool verbose=0, double eps=1e-5, int random_state=0, bool normalized_stress = false) {
  double min_stress = 1e18;
  double stress;
  for (int i = 0; i < n_iterations; i++) {
    stress = smacof_single(dissimilarities, x, x_inter, n_samples, init, metric, n_components, max_iter, verbose, eps, random_state, normalized_stress);
    if (stress < min_stress) {
      for (int k = 0; k < 2*n_samples; k++) {
        x[k] = x_inter[k];
      }
      min_stress = stress;
    }
  }
}

void outputCSV( double *embedding ) {

	// open the file
	FILE *fp = NULL;
	if( (fp = fopen( "scikit_result.csv", "w" )) == NULL ) {

		printf("ERROR: Can't open points output file %s\n", "scikit_result.csv" );
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
	
	FILE *fp = fopen( "distmat_full.csv", "r" );
	if( fp == NULL ) {
		printf( "ERROR cannot open %s\n", "distmat_full.csv" );
		exit( 0 );
	}

	// get dataset statistics
	int line_num = 0;
	
	while( fgets( line, 65536, fp ) != NULL ) {

		// count the number of points (for every line)
		line_num++;
	}

	fclose( fp );
	float n_float = sqrt(line_num);
	n = (int) n_float;
	
	double* distmat = new double[line_num];
	
	fp = fopen( "distmat_full.csv", "r" );
	if( fp == NULL ) {
		printf( "ERROR cannot open %s\n", "distmat_full.csv" );
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
  //for (int g = 0; g < 2*n; g++) {
    //x_inter[g] = 0.1*g;
  //}
  double* x = new double[2*n];
  //for (int g = 0; g < 2*n; g++) {
    //x[g] = 0.1*g;
  int n_iterations = 10;
  bool init = false;
  bool metric=true;
  int n_samples = n;
  int n_components=2; 
  int max_iter=300;
  bool verbose=0;
  double eps=1e-10; 
  int random_state=0;
  bool normalized_stress = false;  
  
  smacof(distmat, x, x_inter, n_iterations, n_samples, init, metric, n_components, max_iter, verbose, eps, random_state, normalized_stress);
  outputCSV(x);
  
  return 0;
}
