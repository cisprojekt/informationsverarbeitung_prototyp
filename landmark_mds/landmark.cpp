// glimmer.cpp : Console program to compute Glimmer CPU MDS on a set of input coordinates
//				
//				Stephen Ingram (sfingram@cs.ubc.ca) 02/08
//

#include <Eigen/Dense>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <string.h>


/*
	32 bit random number generation (default is 16 bit)
*/
int myrand( ) {

	unsigned int n = (unsigned int)rand();
	unsigned int m = (unsigned int)rand();

	return ((int)((n << 16) + m));
}

float max( float a, float b) {

	return (a < b)?b:a;
}
float min( float a, float b) {

	return (a < b)?a:b;
}

/*
	Output the embedding coordinates to a CSV file
*/
void outputCSV( const char *filename, float *embedding ) {

	// open the file
	FILE *fp = NULL;
	if( (fp = fopen( filename, "w" )) == NULL ) {

		printf("ERROR: Can't open points output file %s\n", filename );
		exit( 0 );
	}

	// output header
	fprintf(fp,"X,Y\nDOUBLE,DOUBLE\n"); 

	// output data to file
	for( int i = 0; i < N; i++ ) {
		for( int j = 0; j < n_embedding_dims; j++ ) {
			fprintf(fp, "%f",  embedding[(i*n_embedding_dims)+j] );
			if( j < n_embedding_dims-1) 
				fprintf(fp, "," );
		}
		fprintf(fp,"\n");
	}

	// close the file
	fclose( fp );
}


/*
	main function
*/
int main(int argc, char* argv[])
{

	
	char line[65536];	// line of input buffer
	char item[512];		// single number string
	
	FILE *fp = fopen( "distmat.csv", "r" );
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
	printf("linenum: %d\n",line_num);
	float N_float = 0.5 + sqrt(0.25+2*line_num);
	N = (int) N_float;
	
	float* distmat = new float[line_num];
	
	fp = fopen( "distmat.csv", "r" );
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
					distmat[k++] = (float) atof( item );
					done++;
				}
				else if( line[i] != ' ' ) {

					item[j++] = line[i];
				}
				i++;
		}
		
	}
	
	


	// begin timing -------------------------------------BEGIN TIMING
	clock_t start_time1 = clock();

	// allocate embedding and associated data structures
	coords	= (float *)malloc(sizeof(float)*2*N);

	landmark_mds(distmat, lands, dim)
	
	
	

	clock_t start_time2 = clock();

	printf("Anzahl Punkte: %d\nbenötigte Zeit %f\n", N, static_cast<float>(start_time2-start_time1)/CLOCKS_PER_SEC);
	printf("stop_iteration %d\n", stop_iteration);

	if (strcmp(argv[1],"NONE")) {
		outputCSV(argv[1],coords);
	}
	// quit
	return 0;
}

