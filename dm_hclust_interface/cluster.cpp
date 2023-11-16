#include <math.h>
#include <string.h>
#include <string>
#include <math.h>
#include <vector>
#include <stdio.h>
#include <algorithm>

// include the interface for
// Daniel Müllner's fastclust algorithm
#include "fastcluster.h"

// stores points as x and y coordinates
class Point {
    public:
        double x;
        double y;

        Point(const double xCord, const double yCord) {
            x = xCord;
            y = yCord;
        }

        Point(const Point& p) {
            x = p.x;
            y = p.y;
        }
};

// calculate euclidian distance between two points
double distance(const Point& p1, const Point& p2) {
    double dist =  sqrt(pow(p1.x-p2.x, 2)+pow(p1.y-p2.y, 2));
    return dist;
}

int main(int argc, char** argv)
{

    int i,j,k;

    // parse command line input
    std::string opt_infile;
    int optMethod = HCLUST_METHOD_SINGLE; // default method
    const char* usemsg = "Usage: hclust-demo <file> [-m (single|complete|average|median)]\n";
    for (i=1; i<argc; i++) {

        // check if method is given in command line
        if (0 == strcmp(argv[i], "-m")) {
            i++;
            if (i<argc) {
                if (0 == strcmp(argv[i], "single")) optMethod = HCLUST_METHOD_SINGLE;
                else if (0 == strcmp(argv[i], "complete")) optMethod = HCLUST_METHOD_COMPLETE;
                else if (0 == strcmp(argv[i], "average")) optMethod = HCLUST_METHOD_AVERAGE;
                else if (0 == strcmp(argv[i], "median")) optMethod = HCLUST_METHOD_MEDIAN;
                else {
                    fputs(usemsg, stderr);
                    return 1;
                }
            } else {
                fputs(usemsg, stderr);
                return 1;
            }
        } else {
            opt_infile = argv[i];
        }
    }
    if (opt_infile == "") {
        fputs(usemsg, stderr);
        return 1;
    }

    // read points from input file
    double x, y;
    std::vector<Point> points;
    FILE* f = fopen(opt_infile.c_str(), "r");
    if (!f) {
        fprintf(stderr, "Cannot open '%s'\n", opt_infile.c_str());
        return 2;
    }
    int nPoints = 0;

    // read until end of file stream has been reached
    while (!feof(f)) {
        nPoints++;
        k = fscanf(f, "%lf, %lf", &x, &y);
        if (k != 2) {
            fprintf(stderr, "Input does not match x and y coordinates");
            return 3;
        }
        points.push_back(Point(x, y));
    }
    fclose(f);

    // calculate condensed distance matrix
    double* distMat = new double[(nPoints*(nPoints-1))/2];
    k = 0;
    for (i=0; i<nPoints; i++) {
        for (j=i+1; j<nPoints; j++) {
            distMat[k] = distance(points[i], points[j]);
            k++;
        }
    } 

    // clustering call
    int *merge = new int[2*(nPoints-1)];
    double* height = new double[nPoints-1];
    hclust_fast(nPoints, distMat, optMethod, merge, height);

    int* labels = new int[nPoints];
    cutree_k(nPoints, merge, 5, labels);

    // print results
    for (i=0; i<nPoints; i++) {
        printf("Point %i (%f, %f) is in cluster %i\n", i, points[i].x, points[i].y, labels[i]);
    }

    // clean up
    delete[] distMat;
    delete[] merge;
    delete[] height;
    delete[] labels;

    return 0;
}