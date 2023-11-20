#include <math.h>
#include <string.h>
#include <string>
#include <math.h>
#include <vector>
#include <stdio.h>
#include <algorithm>
#include <emscripten.h>
#include <msgpack.hpp>
#include <iostream>

// include the interface for
// Daniel Müllner's fastclust algorithm
#include "fastcluster.h"

// map for messagepack
struct Position {
    public:
        MSGPACK_DEFINE_MAP(positionId, x, y);
        std::string positionId;
        int x;
        int y;
};

// stores points as x and y coordinates
class Point {
    public:
        int x;
        int y;

        Point(const int xCord, const int yCord) {
            x = xCord;
            y = yCord;
        }

        Point(const Point& p) {
            x = p.x;
            y = p.y;
        }
};

// calculate euclidian distance between two points
EMSCRIPTEN_KEEPALIVE
double distance(const Point& p1, const Point& p2) {
    double dist =  sqrt(pow(p1.x-p2.x, 2)+pow(p1.y-p2.y, 2));
    return dist;
}

// cluster the points using the fastcluster repository
EMSCRIPTEN_KEEPALIVE
extern "C" char* processPoint(char* inputPointer, int inputSize, char* outputPointer, int numberPoints)
{
    // decode given data (points) using messagepack unpacker
    msgpack::object_handle objectHandle = msgpack::unpack(inputPointer, inputSize);
    msgpack::object object = objectHandle.get();

    // initialize positions vector which contains
    // the converted unpacked data
    std::vector<Position> positions;
    object.convert(positions);

    // initialize vector to store cluster IDs
    std::vector<int> positionIds;

    // initialize vector to store points
    std::vector<Point> points;

    // move points from the given vector to
    // a local points vector, also converts them
    // to type Point
    for (auto& position : positions) {
        try {
            points.push_back(Point(position.x, position.y));
        } catch (std::exception& e) {
            std::cout << e.what() << std::endl;
        }
    }

    int i,j,k;

    int nPoints = numberPoints;

    int optMethod = HCLUST_METHOD_AVERAGE; // default method

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
    std::vector<std::tuple<int, int>> pointLabels;
    for (i=0; i<nPoints; i++) {
        printf("Point %i (%f, %f) is in cluster %i\n", i, points[i].x, points[i].y, labels[i]);
        positionIds.push_back(labels[i]);
    }

    // clean up
    delete[] distMat;
    delete[] merge;
    delete[] height;
    delete[] labels;

    // pack processed data using messagepack packer
    msgpack::sbuffer sbuf;
    msgpack::pack(sbuf, positionIds);

    // get output pointer value by buffer size
    // and return buffer data
    *outputPointer = sbuf.size();
    return sbuf.data();
}
