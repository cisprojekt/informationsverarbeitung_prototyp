#include <cmath>
#include <string>
#include <emscripten.h>

extern "C" {
    EMSCRIPTEN_KEEPALIVE
    float calculateEuclideanDistance(float* vector1, float* vector2, int dimension) {

        float sumOfSquares = 0.0;
        for (size_t i = 0; i < dimension; i++) {
            float diff = vector2[i] - vector1[i];
            sumOfSquares += diff * diff;
        }

        float distance = std::sqrt(sumOfSquares);
        return distance;
    }

    EMSCRIPTEN_KEEPALIVE
    float* calculateDistanceMatrix(float* array, int num_vectors, int dimension) {

        float** distanceArray = new float*[num_vectors];

        for (size_t i = 0; i < num_vectors; ++i) {
            distanceArray[i] = new float[i + 1];
        }
        for (size_t i = 0; i < num_vectors; i++) {
            for (size_t j = 0; j < i+1; j++) {
                float distance = calculateEuclideanDistance(array+i*dimension, array+j*dimension, dimension);
                distanceArray[i][j] = distance;
            }
        }
        //flatten the array
        float* flatArray = new float[num_vectors * (num_vectors + 1) / 2];
        int index = 0;
        for (size_t i = 0; i < num_vectors; i++) {
            for (size_t j = 0; j < i+1; j++) {
                flatArray[index] = distanceArray[i][j];
                index++;
            }
        }
        delete[] distanceArray;

        return flatArray;
    }




    EMSCRIPTEN_KEEPALIVE
    int calculateHammingDistance(const std::string& str1, const std::string& str2) {
        if (str1.length() != str2.length()) {
            // Handle error: strings must have the same length
            return -1;
        }

        int distance = 0;
        for (size_t i = 0; i < str1.length(); i++) {
            if (str1[i] != str2[i]) {
                distance++;
            }
        }

        return distance;
    }
}




