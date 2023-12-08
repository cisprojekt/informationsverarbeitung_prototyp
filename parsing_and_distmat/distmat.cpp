#include <cmath>
#include <string>
#include <emscripten.h>

//compiled with:
//em++ -o distmat.js distmat.cpp -s "EXPORTED_FUNCTIONS=['_malloc', '_free', '_calculateEuclideanDistance', '_calculateEuclideanDistanceMatrix', '_calculateHammingDistance', '_calculateHammingDistanceMatrix']" -s "EXPORTED_RUNTIME_METHODS=['ccall']" -sALLOW_MEMORY_GROWTH -s MAXIMUM_MEMORY=1gb -s ABORTING_MALLOC=0

extern "C" {
    EMSCRIPTEN_KEEPALIVE
    float calculateEuclideanDistance(float* vector1, float* vector2, int string_length) {

        float sumOfSquares = 0.0;
        for (size_t i = 0; i < string_length; i++) {
            float diff = vector2[i] - vector1[i];
            sumOfSquares += diff * diff;
        }

        float distance = std::sqrt(sumOfSquares);
        return distance;
    }

    EMSCRIPTEN_KEEPALIVE
    float* calculateEuclideanDistanceMatrix(float* array, int num_strings, int string_length) {

        float** distanceArray = new float*[num_strings];

        for (size_t i = 0; i < num_strings; ++i) {
            distanceArray[i] = new float[i + 1];
        }
        for (size_t i = 0; i < num_strings; i++) {
            for (size_t j = 0; j < i+1; j++) {
                float distance = calculateEuclideanDistance(array+i*string_length, array+j*string_length, string_length);
                distanceArray[i][j] = distance;
            }
        }
        //flatten the array
        float* flatArray = new float[num_strings * (num_strings + 1) / 2];
        int index = 0;
        for (size_t i = 0; i < num_strings; i++) {
            for (size_t j = 0; j < i+1; j++) {
                flatArray[index] = distanceArray[i][j];
                index++;
            }
        }
        //free memory
        for (size_t i = 0; i < num_strings; ++i) {
            delete[] distanceArray[i];
        }
        delete[] distanceArray;

        return flatArray;
    }

    EMSCRIPTEN_KEEPALIVE
    int calculateHammingDistance(char* str1, char* str2, int string_length) {

        int distance = 0;
        for (size_t i = 0; i < string_length; i++) {
            if (str1[i] != str2[i]) {
                distance++;
            }
        }

        return distance;
    }


    EMSCRIPTEN_KEEPALIVE
    int* calculateHammingDistanceMatrix(char** array, int num_strings, int string_length) {

        int** distanceArray = new int*[num_strings];

        for (size_t i = 0; i < num_strings; ++i) {
            distanceArray[i] = new int[i + 1];
        }
        for (size_t i = 0; i < num_strings; i++) {
            for (size_t j = 0; j < i+1; j++) {
                int distance = calculateHammingDistance(array[i], array[j], string_length);
                distanceArray[i][j] = distance;
            }
        }
        //flatten the array
        int* flatArray = new int[num_strings * (num_strings + 1) / 2];
        int index = 0;
        for (size_t i = 0; i < num_strings; i++) {
            for (size_t j = 0; j < i+1; j++) {
                flatArray[index] = distanceArray[i][j];
                index++;
            }
        }
        //free memory
        for (size_t i = 0; i < num_strings; ++i) {
            delete[] distanceArray[i];
        }
        delete[] distanceArray;

        return flatArray;
    }
}




