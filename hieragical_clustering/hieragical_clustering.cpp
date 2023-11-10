#include <iostream>
#include <cstring>
#include <memory>
#include "hieragical_clustering.h"
using namespace std;

int main() {
    cout << "Hello World";
    return 0;
    }
    //Check if node i was visited.
int is_visited(unsigned char* bitset, int i){

    return bitset[i >> 3] & (1 << (i & 7)); // i & 7 is modulo 8 works only for modulo 2^n
}

    //Mark node i as visited.
void set_visited(unsigned char* bitset, int i){

    bitset[i >> 3] |= 1 << (i & 7); // i & 7 is modulo 8 works only for modulo 2^n
}


void cluster_monocrit(double** Z, double* MC, int* T,
                            double cutoff, int n) {
    /*
    Form flat clusters by monocrit criterion.

    Parameters
    ----------
    Z : ndarray
        The linkage matrix.
    MC : ndarray
        The monotonic criterion array.
    T : ndarray
        The array to store the cluster numbers. The i'th observation belongs to
        cluster `T[i]`.
    cutoff : double
        Clusters are formed when the MC values are less than or equal to
        `cutoff`.
    n : int
        The number of observations.
    */

    int k, i_lc, i_rc, root, n_cluster = 0, cluster_leader = -1;
    unique_ptr<int[]> curr_node(new int[n]);
    int visited_size = (((n * 2) - 1) >> 3) + 1;
    unique_ptr<unsigned char[]> visited(new unsigned char[visited_size]);

    //here was some memory error catdch code

    k = 0;
    curr_node[0] = 2 * n - 2;
    while (k >= 0){
        root = curr_node[k] - n;
        i_lc = (int) Z[root][0];
        i_rc = (int) Z[root][1];

        if (cluster_leader == -1 && MC[root] <= cutoff) {  // found a cluster
            cluster_leader = root;
            n_cluster += 1;
        }
        if (i_lc >= n && !is_visited(visited.get(), i_lc)) {
            set_visited(visited.get(), i_lc);
            k += 1;
            curr_node[k] = i_lc;
            continue;
            }

        if (i_rc >= n && !is_visited(visited.get(), i_rc)) {
            set_visited(visited.get(), i_rc);
            k += 1;
            curr_node[k] = i_rc;
            continue;
        }
        if (i_lc < n) {
            if (cluster_leader == -1) n_cluster += 1;// singleton cluster
            T[i_lc] = n_cluster;
        }
        if (i_rc < n){
            if (cluster_leader == -1) n_cluster += 1;// singleton cluster
            T[i_rc] = n_cluster;
        }
        if (cluster_leader == root) { // back to the leader
            cluster_leader = -1;
        }
        k -= 1;
    }
}

void cluster_maxclust_monocrit(double** Z, double* MC, int* T,
                                    int n, int max_nc) {
    /*                                  
    Form flat clusters by maxclust_monocrit criterion.

    Parameters
    ----------
    Z : ndarray
        The linkage matrix.
    MC : ndarray
        The monotonic criterion array.
    T : ndarray
        The array to store the cluster numbers. The i'th observation belongs to
        cluster `T[i]`.
    n : int
        The number of observations.
    max_nc : int
        The maximum number of clusters.
    */ 
    int i, k, i_lc, i_rc, root, nc, lower_idx, upper_idx;
    double thresh;
    int* curr_node = new int[n];
    int visited_size = (((n * 2) - 1) >> 3) + 1;
    unsigned char* visited = new unsigned char[visited_size];

    //here was some memory error catdch code

    lower_idx = 0;
    upper_idx = n - 1;
    while(upper_idx - lower_idx > 1) {
        i = (lower_idx + upper_idx) >> 1;
        thresh = MC[i];

        memset(visited, 0, visited_size); //muss durch C++ func ersetz werden!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        nc = 0;
        k = 0;
        curr_node[0] = 2 * n - 2;

        while(k >= 0) {
            root = curr_node[k] - n;
            i_lc = int (Z[root][0]);
            i_rc = int (Z[root][1]);

            if (MC[root] <= thresh) {  // this subtree forms a cluster
                nc += 1;
                if (nc > max_nc) { // illegal
                    break;
                k -= 1;
                }
                set_visited(visited, i_lc);
                set_visited(visited, i_rc);
                continue;
            }
            if (!is_visited(visited, i_lc)) {
                set_visited(visited, i_lc);
                if (i_lc >= n) {
                    k += 1;
                    curr_node[k] = i_lc;
                    continue;
                }
                else{  // singleton cluster
                    nc += 1;
                    if (nc > max_nc) break;
                }
            }

            if (!is_visited(visited, i_rc)) {
                set_visited(visited, i_rc);
                if (i_rc >= n){
                    k += 1;
                    curr_node[k] = i_rc;
                    continue;
                    }
                else{  // singleton cluster
                    nc += 1;
                    if (nc > max_nc ){
                        break;
                    }
                }
            k -= 1;
            }

        if (nc > max_nc) lower_idx = i;
        else upper_idx = i;
        
        }

    delete visited;
    cluster_monocrit(Z, MC, T, MC[upper_idx], n);
    }
}
