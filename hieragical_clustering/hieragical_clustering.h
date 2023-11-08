

int is_visited(unsigned char* bitset, int i);

void set_visited(unsigned char* bitset, int i);

void cluster_maxclust_monocrit(double** Z, double* MC, int* T,
                                    int n, int max_nc);

void cluster_monocrit(double** Z, double* MC, int* T,
                            double cutoff, int n);                                    

