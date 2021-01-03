#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <omp.h>
#include <math.h>
#include "statistics.h"

#define DATATYPE_VALUE float
#define DATATYPE_SAMPLESIZE unsigned int
#define DATATYPE_GENESIZE unsigned int
#define CONV_THRESHOLD 50
#define MEMBLOCK 1024
#define BIGNINT -1000000
#define MASK_FDR 0b00000001
#define MASK_FILTER 0b00000010
#define MAX_FILES 255
#define MAX_SAMPLESIZE 1023
#define DEFAULT_FDR 0.05


//new UB8,increase max normal sample size to 1023
//Use a 64-bit int to represent i, j, nx, cx and state
// Bits 0-6 for state, 0-3 for '<' and '>' test, 2-3 for '=' test
// 0 - if   count < threshold in normal group (sample 1)
// 1 - if m-count < threshold in normal group (sample 1)
// 2 - if   count < threshold in case group (sample 2)
// 3 - if m-count < threshold in case group (sample 2)
// 4 - if count0 < threshold in normal group (sample 1)
// 5 - if count0 < threshold in case group (sample 2)
//  6-15, nx, count of <|> (smaller one) in normal group, capacity 2^10-1 = 1023
// 16-25, cx, ibid in case group, sample size 1023*2=2046
// 26-44, gi, gene index, capacity, 2^19-1 = 524287 gene number
// 45-63, gj, gene index 
typedef unsigned long int UB8;

bool VERBOSE;
double EPSILON;

//for one sample only
struct pair
{
    DATATYPE_GENESIZE h, l;
    DATATYPE_SAMPLESIZE count;
};


struct gene_state
{
    unsigned char state;  //rightmost two bits for in the same or reversed pairs
    DATATYPE_GENESIZE ng; // control group, greater than the current gene
    DATATYPE_GENESIZE nl; // control group, less than the current gene
    DATATYPE_GENESIZE cg; // case group, greater than the current gene
    DATATYPE_GENESIZE cl; // case group, less than the current gene
    double p;             //hyergeometric test value
};


struct gene_state2
{
    unsigned char state;  //rightmost two bits for in the same or reversed pairs
    DATATYPE_GENESIZE ng; // control group, greater than the current gene
    DATATYPE_GENESIZE nl; // control group, less than the current gene
    float g2l;            // reversed pair, the number of genes which are greater than the current gene but lower in the case group
    float l2g;            // similar as above, but in the reversed order
    double p;             //hyergeometric test value
};

//Stable gene-pair ordering for one sample
int stable_pairs_one(DATATYPE_GENESIZE n,        //number of genes
                     DATATYPE_SAMPLESIZE m,      //sample size
                     DATATYPE_VALUE data[n * m], //data matrix
                     DATATYPE_SAMPLESIZE max,    //threshold, max number of exceptions
                     DATATYPE_GENESIZE max0,     //max allowed number of equal pairs
                     int nthreads,
                     struct pair *pairs[nthreads],
                     DATATYPE_GENESIZE count_pairs[nthreads],
                     int *exceptions[nthreads] //for FDR control, number of pairs with the specified number of exceptions
);

/* Stable gene-pair ordering for one  sample, using UB8 */
int stable_pairs_one2(DATATYPE_GENESIZE n,        //number of genes
                      DATATYPE_SAMPLESIZE m,      //sample size
                      DATATYPE_VALUE data[n * m], //data matrix
                      DATATYPE_SAMPLESIZE max,    //threshold, max number of exceptions
                      DATATYPE_GENESIZE max0,     //max allowed number of equal pairs
                      int nthreads,
                      UB8 *pairs[nthreads],
                      DATATYPE_GENESIZE count_pairs[nthreads],
                      int *exceptions[nthreads] //for FDR control, number of pairs with the specified number of exceptions
);
/* do intersect genepairs*/
int stable_pairs_one3(DATATYPE_GENESIZE n,        //number of genes
                      DATATYPE_SAMPLESIZE m,      //sample size
                      DATATYPE_VALUE data[n * m], //data matrix
                      DATATYPE_SAMPLESIZE max,    //threshold, max number of exceptions
                      DATATYPE_GENESIZE max0,     //max allowed number of equal pairs
                      int nthreads,
                      UB8 *pairs[nthreads],
                      DATATYPE_GENESIZE count_pairs[nthreads],
                      int *exceptions[nthreads] //for FDR control, number of pairs with the specified number of exceptions
);
/* Stable gene-pair ordering for two samples,
 * the second sample is one column, individual sample,
 * stable pairs for the first sample is already given*/
int stable_pairs_ind(DATATYPE_GENESIZE n, //number of genes
                     int nthreads,
                     UB8 *pairs[nthreads],
                     DATATYPE_GENESIZE count_pairs[nthreads],
                     DATATYPE_VALUE column[n] //data matrix
);

int filter_gene_orig(DATATYPE_GENESIZE n, //number of genes
                     int nthreads,
                     UB8 *pairs[nthreads],
                     DATATYPE_GENESIZE count_pairs[nthreads],
                     struct gene_state2 *states[n],
                     double alpha, //FDR alpha level for regulation direction
                     int max_cycles,
                     int conv_threshold);

int filter_gene_orig_new(DATATYPE_GENESIZE n, //number of genes
                    int nthreads,
                    UB8 *pairs[nthreads],
                    DATATYPE_GENESIZE count_pairs[nthreads],
                    struct gene_state2 *states[n],
                    double alpha, //FDR alpha level for regulation direction
                    int max_cycles,
                    int conv_threshold,
                    float *normal_rank,
                    int *case_rank);

int filter_gene_dirs(DATATYPE_GENESIZE n, //number of genes
                     int nthreads,
                     UB8 *pairs[nthreads],
                     DATATYPE_GENESIZE count_pairs[nthreads],
                     struct gene_state *states[n],
                     double alpha, //FDR alpha level for regulation direction
                     int max_cycles,
                     int conv_threshold);

int filter_gene_dirs_new(DATATYPE_GENESIZE n, //number of genes
                         int nthreads,
                         UB8 *pairs[nthreads],
                         DATATYPE_GENESIZE count_pairs[nthreads],
                         struct gene_state *states[n],
                         double alpha, //FDR alpha level for regulation direction
                         int max_cycles,
                         int conv_threshold,
                         float *normal_rank,
                         int *case_rank);

