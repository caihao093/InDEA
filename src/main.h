#include "indea.h"
#include <getopt.h>
#include <unistd.h>

#define STABLE "stable_pairs_"
#define CONCORDANT "concordant_pairs_"
#define REVERSED "reversed_pairs_"
#define UP "up_regulated_"
#define DOWN "down_regulated_"
#define GENE_STATE "gene_state_"
#define SIMULATION "simulation_"
#define SIMILARITY "similarity_"
#define BAYESIAN "bayesian_"
#define SAMPLED_GENES "sampled_genes_"
#define EXT ".dat"

float MEMORY_USAGE = -1.;

char *make_file_name(char *stem, int index, char *ext, char *fname);
char *make_file_name2(char *stem, int i1, char *in, int i2, char *ext, char *fn);
char *make_file_name3(char *stem, int i1, char *in1, int i2, char *in2, int i3, char *ext, char *fn);
int read_data(char *file_name, int n_genes, int sample_size, DATATYPE_VALUE *data);
int extract_column(unsigned int n, unsigned int m, DATATYPE_VALUE data[n * m],
                   unsigned int c, DATATYPE_VALUE column[n]);
//get_rank:
//输入：一个浮点型，长度为n_genes的数组及其长度
//输出：一个整型数组，记录输入数组的秩次，长度为n_genes
int get_rank(float *mean_value, int n_genes, float *rank);

//get_mean:输入一个n_genes行，sample_size列的矩阵及其维度
//输出每一行的均值组成的长度为n_genes的均值数组
int get_mean(float *data, int n_genes, int sample_size, float *mean_value);

//extract_data:
//输入：一个有n_genes行sample_size列的浮点型数据矩阵data及其维度，
////截取样本数据的起始样本from_sample，要截取的样本数needed_samples
//输出：n_genes行needed_samples列，从data中截取出来相应位置的浮点型数据矩阵
int extract_data(float *data,
                 int n_genes,
                 int sample_size,
                 int from_sample,
                 int needed_samples,
                 float *extracted_data);

int get_normal_rank(float *data,
                           int n_genes,
                           int sample_size,
                           float *normal_rank);
int get_case_rank(float *column, int n_genes, int *case_rank);
float get_threshod_rd(UB8 **pairs,
                      int n_threads,
                      DATATYPE_GENESIZE count_pairs[n_threads],
                      int *exceptions[n_threads],
                      int max1,
                      float *normal_rank);


/******************************************************************************
 * Sample Type : individual samples, each column is a sample
 * Job types 1 and 2: select concordant and reversed genes pairs
 *  Input: Data Files, FDR Levels (or Max Exception Number), max_equals, mode
 * Output: concordant and reversed gene pairs 
 ******************************************************************************/

/* New Output Format */
int select_consistent_pairs_ind2(int n_files,
                                 char *files[n_files],
                                 int sample_sizes[n_files],
                                 float fdr[n_files],
                                 int n_genes,
                                 int max_equals,
                                 char job,
                                 char pair_mode,
                                 float alpha,
                                 int max_cycles,
                                 int ori_cycles,
                                 int max_threshold,
                                 char algorithm);
