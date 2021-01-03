#include "main.h"

char *make_file_name(char *stem, int index, char *ext, char *fname)
{
    int len = snprintf(NULL, 0, "%s%d%s", stem, index, ext);
    fname = realloc(fname, sizeof(char) * (len + 1));
    snprintf(fname, len + 1, "%s%d%s", stem, index, ext);
    return fname;
}

char *make_file_name2(char *stem, int i1, char *in, int i2, char *ext, char *fn)
{
    int len = snprintf(NULL, 0, "%s%d%s%d%s", stem, i1, in, i2, ext);
    fn = realloc(fn, sizeof(char) * (len + 1));
    snprintf(fn, len + 1, "%s%d%s%d%s", stem, i1, in, i2, ext);
    return fn;
}

char *make_file_name3(char *stem, int i1, char *in1, int i2, char *in2, int i3, char *ext, char *fn)
{
    int len = snprintf(NULL, 0, "%s%d%s%d%s%d%s", stem, i1, in1, i2, in2, i3, ext);
    fn = realloc(fn, sizeof(char) * (len + 1));
    snprintf(fn, len + 1, "%s%d%s%d%s%d%s", stem, i1, in1, i2, in2, i3, ext);
    return fn;
}

int read_data(char *file_name, int n_genes, int sample_size, DATATYPE_VALUE *data)
{
    int j, k, status;
    FILE *file = fopen(file_name, "r");

    if (file == NULL)
    {
        fprintf(stderr, "ERROR: file %s cannot be open to read.\n", file_name);
        exit(EXIT_FAILURE);
    }

    for (j = 0; j < n_genes; j++)
        for (k = 0; k < sample_size; k++)
            status = fscanf(file, "%g", &data[j * sample_size + k]);
    fclose(file);
    return EXIT_SUCCESS;
}

int extract_column(unsigned int n, unsigned int m, DATATYPE_VALUE data[n * m],
                   unsigned int c, DATATYPE_VALUE column[n])
{
    unsigned int i;
    for (i = 0; i < n; i++)
        column[i] = data[i * m + c];
    return EXIT_SUCCESS;
}

int get_mean(float *data, int n_genes, int sample_size, float *mean_value)
{
    for (int nr_i = 0; nr_i < n_genes; nr_i++)
    {
        float sum = 0;
        for (int nr_j = 0; nr_j < sample_size; nr_j++)
        {
            sum += data[nr_i * sample_size + nr_j];
        }
        mean_value[nr_i] = sum / sample_size;
    }
}

int get_rank(float *mean_value, int n_genes, float *rank)
{
    for (int index_i = 0; index_i < n_genes; index_i++)
    {
        float count_rank_n = 0;
        for (int index_j = 0; index_j < n_genes; index_j++)
        {

            if (mean_value[index_j] < mean_value[index_i])
                count_rank_n++;
        }
        rank[index_i] = count_rank_n + 1;
    }
}

int extract_data(float *data,
                 int n_genes,
                 int sample_size,
                 int from_sample,
                 int needed_samples,
                 float *extracted_data)
{
    for (int i = 0; i < n_genes; i++)
    {
        for (int col = 0; col < needed_samples; col++)

        {
            extracted_data[i * needed_samples + col] = data[i * sample_size + (col + from_sample)];
        }
    }
}


int get_normal_rank(float *data,
                           int n_genes,
                           int sample_size,
                           float *normal_rank)
{
    float *data_mean = malloc(sizeof(float) * n_genes);
    get_mean(data, n_genes, sample_size, data_mean);
    get_rank(data_mean, n_genes, normal_rank);
}
//get case_rank:get rank for our dealing sample(column)
int get_case_rank(float *column, int n_genes, int *case_rank)
{
    for (int index_i = 0; index_i < n_genes; index_i++)
    {
        int count_rank = 0;
        for (int index_j = 0; index_j < n_genes; index_j++)
        {
            if (column[index_i] > column[index_j])
            {
                count_rank++;
            }
        }
        case_rank[index_i] = count_rank + 1;
    } //end of get case_rank
}

int compare(const void *a, const void *b)
{
    float *x = (float *)a;
    float *y = (float *)b;
    if (*x < *y)
        return -1;
    else if (*x > *y)
        return 1;
    return 0;
}

float get_threshod_rd(UB8 **pairs,
                      int n_threads,
                      DATATYPE_GENESIZE count_pairs[n_threads],
                      int *exceptions[n_threads],
                      int max1,
                      float *normal_rank)
{

    //get number of stable_pairs
    int num_pairs = 0;
    for (int i = 0; i <= max1; i++)
    {
        num_pairs += exceptions[0][i];
    }
    fprintf(stdout, "number of stable pairs: %d\n", num_pairs);

    /*filter stable pairs: remove pairs with lower %5 rank difference*/
    //step1:calculate rank differences
    float *rank_differences = malloc(sizeof(float) * num_pairs);
    int rank_index = 0;
    for (int indexi = 0; indexi < n_threads; indexi++)
    {
        for (int indexj = 0; indexj < count_pairs[indexi]; indexj++)
        {
            UB8 current = pairs[indexi][indexj];
            if (current != 0)
            {
                unsigned int vi = ((unsigned int)(current >> 22)) & 0x1FFFFF; //clear left bits
                unsigned int vj = ((unsigned int)(current >> 43)) & 0x1FFFFF;

                float rd_tmp;
                rd_tmp = fabs(normal_rank[vi] - normal_rank[vj]);
                rank_differences[rank_index] = rd_tmp;
                rank_index++;
            }
        }
    }

    qsort(rank_differences, num_pairs, sizeof(rank_differences[0]), compare);
    int threshod_rd_index = (int)(num_pairs * 0.15);
    fprintf(stdout, "threshod_rd_index: %d\n", threshod_rd_index);

    return rank_differences[threshod_rd_index];
}


/******************************************************************************
 * Sample Type : individual samples, each column is a sample
 * Job types 1 and 2: select concordant and reversed genes pairs
 *  Input: Data Files, FDR Levels (or Max Exception Number), max_equals, mode
 * Output: concordant and reversed gene pairs 
 ******************************************************************************/

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
                                 char algorithm)
{
    int i1, i2, end1, step1;
    switch (pair_mode)
    {
    case 0:
        end1 = 1;
        step1 = 1;
        break;

    case 1:
        end1 = n_files - 1;
        step1 = 1;
        break;

    case 2:
        end1 = n_files - 1;
        step1 = 2;
        break;

    } //end switch

    for (i1 = 0; i1 < end1; i1 += step1)
    {
        /* load data1 into memory */
        int size1 = sample_sizes[i1];
        char *file1 = files[i1];
        DATATYPE_VALUE *data1 = malloc(sizeof(DATATYPE_VALUE) * n_genes * size1);
        if (data1 == NULL)
        {
            fprintf(stderr, "ERROR: memory is not allocated for data file %d,  %s.\n", i1, file1);
            break;
            //exit(EXIT_FAILURE);
        }
        read_data(file1, n_genes, size1, data1);

        /* get normal rank */
        float *normal_rank = malloc(sizeof(float) * n_genes);
        if (algorithm == 0 )
        {
            //get normal rank
            fprintf(stdout, "# Get normal rank...\n");
            get_normal_rank(data1, n_genes, size1, normal_rank);
        }
        fprintf(stdout, "# Normal rank done...\n");
        /* max exceptions */
        float fdr1 = fdr[i1];
        int mode1 = 0, max1;
        if (fdr1 >= 1. || fdr1 == 0.)
        { //selection based on max exceptions
            max1 = (int)fdr1;
            max1 = max1 < 0 ? 0 : max1;
            mode1 = 1;
        }
        else
        {
            max1 = (size1 - 1) / 2; //WARNING: boundary issue
            mode1 = 0;
        } //end if...else...

        /* select stable pairs in sample i1 */
        int n_threads;
#pragma omp parallel shared(n_threads)
        {
#pragma omp single
            n_threads = omp_get_num_threads();
        }

       
        UB8 **pairs = NULL;
        int **exceptions1 = NULL;
        DATATYPE_GENESIZE *count_pairs = NULL;
        pairs = malloc(sizeof(UB8 *) * n_threads);
        exceptions1 = malloc(sizeof(int *) * n_threads);
        count_pairs = malloc(sizeof(int *) * n_threads);

        fprintf(stdout, "# Start select stable pairs.........\n");

        stable_pairs_one2(n_genes, size1, data1, max1, max_equals,
                          n_threads, pairs, count_pairs, exceptions1);
        fprintf(stdout, "# Select stable pairs done..........\n");

        struct FDR result1;
        int j, k;

        /* FDR control and print out stable pairs */
        /* FDR control mode*/
        if (mode1 == 0)
        {
            result1 = bh_ctrl(size1, max1, exceptions1[0], n_genes * (n_genes - 1) / 2, fdr1);
            max1 = result1.index;
        }
        else
        {
            result1 = bh_eval(size1, max1, exceptions1[0], n_genes * (n_genes - 1) / 2);
            result1.p = bino_p(size1, result1.index);
        }

        /* print ouf statistics */
        fprintf(stdout, "# Data Set %d - Significantly Stable Gene Pairs \n", i1);
        fprintf(stdout, "#  FDR control level or max exception numbers: %g\n", fdr1);
        fprintf(stdout, "#  FDR actual level: %g\n", result1.value);
        fprintf(stdout, "#  Max. exception numbers: %d\n", max1);
        fprintf(stdout, "#  Max. binomial test P value: %g\n", result1.p);
        for (j = 0; j <= max1; j++)
            fprintf(stdout, "#  Number of pairs with %d exceptions: %d\n", j, exceptions1[0][j]);

        char *out_name1 = NULL, *out_name3 = NULL, *out_name4 = NULL;
        FILE *out_file1 = NULL, *out_file3 = NULL, *out_file4 = NULL;
        
        /* write stable pairs to files */
        for (j = 0; j < n_threads; j++)
            for (k = 0; k < count_pairs[j]; k++)
            {
                UB8 current = pairs[j][k];
                if (current == 0)
                    continue;

                unsigned char state;
                unsigned char e1 = ((unsigned char)(current >> 4)) & 0b00000001;
                //c1 表示nx，扩展2位后需要unsigned short int才放得下
                // unsigned char c1 = (unsigned char)(current >> 6);
                unsigned short int c1 = (unsigned int)(current >> 6);

                state = ((c1 <= max1 && e1 == 1) ? 1 : 0);
                state = state | (((size1 - c1 <= max1 && e1 == 1) ? 1 : 0) << 1);

                // unsigned int vi = ((unsigned int)(current >> 22)) & 0x1FFFFF; //clear left bits
                // unsigned int vj = ((unsigned int)(current >> 43)) & 0x1FFFFF;
                //由于nx,cx的位数扩充以及gi,gj的位数缩减，移位数目和与操作的1数目需要相应修改
                unsigned int vi = ((unsigned int)(current >> 26)) & 0x7FFFF; //clear left bits
                unsigned int vj = ((unsigned int)(current >> 45)) & 0x7FFFF;
                

                if (state == 0)
                {
                    current = 0;
                }
                else
                {
                    current = (UB8)state;
                    //由于nx,cx的位数扩充以及gi,gj的位数缩减，移位数目需要相应修改
                    current = current | ((UB8)vi << 26);
                    current = current | ((UB8)vj << 45);
                }
                pairs[j][k] = current;

            } //end for

        /*loop2 control*/
        int begin2, step2;
        switch (pair_mode)
        {
        case 0:
        case 1:
            begin2 = i1 + 1;
            step2 = 1;
            break;

        case 2:
            begin2 = 1;
            step2 = 2;
            break;

        } //end switch
        /* loop over the rest of the data files */
        for (i2 = begin2; i2 < n_files; i2 += step2)
        {
            if (VERBOSE)
                fprintf(stdout, "#Processing data file %d...\n", i2);
            /* load data2 into memory */
            int size2 = sample_sizes[i2];
            char *file2 = files[i2];
            DATATYPE_VALUE *data2 = malloc(sizeof(DATATYPE_VALUE) * n_genes * size2);
            if (data2 == NULL)
            {
                fprintf(stderr, "ERROR: memory is not allocated for data file %d,  %s.\n", i2, file2);
                break;
                //exit(EXIT_FAILURE);
            }
            read_data(file2, n_genes, size2, data2);

            int i3;
            DATATYPE_VALUE *column;
            column = malloc(sizeof(DATATYPE_VALUE) * n_genes);

            /* table to store all gene directions */
            // bits 0,1 for state in the RandComp V2.0 algo. 00=0, flat; 10=2, up; 01=1, down
            // bits 2,3 for state in the original algo. 00xx, flat; 10xx=4+, up; 01xx=4+, down
            unsigned char *genes_result;
            genes_result = malloc(sizeof(unsigned char) * n_genes * size2);
            for (j = 0; j < n_genes; j++)
                for (k = 0; k < size2; k++)
                    genes_result[j * size2 + k] = 0;

            for (i3 = 0; i3 < size2; i3++)
            {
                fprintf(stdout, "# Processing column %d...\n", i3);

                /* column i3*/
                extract_column(n_genes, size2, data2, i3, column);

                /* consistent pairs */
                fprintf(stdout, "ngenes = %d\n", n_genes);
                stable_pairs_ind(n_genes, n_threads, pairs, count_pairs, column);
                

                int count_same = 0, count_reverse = 0;

                /* write stable pairs to files */
                for (j = 0; j < n_threads; j++)
                    for (k = 0; k < count_pairs[j]; k++)
                    {
                        UB8 current = pairs[j][k];
                        unsigned char state = ((unsigned char)current) & 0b1111;
                        //由于nx,cx的位数扩充以及gi,gj的位数缩减，移位数目需要相应修改
                        unsigned int vi = ((unsigned int)(current >> 26)) & 0x7FFFF; //clear left bits
                        unsigned int vj = ((unsigned int)(current >> 45)) & 0x7FFFF;

                        /* keep only the consistent pairs */
                        switch (state)
                        {
                        case 5: // 1+4
                            if (job == 1 || VERBOSE)
                                fprintf(out_file3, "%d\t%d\n", vi, vj);
                            count_same += 1;
                            break;
                        case 10: // 2+8
                            if (job == 1 || VERBOSE)
                                fprintf(out_file3, "%d\t%d\n", vj, vi);
                            count_same += 1;
                            break;
                        case 6: // 2+4
                            if (job == 1 || VERBOSE)
                                fprintf(out_file4, "%d\t%d\n", vj, vi);
                            count_reverse += 1;
                            break;
                        case 9: // 1+8
                            if (job == 1 || VERBOSE)
                                fprintf(out_file4, "%d\t%d\n", vi, vj);
                            count_reverse += 1;
                            break;
                        } //end switch
                    }     //end for
                fprintf(stdout, "#  Number of concordant pairs between data set %d column %d and %d: %d\n", i2, i3, i1, count_same);
                fprintf(stdout, "#  Number of reversed pairs between data set %d column %d and %d: %d\n", i2, i3, i1, count_reverse);
                /* clean up */

                
                if (algorithm == 0)
                {
                    printf("# running with the InDEA method..........\n");
                    //state for each gene, not pairs
                    struct gene_state **genes; //an array of pointers to struct
                    //allocate an array of n pointers to struct gene_states
                    genes = (struct gene_state **)malloc(sizeof(struct gene_state *) * n_genes);
                    //initialization
                    for (j = 0; j < n_genes; j++)
                        genes[j] = (struct gene_state *)NULL;

                    //get case rank
                    fprintf(stdout, "Get case rank for column %d...\n", i3);
                    int *case_rank = malloc(sizeof(int) * n_genes * size2);
                    get_case_rank(data2, n_genes, case_rank);

                    //call
                    fprintf(stdout, "#  Filter dysregulated genes...\n");
                    filter_gene_dirs_new(n_genes, n_threads, pairs, count_pairs,
                                         genes, alpha, max_cycles, max_threshold, normal_rank, case_rank);

                    /* write dysregulated genes into files*/
                    for (j = 0; j < n_genes; j++)
                        if (genes[j] != NULL)
                        {
                            genes_result[j * size2 + i3] = genes[j]->state;
                            free(genes[j]);
                        }
                    free(genes);
                } //end if algo.

            } //end for i3

            /*Print out genes_result */
            /* open file to write dysregulated genes */
            out_name1 = malloc(sizeof(char *) * 10);
            out_name1 = make_file_name(GENE_STATE, i2, EXT, out_name1);
            out_file1 = fopen(out_name1, "w");
            if (out_file1 == NULL)
            {
                out_name1 = "stderr";
                out_file1 = stderr;
            }
            fprintf(stdout, "#  Writing gene states to file: %s\n", out_name1);
            free(out_name1);
            for (j = 0; j < n_genes; j++)
            {
                for (k = 0; k < size2; k++)
                    fprintf(out_file1, "%d\t", genes_result[j * size2 + k]);
                fprintf(out_file1, "\n");
            }
            if (out_file1 != stderr)
                fclose(out_file1);

        } //end inner for
        /*clean up*/
        /* clean up*/
        for (j = 0; j < n_threads; j++)
            free(pairs[j]);
        free(pairs);
        free(count_pairs);
        free(data1);

    } //end outer for
    return EXIT_SUCCESS;

} // end of select_stable_pairs_ind


int main(int argc, char **argv)
{

    /* default values */
    //algorithm: 0---InDEA;
    char algorithm = 0;

    VERBOSE = false;
    EPSILON = FLT_EPSILON;
    int max_equals = 256;
    float fdr_level = 0.05;
    char sample_type = 1;
    char job_type = 2;
    char pair_mode = 0;
    int times = 1;

    int max_cycles = 128;
    int ori_cycles = 2;
    int conv_threshold = 50;
    char *cname = "";
    int n_changes = 0;

    int n_files, n_genes, sample_sizes[MAX_FILES];
    float fdr_levels[MAX_FILES];
    char *files[MAX_FILES];

    /* check number of files */
    if (optind < argc)
    {
        n_files = atoi(argv[optind++]);
        if (n_files <= 0)
        {
            fprintf(stderr, "ERROR: n_files should be an non-negative integer.\n");
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        fprintf(stderr, "ERROR: too few arguments.\n\n");
        exit(EXIT_FAILURE);
    }

    /* check number of genes */
    if (optind < argc)
    {
        n_genes = atoi(argv[optind++]);
        if (n_genes <= 0)
        {
            fprintf(stderr, "ERROR: n_genes should be an non-negative integer.\n\n");
            exit(EXIT_FAILURE);
        }
    }
    else
    {
        fprintf(stderr, "ERROR: too few arguments.\n\n");
        exit(EXIT_FAILURE);
    }
    fprintf(stdout, "n_genes = %d\n", n_genes);

    int i;
    char *file_name;

    /* check file names  */
    for (i = 0; i < n_files; i++)
    {
        if (optind < argc)
        {
            file_name = argv[optind++];
            if (access(file_name, R_OK) == -1)
            {
                fprintf(stderr, "ERROR: data file %s does not exist or not readable.\n", file_name);
                exit(EXIT_FAILURE);
            }
            else
                files[i] = file_name;
        }
        else
        {
            fprintf(stderr, "ERROR: number of file names does not match with n_files.\n\n");
            exit(EXIT_FAILURE);
        } //end if...else...

    } //end for

    /* check sample size */
    int size_current;
    for (i = 0; i < n_files; i++)
    {
        if (optind < argc)
        {
            size_current = atoi(argv[optind++]);
            if (size_current <= 0 || size_current > MAX_SAMPLESIZE)
            {
                fprintf(stderr, "ERROR: sample size %d is too small(<=0) or too big(>65535).\n", size_current);
                exit(EXIT_FAILURE);
            }
            else
                sample_sizes[i] = size_current;
        }
        else
        {
            fprintf(stderr, "ERROR: number of sample sizes does not match with n_files.\n\n");
            exit(EXIT_FAILURE);
        } //end if...else...
    }     //end for

    /* check fdr_level or exception number */
    float fdr_current;
    i = 0;
    while (optind < argc && i < n_files)
    {
        fdr_current = atof(argv[optind++]);
        if (fdr_current < 0.)
        {
            fprintf(stderr, "ERROR: FDR level or exception threshold  %f should be non-negative.\n", fdr_current);
            exit(EXIT_FAILURE);
        }
        else
            fdr_levels[i++] = fdr_current;
    } //end while
    while (i < n_files)
        fdr_levels[i++] = DEFAULT_FDR;

    /*do our job*/
    select_consistent_pairs_ind2(n_files, files, sample_sizes, fdr_levels, n_genes, max_equals,
                                 job_type, pair_mode, fdr_level, max_cycles, ori_cycles, conv_threshold, algorithm);

    return EXIT_SUCCESS;
}