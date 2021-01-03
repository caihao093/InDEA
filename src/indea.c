#include "indea.h"

#pragma omp declare simd
int less(DATATYPE_VALUE a, DATATYPE_VALUE b)
{
    if (a < b)
        return 1;
    else
    {
        if (fabs(a - b) < EPSILON)
            return random_bin();
        else
            return 0;
    }
}

#pragma omp declare simd
int equal(DATATYPE_VALUE a, DATATYPE_VALUE b)
{
    if (fabs(a - b) < EPSILON)
        return 1;
    else
        return 0;
}

#define MEMCHECK(var)                                                               \
    if (var == NULL)                                                                \
    {                                                                               \
        fprintf(stderr, "# Memory Allocation Error! More Memory May Be Needed!\n"); \
        exit(-1);                                                                   \
    }

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
)
{
    random_init();

#pragma omp parallel
    {

        int ithread = omp_get_thread_num();
        struct pair *new_array;

        int count_memblock = 0;
        count_pairs[ithread] = 0;
        int count, count0;
        DATATYPE_SAMPLESIZE k;
        exceptions[ithread] = malloc(sizeof(int) * (max + 1));
        pairs[ithread] = malloc(sizeof(struct pair) * MEMBLOCK);
        for (k = 0; k < max + 1; k++)
            exceptions[ithread][k] = 0;

        unsigned char state;
        DATATYPE_GENESIZE i;
#pragma omp for
        for (i = 0; i < n; i++)
        {
            DATATYPE_GENESIZE j;
            for (j = i + 1; j < n; j++)
            {
                //fprintf(stderr, "thread=%d, counts=%d\n",ithread, count_pairs[ithread]);
                if (count_pairs[ithread] >= MEMBLOCK * count_memblock)
                {
                    count_memblock += 1;
                    new_array = realloc(pairs[ithread], sizeof(struct pair) * MEMBLOCK * count_memblock);
                    if (new_array) //avoid memory leak
                        pairs[ithread] = new_array;
                    else
                    {
                        fprintf(stderr, "# Memory is not allocated on thread %d.\n", ithread);
                        exit(-1);
                    }
                } //end if

                count = 0;
                count0 = 0;

#pragma omp simd
                for (k = 0; k < m; k++) // count1 [i]<[j]
                    count += less(data[i * m + k], data[j * m + k]);
#pragma omp simd
                for (k = 0; k < m; k++) // count1 [i]<[j]
                    count0 += equal(data[i * m + k], data[j * m + k]);

                state = ((count <= max && count0 <= max0) ? 1 : 0);
                state = state | (((m - count <= max && count0 <= max0) ? 1 : 0) << 1);
                if (count <= max && count0 <= max0) // [i]>[j]
                {
                    //insert to the pairs
                    pairs[ithread][count_pairs[ithread]].h = i;
                    pairs[ithread][count_pairs[ithread]].l = j;
                    pairs[ithread][count_pairs[ithread]].count = count;
                    count_pairs[ithread] += 1;
                    exceptions[ithread][count] += 1;
                }                                            //end if
                else if (m - count <= max && count0 <= max0) // [i] <[j]
                {
                    //insert to the pairs
                    pairs[ithread][count_pairs[ithread]].h = j;
                    pairs[ithread][count_pairs[ithread]].l = i;
                    pairs[ithread][count_pairs[ithread]].count = m - count;
                    count_pairs[ithread] += 1;
                    //could be improved
                    exceptions[ithread][m - count] += 1;
                } //end if
            }     //end for j
        }         //end for i
    }             //end parallel
    //sum all the counts into exceptions[0]

    DATATYPE_SAMPLESIZE k;
    int ithread;
    for (k = 0; k < max + 1; k++)
        for (ithread = 1; ithread < nthreads; ithread++)
            exceptions[0][k] += exceptions[ithread][k];

    return EXIT_SUCCESS;
}

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
)
{
    random_init();

#pragma omp parallel
    {

        int ithread = omp_get_thread_num();
        UB8 *new_array = NULL, current = 0;
        int count_memblock = 0;
        count_pairs[ithread] = 0;
        // unsigned char count, count0;
        //由于nx,cx 扩充了两位，所以count和count0需要用short int才放得下
        unsigned short int count, count0;

        DATATYPE_SAMPLESIZE k;
        exceptions[ithread] = malloc(sizeof(int) * (max + 1));
        pairs[ithread] = malloc(sizeof(UB8) * MEMBLOCK);
        for (k = 0; k < max + 1; k++)
            exceptions[ithread][k] = 0;

        unsigned char state;
        DATATYPE_GENESIZE i;
#pragma omp for
        for (i = 0; i < n; i++)
        {
            DATATYPE_GENESIZE j;
            for (j = i + 1; j < n; j++)
            {
                if (count_pairs[ithread] >= MEMBLOCK * count_memblock)
                {
                    count_memblock += 1;
                    new_array = realloc(pairs[ithread], sizeof(UB8) * MEMBLOCK * count_memblock);
                    if (new_array) //avoid memory leak
                        pairs[ithread] = new_array;
                    else
                    {
                        fprintf(stderr, "# Memory is not allocated on thread %d.\n", ithread);
                        exit(-1);
                    }
                } //end if

                count = 0;
                count0 = 0;

#pragma omp simd
                for (k = 0; k < m; k++) // count [i]<[j]
                    count += less(data[i * m + k], data[j * m + k]);
#pragma omp simd
                for (k = 0; k < m; k++) // count0 [i]==[j]
                    count0 += equal(data[i * m + k], data[j * m + k]);

                state = ((count <= max && count0 <= max0) ? 1 : 0);
                state = state | ((((m - count) <= max && count0 <= max0) ? 1 : 0) << 1);

                if (state != 0)
                {
                    if (state == 1)
                        exceptions[ithread][count] += 1;
                    else
                        exceptions[ithread][m - count] += 1;

                    state = state | (((count0 <= max0) ? 1 : 0) << 4);
                    //insert to the pairs
                    current = (UB8)state;
                    current = current | ((UB8)count << 6);

                    // current = current | ((UB8)i << 22);
                    // current = current | ((UB8)j << 43);
                    // 由于gi, gj都由原来的21位缩减到19位以及nx, cx位数的扩充，所以需要改变相应移位数
                    current = current | ((UB8)i << 26);
                    current = current | ((UB8)j << 45);

                    pairs[ithread][count_pairs[ithread]] = current;
                    count_pairs[ithread] += 1;
                }
                else
                    pairs[ithread][count_pairs[ithread]] = 0;

            } //end for j
        }     //end for i
    } //end parallel

    /* sum all the counts into exceptions[0] */
    DATATYPE_SAMPLESIZE k;
    int ithread;
    for (k = 0; k < max + 1; k++)
        for (ithread = 1; ithread < nthreads; ithread++)
            exceptions[0][k] += exceptions[ithread][k];

    return EXIT_SUCCESS;
}


/* Stable gene-pair ordering for two samples,
 * the second sample is one column, individual sample,
 * stable pairs for the first sample is already given*/
int stable_pairs_ind(DATATYPE_GENESIZE n, //number of genes
                     int nthreads,
                     UB8 *pairs[nthreads],
                     DATATYPE_GENESIZE count_pairs[nthreads],
                     DATATYPE_VALUE column[n] //data matrix
)
{
#pragma omp parallel
    {
        int ithread = omp_get_thread_num();
        unsigned int index, i, j;
        unsigned char state;
        UB8 current;
        for (index = 0; index < count_pairs[ithread]; index++)
        {
            current = pairs[ithread][index];
            state = ((unsigned char)current) & 0b11;

            // 由于gi, gj都由原来的21位缩减到19位以及nx, cx位数的扩充，所以需要改变相应移位数以及与运算的1数目
            i = ((unsigned int)(current >> 26)) & 0x7FFFF; //clear left bits
            j = ((unsigned int)(current >> 45)) & 0x7FFFF;

            // fprintf(stdout,"i = %d\tj = %d\n",i,j);
            switch (state)
            {
            case 1:                                  // [i] > [j] in normal
                if (less(column[i], column[j]) == 1) // [i] < [j] in case
                    state = 9;                       //1 + 8;
                else if (equal(column[i], column[j]) != 1)
                    state = 5; // 1 + 4;
                break;
            case 2:                                  // [i] < [j] in normal
                if (less(column[i], column[j]) == 1) // [i] < [j] in case
                    state = 10;                      //2 + 8;
                else if (equal(column[i], column[j]) != 1)
                    state = 6; // 2 + 4;
                break;
            } // end switch
            current = (UB8)state;

            // 由于gi, gj都由原来的21位缩减到19位以及nx, cx位数的扩充，所以需要改变相应移位数
            current = current | ((UB8)i << 26);
            current = current | ((UB8)j << 45);

            pairs[ithread][index] = current;
        } // end for
    }     //end parallel

    return EXIT_SUCCESS;
}

/****************************************************************************
 * Filter dys-regulated genes, judge directions, RanComp V2.0 
 * Input: consistent pairs
 * Output: up and down-regulated genes 
 *****************************************************************************/
#define STATE_INC(i, ngv, nlv, cgv, clv)               \
    if (states[i] != NULL)                             \
    {                                                  \
        states[i]->ng += ngv;                          \
        states[i]->nl += nlv;                          \
        states[i]->cg += cgv;                          \
        states[i]->cl += clv;                          \
    }                                                  \
    else                                               \
    {                                                  \
        states[i] = malloc(sizeof(struct gene_state)); \
        states[i]->ng = ngv;                           \
        states[i]->nl = nlv;                           \
        states[i]->cg = cgv;                           \
        states[i]->cl = clv;                           \
        states[i]->state = 0;                          \
    }
#define STATE_UPDATE(iv, jv)       \
    {                              \
        if (states[j]->state == 0) \
            states[i]->iv += 1;    \
        if (states[i]->state == 0) \
            states[j]->jv += 1;    \
    }


int filter_gene_dirs_new(DATATYPE_GENESIZE n, //number of genes
                         int nthreads,
                         UB8 *pairs[nthreads],
                         DATATYPE_GENESIZE count_pairs[nthreads],
                         struct gene_state *states[n],
                         double alpha, //FDR alpha level for regulation direction
                         int max_cycles,
                         int conv_threshold,
                         float *normal_rank,
                         int *case_rank)
{
    unsigned int i, j;

    //initialization of gene states
    for (i = 0; i < nthreads; i++)
        for (j = 0; j < count_pairs[i]; j++)
        {
            if (pairs[i][j] != 0)
            {
                UB8 current = pairs[i][j];
                unsigned char state = ((unsigned char)current) & 0b1111;

                // 由于gi, gj都由原来的21位缩减到19位以及nx, cx位数的扩充，所以需要改变相应移位数以及与运算的1数目
                unsigned int c1 = ((unsigned int)(current >> 26)) & 0x7FFFF; //clear left bits
                unsigned int c2 = ((unsigned int)(current >> 45)) & 0x7FFFF;

                switch (state)
                {
                case 5: // 1 + 4, [i]>[j] in both
                    STATE_INC(c1, 0, 1, 0, 0)
                    STATE_INC(c2, 1, 0, 0, 0)
                    STATE_INC(c1, 0, 0, 0, 1)
                    STATE_INC(c2, 0, 0, 1, 0)
                    break;
                case 10: // 2 + 8, [i]<[j] in both
                    STATE_INC(c1, 1, 0, 0, 0)
                    STATE_INC(c2, 0, 1, 0, 0)
                    STATE_INC(c1, 0, 0, 1, 0)
                    STATE_INC(c2, 0, 0, 0, 1)
                    break;
                // Stable in the two groups, different direction
                case 6: // 2 + 4
                    STATE_INC(c1, 1, 0, 0, 0)
                    STATE_INC(c2, 0, 1, 0, 0)
                    STATE_INC(c1, 0, 0, 0, 1)
                    STATE_INC(c2, 0, 0, 1, 0)
                    break;
                case 9: // 1 + 8
                    STATE_INC(c1, 0, 1, 0, 0)
                    STATE_INC(c2, 1, 0, 0, 0)
                    STATE_INC(c1, 0, 0, 1, 0)
                    STATE_INC(c2, 0, 0, 0, 1)
                    break;
                case 0: //other, remove
                    break;
                } //end switch
            }     //end if
        }         //end double-for

    //screening for up or down direction using the two-tailed hypergeometric test
    int cycles = 0;
    double *pvalues;
    double p_upper;
    int gene_up = BIGNINT, gene_down = BIGNINT, gene_flat = BIGNINT;
    int gene_up_pre, gene_down_pre, gene_flat_pre;

    do
    {
        int count = 0;
//two-tailed hypergeometric test
#pragma omp parallel for reduction(+ \
                                   : count)
        for (i = 0; i < n; i++)
            if (states[i] != NULL)
            {
                states[i]->p = fisher_test((int)states[i]->ng, (int)states[i]->nl,
                                           (int)states[i]->cg, (int)states[i]->cl);
                count += 1;
            } //end if

        // copy p values to pvalues
        fprintf(stdout, "#  Cycle %d, %d genes have been tested using the Fisher exact test (aka. two-tail hypergeometric test).\n", cycles, count);

        pvalues = malloc(sizeof(double) * count);
        count = 0;
        for (i = 0; i < n; i++)
            if (states[i] != NULL)
                pvalues[count++] = states[i]->p;

        //FDR control
        p_upper = bh_threshold(count, pvalues, alpha);
        fprintf(stdout, "#  FDR = %f controlled P value threshold, %.17g. \n", alpha, p_upper);
        free(pvalues);

        //assign states
        gene_up_pre = gene_up;
        gene_down_pre = gene_down;
        gene_flat_pre = gene_flat;
        gene_up = 0;
        gene_down = 0;
        gene_flat = 0;
#pragma omp parallel for reduction(+ \
                                   : gene_up, gene_down, gene_flat)
        for (i = 0; i < n; i++)
            if (states[i] != NULL)
            {
                if ((states[i]->p) < p_upper)
                { //in case group, g/l ratio is higher, i is down-regulated.
                    if (((states[i]->ng + 1.0e-5) / (states[i]->nl + 1.0e-5)) <
                        ((states[i]->cg + 1.0e-5) / (states[i]->cl + 1.0e-5)))
                    {
                        states[i]->state = 1;
                        gene_down += 1;
                    }
                    else
                    {
                        states[i]->state = 2;
                        gene_up += 1;
                    }
                }
                else
                {
                    states[i]->state = 0;
                    gene_flat += 1;
                }
            } //end if for
        fprintf(stdout, "#  %d genes have been tested as up-regulated.\n", gene_up);
        fprintf(stdout, "#  %d genes have been tested as down-regulated.\n", gene_down);
        fprintf(stdout, "#  %d genes have no particular direction.\n", gene_flat);

//count ng, nl, cg and cl, if gene_state is 0 (flat).
//clear first
#pragma omp parallel for
        for (i = 0; i < n; i++)
            if (states[i] != NULL)
            {
                states[i]->ng = 0;
                states[i]->nl = 0;
                states[i]->cg = 0;
                states[i]->cl = 0;
            }

#pragma omp parallel
        {
            int ithread = omp_get_thread_num();
            unsigned int index, i, j;
            unsigned char state;
            UB8 current;
            for (index = 0; index < count_pairs[ithread]; index++)
            {
                current = pairs[ithread][index];
                state = ((unsigned char)current) & 0b1111;

                // 由于gi, gj都由原来的21位缩减到19位以及nx, cx位数的扩充，所以需要改变相应移位数以及与运算的1数目
                i = ((unsigned int)(current >> 26)) & 0x7FFFF; //clear left bits
                j = ((unsigned int)(current >> 45)) & 0x7FFFF;

#pragma omp critical
                switch (state)
                {
                case 5: // 1 + 4, [i]>[j]
                    STATE_UPDATE(nl, ng)
                    STATE_UPDATE(cl, cg)
                    break;
                case 10: // 2 + 8, [i]<[j]
                    STATE_UPDATE(ng, nl)
                    STATE_UPDATE(cg, cl)
                    break;
                // Stable in the two groups, different direction
                case 6: // 2 + 4 [i]<[j] --> [i]>[j]
                    //for gene i
                    if (states[j]->state == 0)
                    {
                        states[i]->ng == 0;
                        float a1i6 = 0.0, a2i6 = 0.0, ai6 = 0.0;
                        a1i6 = (float)normal_rank[i] / (float)(n / 2);
                        a2i6 = (fabs(normal_rank[j] - normal_rank[i]) + fabs(case_rank[i] - case_rank[j])) / (n / 2);
                        // ai6 = (a1i6 * a2i6);
                        ai6 = (a1i6 + a2i6) / 2;
                        states[i]->cl += ai6;
                    }
                    //for gene j
                    if (states[i]->state == 0)
                    {
                        states[j]->nl += 1;
                        float a1j6 = 0.0, a2j6 = 0.0, aj6 = 0.0;
                        a1j6 = (float)(n - normal_rank[j]) / (float)(n / 2);
                        a2j6 = (fabs(normal_rank[i] - normal_rank[j]) + fabs(case_rank[j] - case_rank[i])) / (n / 2);
                        // aj6 = (a1j6 * a2j6);
                        aj6 = (a1j6 + a2j6) / 2;
                        states[j]->cg += aj6;
                    }

                    break;
                case 9: // 1 + 8 [i]>[j] --> [i]<[j]
                    //for gene i
                    if (states[j]->state == 0)
                    {
                        states[i]->nl += 1;
                        float a1i9 = 0.0, a2i9 = 0.0, ai9 = 0.0;
                        a1i9 = (float)normal_rank[i] / (float)(n / 2);
                        a2i9 = (fabs(normal_rank[j] - normal_rank[i]) + fabs(case_rank[i] - case_rank[j])) / (n / 2);
                        // ai9 = (a1i9 * a2i9);
                        ai9 = (a1i9 + a2i9) / 2;
                        states[i]->cg += ai9;
                    }
                    //for gene j
                    if (states[i]->state == 0)
                    {
                        states[j]->ng += 1;
                        float a1j9 = 0.0, a2j9 = 0.0, aj9 = 0.0;
                        a1j9 = (float)(n - normal_rank[j]) / (float)(n / 2);
                        a2j9 = (fabs(normal_rank[i] - normal_rank[j]) + fabs(case_rank[j] - case_rank[i])) / (n / 2);
                        // aj9 = (a1j9 * a2j9);
                        aj9 = (a1j9 + a2j9) / 2;
                        states[j]->cl += aj9;
                    }

                    break;
                } // end switch
            }     // end for
        }         //end parallel
        cycles += 1;
    } while (cycles < max_cycles && (abs(gene_flat - gene_flat_pre) > conv_threshold ||
                                     abs(gene_up - gene_up_pre) > conv_threshold ||
                                     abs(gene_down - gene_down_pre) > conv_threshold));

    if (abs(gene_flat - gene_flat_pre) <= conv_threshold &&
        abs(gene_up - gene_up_pre) <= conv_threshold &&
        abs(gene_down - gene_down_pre) <= conv_threshold)
        fprintf(stdout, "#  Convergence has been reached after %d cycles.\n", cycles);
    else
        fprintf(stdout, "#  Max cycles  %d have been exceeded before the convergence.\n", max_cycles);

    return EXIT_SUCCESS;
} //end filter_gene_dirs

