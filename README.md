# InDEA
This program is the implement of the Individual Differential Expression Analysis (InDEA) algorithm, which is an approach to detect individual differential expressed genes (DEGs). The algorithm was designed on the basis of RankCompV2 considering the assumption of gene’s rank position and relative rank difference would influence this gene pair’s contribution in DEG identification. The RankCompV2 algorithm can be found in https://github.com/pathint/reoa. 

# Algorithm
For a gene pair (A, G) of a given gene A within a sample, two REO pattern can be founded, the expression of gene A higher or lower than that of gene G (A > G or A< G).The analysis of InDEA start with determination of the background list. Gene pairs with coincident relative expression orderings (REOs) pattern in more than 95% of all normal samples were included in the background list. Then reversal gene pairs within each individual disease sample were determined,compared with those in the background list.  After that, considering that gene’s rank position and relative rank difference would influence gene pair’s contribution in DEGs identification, for a given gene A, the contribution of reversal gene pairs involved A was calculated. Fisher’s exact test was used to identify whether the gene was dysregulated. To minimize the confound effects of the dysregulation of the partner genes, a filter process by iteratively excluding gene pairs with DEGs from background gene pairs and reversal gene pairs was carried out. Then Fisher’s exact test was applied again to identify new DEGs. The filter process was performed until the DEGs list tend to be stable.

# Install
There are two possible ways to install the program.

1.Install from the Precompiled Execulables. Just copy the binary file (InDEA) to the folder where the executables are located, such as `/usr/local/bin`. Make sure the files have the executable permission. If it does not, use `chmod 755 InDEA` to make the modification. 

2.Compile from the Source Files. Under the src folder, use `gcc` tool to compile the sources. Notably, `-lm` and `-fopenmp` option is needed while using `gcc` command.
e.g. `gcc main.c indea.c statistics.c -fopenmp -lm -o InDEA`

# Usage
Command format: `InDEA n_files n_genes {data_file...} {sample_size...} {fdr_level or exception_number}`

n_files: the number of input data files

n_genes: the number of genes in each data file, aka. total row number.

{data_file...}: list of data file names. All files should have the same number of rows. 

{sample_size...}: list of sample sizes (aka. total column number) of the corresponding data files

{fdr_level or exception_number}: False discovery rate (FDR) level or exception number, which are used to select significantly stable pairs Float numbers between 0. and 1. are taken as FDR levels, stable gene pairs would discovered in Binomial mode; otherwise, they are taken as exception numbers while stable gene pairs would discovered in highly stable mode.

The Default FDR level for identifing DEGs is 0.05

e.g. 
```
InDEA 2 20572 normal.txt disease.txt 207 25 0.05 
```
```
InDEA 2 20572 normal.txt disease.txt 207 25 10 
```

The input data sets should be given as text-based data matrix data files. One file contains the expression matrix of one group of samples. The number of rows corresponds to the total number of genes and the number of columns corresponds to the sample size. The values should be tab or space delimited. It should be noticed that the number of genes within normal sample data and disease sample data should be the same and sequentially consistent.

The output of dysregulated genes is given by the indices (starting from 0). The output matrix (.dat file) contains the dysregulate status of every gene in every disease samples. There are three possible states for each gene: 0, non-dysregulated genes; 1, down-regulated genes; 2, up-regulated genes.


