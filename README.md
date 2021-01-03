# InDEA
This program is the implement of the Individual Differential Expression Analysis (InDEA) algorithm, which is an approach to detect individual dysregulated genes. InDEA can be used in gene expression, methylation and other gene data sets.

# Algorithm
The analysis of InDEA start with Determination of the background list. Gene pairs with REO pattern stable in more than 95% of all normal samples was included in the background list. Then reversal gene pairs within individual disease sample were determined. Gene pairs would be determined as reversal pairs while their REO pattern reversal in disease sample, compared with those in the background list. After that, considering gene’s rank position and relative rank difference would influence gene pair’s contribution in DEG identification, for a given gene A, the contribution of reversal gene pairs involved A was calculated. Fisher’s exact test was used to identify whether the gene was dysregulated in the disease sample by testing the null hypothesis that the proportion of upregulation-supporting reversal gene pairs of a gene was equal to the proportion of downregulation-supporting reversal gene pairs of the gene. To minimize the confound effects of the dysregulation of the partner genes, a filter process by iteratively excluding gene pairs with DEGs from stable gene pairs and reversal gene pairs was carried out. Then Fisher’s exact test was applied again to identify new DEGs. The filter process was performed until the DEG list tend to be stable.

# Install
There are two possible ways to install the program.

1.Install from the Precompiled Execulables. Just copy the binary file to the folder where the executables are located, such as `/usr/local/bin`. Make sure the files have the executable permission. If it does not, use `chmod 755 InDEA` to make the modification. 

2.Compile from the Source Files. Under the src folder, use `gcc` tool to compile the sources. Notably, `-lm` and `-fopenmp` option is needed while using `gcc` command.
e.g. `gcc main.c indea.c statistics.c -fopenmp -lm -o InDEA`

# Usage
Command format: `InDEA n_files n_genes {data_file...} {sample_size...} {fdr_level or exception_number}`

n_files: the number of input data files

n_genes: the number of genes in each data file, aka. total row number.

{data_file...}: list of data file names. All files should have the same number of rows.

{sample_size...}: list of sample sizes (aka. total column number) of the corresponding data files

{fdr_level or exception_number}: FDR level or exception number, which are used to select significantly stable pairs Float numbers between 0. and 1. are taken as FDR levels, stable gene pairs would discovered in Binomial mode; otherwise, they are taken as exception numbers while stable gene pairs would discovered in highly stable mode.

e.g. 
```
InDEA 2 20572 normal.dat disease.dat 207 25 0.05
```
```
InDEA 2 20572 normal.dat disease.dat 207 25 10
```

The input data sets should be given as text-based data matrix data files. One file contains the microarray value matrix of one group of samples. The number of rows corresponds to the total number of gene probes and the number of columns corresponds to the sample size. The values should be tab or space delimited. It should be noticed that the number of genes within normal sample data and disease sample data should be the same.
The output of dysregulated genes is given by the indices (starting from 0). The output matrix (.dat file) contains the dysregulate status of every gene in every disease samples. There are three possible states for each gene: 0, non-dysregulated genes; 1, down-regulated genes; 2, up-regulated genes.


