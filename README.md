# Cancer-genes-disfavoring-T-cell-immunity
Code accompanying publication "Cancer genes disfavoring T cell immunity identified via integrated systems approach"

TCGA Data Analysis

Data analysis and visualization were performed in the NIH Integrated Analysis Portal (NIDAP) using R programs developed on Foundry (Palantir Technologies). Normalized RNA-Seq gene expression data and merged clinical metadata for TCGA Skin Cutaneous Melanoma (SKCM) was downloaded from the FireBrowse portal (https://gdac.broadinstitute.org/) and subsequently log-normalized using limma-voom (1).  Further downstream analyses were limited to metastatic tumor samples having high (above median) CD3E gene expression.  Correlation analysis was performed on all genes to the combined geometric mean of PRF1 and GZMA expression (CYT).  Pathway enrichment of genes with significantly negative Pearson's Correlation Coefficients (p<0.05) was performed using l2p (https://github.com/ccbr/l2p) on GO, KEGG, and various MSigDB pathway databases, including Hallmark (2), C2 and C7 collections.  Survival analysis of samples with high vs. low expression of significantly negatively correlated genes was performed using survminer (3).  

References:

1.	Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth, G.K. (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47.
2.	Liberzon A, Birger C, Thorvaldsdóttir H, Ghandi M, Mesirov JP, Tamayo P. The Molecular Signatures Database (MSigDB) hallmark gene set collection. Cell Syst. 2015 Dec 23;1(6):417-425.
3.	Alboukadel Kassambara, Marcin Kosinski and Przemyslaw Biecek (2019). survminer: Drawing Survival Curves using 'ggplot2'. R package version 0.4.6. https://CRAN.R-project.org/package=survminer

