broad_skcm_normalized_genelevel <- function(SKCM_rnaseqv2_illuminahiseq_rnaseqv2_unc_edu_Level_3_RSEM_genes_normalized_good) {
    
    dg <- SKCM_rnaseqv2_illuminahiseq_rnaseqv2_unc_edu_Level_3_RSEM_genes_normalized_good
    
    dg$Gene <- unlist(lapply(strsplit(dg$gene_id,split='\\|'),function(x) x[[1]]))
    cc <- gsub("\\.","_", colnames(dg))
    print(cc)
    colnames(dg) <- cc

    return(dg)
}

print("template_function_broad_skcm_normalized_genelevel.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
var_SKCM_rnaseqv2_illuminahiseq_rnaseqv2_unc_edu_Level_3_RSEM_genes_normalized_good<-readRDS("var_SKCM_rnaseqv2_illuminahiseq_rnaseqv2_unc_edu_Level_3_RSEM_genes_normalized_good.rds")
invisible(graphics.off())
var_broad_skcm_normalized_genelevel<-broad_skcm_normalized_genelevel(var_SKCM_rnaseqv2_illuminahiseq_rnaseqv2_unc_edu_Level_3_RSEM_genes_normalized_good)
invisible(graphics.off())
saveRDS(var_broad_skcm_normalized_genelevel,"var_broad_skcm_normalized_genelevel.rds")
