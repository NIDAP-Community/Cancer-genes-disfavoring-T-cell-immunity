Samples_high_CD3E_Expression <- function(TCGA_SKCM_harmonized_clinical_data, Normalized_counts_SKCM) {
    
    df <- Normalized_counts_SKCM
    dhc <- TCGA_SKCM_harmonized_clinical_data

    # select only tumor samples
    ind_tp_grp <- grep("TM",dhc$shortLetterCode)
    pt_samples <- dhc$Sample[ind_tp_grp]
    dfpt <- df[,c("Gene",pt_samples)] # raw counts for tumor samples
 
    # select the expression data for the gene of choice
    gn <- c("CD3E")
    ind <- match(gn, dfpt$Gene)  
    dfg <- as.numeric(dfpt[ind,2:dim(dfpt)[2]])

    # classify samples above and below threshold
    above_threshold <- dfg >= median(dfg)
    below_threshold <- dfg < median(dfg)
    label_above_thr <- "median"
    label_below_thr <- "median"
    
    samples_above_threshold <- colnames(dfpt[,-1])[above_threshold]
    samples_below_threshold <- colnames(dfpt[,-1])[below_threshold]

    ind_patients_above_threshold_hclinial <- match(samples_above_threshold,dhc$Sample)
    ind_patients_below_threshold_hclinial <- match(samples_below_threshold,dhc$Sample)

    # select output dataset (1 for above threshold, 0 below threshold)
    select_output_threshold <- TRUE

    # select output dataset (1 for raw counts, 0 for clinical data)
    select_output_data_type <- TRUE

    if (select_output_threshold) {
        if (select_output_data_type) {
            df_out <- dfpt[,c("Gene", samples_above_threshold)]
        } else {
            df_out <- dhc[ind_patients_above_threshold_hclinial,]
        }
    } else {
    if (select_output_data_type) {
            df_out <- dfpt[,c("Gene", samples_below_threshold)]
        } else {
            df_out <- dhc[ind_patients_below_threshold_hclinial,]
        } 
    }

    return(df_out)
}

print("template_function_Samples_high_CD3E_Expression.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
var_TCGA_SKCM_harmonized_clinical_data<-readRDS("var_TCGA_SKCM_harmonized_clinical_data.rds")
var_Normalized_counts_SKCM<-readRDS("var_Normalized_counts_SKCM.rds")
invisible(graphics.off())
var_Samples_high_CD3E_Expression<-Samples_high_CD3E_Expression(var_TCGA_SKCM_harmonized_clinical_data,var_Normalized_counts_SKCM)
invisible(graphics.off())
saveRDS(var_Samples_high_CD3E_Expression,"var_Samples_high_CD3E_Expression.rds")
