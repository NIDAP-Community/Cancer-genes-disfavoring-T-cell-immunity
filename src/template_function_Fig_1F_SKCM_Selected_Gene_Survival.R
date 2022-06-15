Fig_1F_SKCM_Selected_Gene_Survival <- function(  TCGA_SKCM_harmonized_clinical_data,  Normalized_counts_SKCM, Correlation_analysis) {
   
    library(matrixStats)
    library(survival)
    library(survminer)

    df <- Normalized_counts_SKCM
    dhc <- TCGA_SKCM_harmonized_clinical_data
    dg <- Correlation_analysis

    # select only tumor samples
    ind_tp_grp <- grep("TM",dhc$shortLetterCode)
    pt_samples <- dhc$Sample[ind_tp_grp]
    dfpt <- df[,c("Gene",pt_samples)]

    # select the expression data for the gene list of choice
    gn <- dg %>% filter(CorrCoef < 0, Pval < 0.05) 
    gn <- dg %>% filter(CorrCoef < 0, Pval < 0.05) %>% dplyr::pull(Gene)
    ind <- match(gn, dfpt$Gene)  
    dfg <- as.matrix(dfpt[ind,2:dim(dfpt)[2]])
    dfg <- t(scale(t(dfg)))

    # classify samples above and below upper and lower quartiles of aggregated gene expression differences 
    md_thr <- rowMedians(dfg) # median across samples for the list of genes of choice
    diff_thr <- dfg-md_thr # difference between the expression and median threshold
    above_threshold <- colMedians(diff_thr)>0 # samples with difference above threshold
    below_threshold <- colMedians(diff_thr)<0 # samples with difference below threshold 
    label_above_thr <- "above median"
    label_below_thr <- "below median"

    samples_above_threshold <- colnames(dfpt[,-1])[above_threshold]
    samples_below_threshold <- colnames(dfpt[,-1])[below_threshold]

    ind_patients_above_threshold_hclinial <- match(samples_above_threshold,dhc$Sample)
    ind_patients_below_threshold_hclinial <- match(samples_below_threshold,dhc$Sample)

    # prepare the vectors necessary for survival analysis 
    surv_days_all <- ifelse(dhc$"vital_status_cl3p" == "Dead", dhc$"days_to_death_cl3p", dhc$"days_to_last_follow_up_cl3p")
    surv_days <- surv_days_all[c(ind_patients_above_threshold_hclinial, ind_patients_below_threshold_hclinial)]
        cs <- dhc$"vital_status_cl3p"[c(ind_patients_above_threshold_hclinial, ind_patients_below_threshold_hclinial)]
        csn <- as.numeric(cs=="Dead") # if the vital_status is "Dead" attribute 1; else, attribute 0
    surv_censor <- csn
    cluster <- c(rep(0,length(ind_patients_above_threshold_hclinial)), rep(1, length(ind_patients_below_threshold_hclinial)))

    # prepare the dataframe necessary for survival analysis
    ds <- as.data.frame(matrix(NA, ncol=3, nrow=length(cluster)))
    ds[,1] <- surv_days
    ds[,2] <- surv_censor
    ds[,3] <- cluster
    colnames(ds) <- c("surv_days", "surv_censor", "cluster")

    # evaluate survival curves
    fit1 <- survfit(Surv(as.numeric(surv_days),surv_censor) ~ cluster, data = ds) # single variable

    # explicit measures of survival fit
    pval <- surv_pvalue(fit1,ds)
    sfit <- summary(fit1)$table
    print(pval)
    print(sfit)

    p <- survminer::ggsurvplot(fit1, #plot just the KM curves 
            data = ds,
            pval = TRUE, # p-value of the Log-Rank test comparing the groups 
            conf.int = TRUE,  # 95% confidence interval
            risk.table = TRUE, # display the risk table below the graph 
            risk.table.col = "strata",  # control what is displayed in the table
            title = unique(dhc$"project_id"),
            surv.median.line = "hv", # display horizontal and vertical line for median survival
            ggtheme = theme_gray(), # have a background theme
            palette = c("red","black"),  # control the colors for the KM curves and associated data
            ncensor.plot = TRUE, # add censor plot below the risk table to display the censor data
            legend.labs = c(label_above_thr, label_below_thr)
            ) 
    print(p)           

    df_out <- as.data.frame(sfit)
    df_out <- cbind(as.data.frame(rownames(df_out)), df_out)
    colnames(df_out)[c(1,6:7)] <- c("cluster","rmean", "se_rmean")

    select_output <- 1
    if (select_output==1){
        df_out_r <- pval
    } else if (select_output==2){
        df_out_r <- df_out
    } else if (select_output==3){
            df_out_r <- ds
    } 

    return(df_out_r)
}

print("template_function_Fig_1F_SKCM_Selected_Gene_Survival.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
var_TCGA_SKCM_harmonized_clinical_data<-readRDS("var_TCGA_SKCM_harmonized_clinical_data.rds")
var_Normalized_counts_SKCM<-readRDS("var_Normalized_counts_SKCM.rds")
var_Correlation_analysis<-readRDS("var_Correlation_analysis.rds")
invisible(graphics.off())
var_Fig_1F_SKCM_Selected_Gene_Survival<-Fig_1F_SKCM_Selected_Gene_Survival(var_TCGA_SKCM_harmonized_clinical_data,var_Normalized_counts_SKCM,var_Correlation_analysis)
invisible(graphics.off())
saveRDS(var_Fig_1F_SKCM_Selected_Gene_Survival,"var_Fig_1F_SKCM_Selected_Gene_Survival.rds")
