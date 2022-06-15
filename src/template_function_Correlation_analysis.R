Correlation_analysis <- function(Samples_high_CD3E_Expression) {
    
    df <- Samples_high_CD3E_Expression
    
    #evaluate geometric or arithmentic mean
    gl <- trimws(unlist(strsplit(c("GZMA,PRF1"), ",")), which=c("both")) # unpack the gene list provided by the user and remove white spaces    
    ind_gn <- match(gl, df$Gene)
    dfgn <- df[ind_gn,2:dim(df)[2]]

    dx <- colSums(dfgn)/2 #Calculate geometric mean with input log values 
    dg <- df[-ind_gn,]
    dy <- df[-ind_gn,2:dim(df)[2]]

    # evaluate correlation
    df_out <- as.data.frame(as.matrix(NA, nrow=dim(dg)[1], ncol=5))
    
    for (i in 1:dim(dg)[1]){
        cc <- cor.test(dx, as.numeric(dy[i,]), method="pearson")
        df_out[i,1] <- dg$Gene[i]
        df_out[i,2] <- round(cc$estimate,3)
        df_out[i,3] <- cc$p.value
        df_out[i,4] <- round(cc$conf.int[1],3)
        df_out[i,5] <- round(cc$conf.int[2],3)
    }
    colnames(df_out) <- c("Gene", "CorrCoef", "Pval","LowerCI", "UpperCI")

    return(df_out)

}

print("template_function_Correlation_analysis.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
var_Samples_high_CD3E_Expression<-readRDS("var_Samples_high_CD3E_Expression.rds")
invisible(graphics.off())
var_Correlation_analysis<-Correlation_analysis(var_Samples_high_CD3E_Expression)
invisible(graphics.off())
saveRDS(var_Correlation_analysis,"var_Correlation_analysis.rds")
