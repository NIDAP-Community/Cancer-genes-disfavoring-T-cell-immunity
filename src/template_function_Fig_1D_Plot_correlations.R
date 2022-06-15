Fig_1D_Plot_correlations <- function(Correlation_analysis) {
    library(tidyverse)
    library(ggbeeswarm)
    library(ggrepel)
    
    df <- Correlation_analysis

    df <- df %>% mutate(color= case_when(Gene %in% c("MYC","CD8A","IFNG") ~ "red", TRUE~"black")) %>% arrange(color)
    highlight_df <- df %>% filter(color=="red")

    g <- ggplot(df,aes(x=1,y=CorrCoef,color=color)) +
         geom_quasirandom(method='pseudorandom',bandwidth=0.001) +
         theme_bw() +
         labs(x=expression(CD3E^HI~melanoma), y=expression(atop("Correlation Coefficient",paste("Correlation to CYT")))) +
         scale_color_manual("color", values = c("black" = "black", "red" = "red")) +
         geom_text_repel(data=subset(df, Gene %in% c("MYC")), aes(label=Gene), nudge_x=0.05, nudge_y = -0.10) +
         geom_text_repel(data=subset(df, Gene %in% c("CD8A","IFNG")), aes(label=Gene), nudge_x=0.05, nudge_y = 0.05) +
         theme(plot.margin = margin(0.5,7,0.5,7, "cm"), legend.position="none")
    print(g)
    
    return(highlight_df)
}

print("template_function_Fig_1D_Plot_correlations.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
var_Correlation_analysis<-readRDS("var_Correlation_analysis.rds")
invisible(graphics.off())
var_Fig_1D_Plot_correlations<-Fig_1D_Plot_correlations(var_Correlation_analysis)
invisible(graphics.off())
saveRDS(var_Fig_1D_Plot_correlations,"var_Fig_1D_Plot_correlations.rds")
