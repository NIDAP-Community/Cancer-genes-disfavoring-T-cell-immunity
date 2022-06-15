Fig_1E_Pathways_selected <- function(L2P_pathways) {

    df <-  L2P_pathways

    df %>% mutate(pathwayname = toupper(pathwayname)) -> df1
    
    pathway_select <- c("KRIGE_RESPONSE_TO_TOSEDOSTAT_24HR_DN","RIBOSOME BIOGENESIS","MYC_UP.V1_UP","SMID_BREAST_CANCER_LUMINAL_A_UP","BILD_MYC_ONCOGENIC_SIGNATURE","GTP METABOLIC PROCESS","REACTOME_EUKARYOTIC_TRANSLATION_INITIATION","ORGANONITROGEN COMPOUND METABOLIC PROCESS","CHARAFE_BREAST_CANCER_LUMINAL_VS_BASAL_UP")

    df1 %>% filter(pathwayname %in% pathway_select) -> df1
    df1 <- unique(df1)
    df1 %>% arrange(desc(fdr)) -> df1

    path_rename <- list(
        "KRIGE_RESPONSE_TO_TOSEDOSTAT_24HR_DN" = "Response to Tosedostat (24hr)",
        "ORGANONITROGEN COMPOUND METABOLIC PROCESS" = "Organonitrogen Metabolic Process",
        "RIBOSOME BIOGENESIS"= "Ribosome Biogenesis",
        "MYC_UP.V1_UP" = "MYC Up",
        "CHARAFE_BREAST_CANCER_LUMINAL_VS_BASAL_UP" = "Luminal Breast Cancer",
        "BILD_MYC_ONCOGENIC_SIGNATURE" = "MYC Oncogenic Signature",
        "REACTOME_EUKARYOTIC_TRANSLATION_INITIATION" = "Translation Initiation",
        "GTP METABOLIC PROCESS" = "GTP Hydrolysis"
    )
    df1 %>% mutate("pathway" = unlist(path_rename[pathwayname])) -> df1
    df1 %>% mutate(fdr = -log(fdr)) %>% arrange(desc(fdr))-> df1

    g <- ggplot(df1, aes(x=reorder(pathway, fdr), y=fdr)) + 
            geom_bar(stat = "identity") +
            theme_bw() + 
            theme(axis.title.y = element_blank()) +
            labs(y=expression(paste("-lo",g[10],"(FDR adjusted ",italic("p")," value)"))) +
            coord_flip() 
    
    print(g)
    return(df)
        
}

print("template_function_Fig_1E_Pathways_selected.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
var_L2P_pathways<-readRDS("var_L2P_pathways.rds")
invisible(graphics.off())
var_Fig_1E_Pathways_selected<-Fig_1E_Pathways_selected(var_L2P_pathways)
invisible(graphics.off())
saveRDS(var_Fig_1E_Pathways_selected,"var_Fig_1E_Pathways_selected.rds")
