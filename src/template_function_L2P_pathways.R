L2P_pathways <- function(Correlation_analysis) {

    # image: png
    suppressMessages(library(l2p))
    suppressMessages(library(magrittr))
    suppressMessages(library(dplyr))
    suppressMessages(library(tidyverse))
    suppressMessages(library(ggplot2))
    suppressMessages(library(RCurl))
    suppressMessages(library(RColorBrewer))
    
    Correlation_analysis %>% dplyr::select("Gene","CorrCoef","Pval") -> genesmat
    genes_universe = as.vector(unique(unlist(genesmat["Gene"])))
    compsum <- 0.1*length(genes_universe)
    genesmat %>% dplyr::filter(CorrCoef < 0) -> genesmat
    genesmat %>% dplyr::arrange(desc(`Pval`)) -> genesmat
    genesmat %>% dplyr::filter(!is.na(`Pval`)) -> genesmat

    sort_descending = FALSE
    if (sort_descending) {
        genes_to_include = head(genesmat["Gene"], compsum)
    } else {
        genes_to_include = tail(genesmat["Gene"], compsum)
    }
    genes_to_include <- as.vector(unique(unlist(genes_to_include)))
    
    gene_set_sources_to_include = c("GO","C2","C6","C7","H")
    categories_string <- paste(gene_set_sources_to_include, collapse=",")

    use_built_in_gene_universe = FALSE
    if (use_built_in_gene_universe) {
        x <- l2pwcats(genes_to_include, categories_string)
        print("Using built-in gene universe.")
    } else {
        x <- l2puwcats(genes_to_include, genes_universe, categories_string)
        print("Using all genes in differential expression analysis as gene universe.")
    }
    x %>%
        dplyr::arrange(pval) %>% dplyr::mutate(hitsPerc=(pwhitcount/(pwnohitcount+pwhitcount))*100) %>% dplyr::mutate(pathtotal=pwnohitcount+pwhitcount) %>% dplyr::filter(ratio >= 0) %>%
        dplyr::select("pathwayname", "category", "pathwayaccessionidentifier", "pval", "fdr", "pwhitcount", "genesinpathway", "pwnohitcount","pathtotal","hitsPerc", "inputcount", "pwuniverseminuslist","ratio")  %>% dplyr::rename(diff_ratio = ratio) %>% dplyr::filter(pval < 0.05)%>% dplyr::filter(pwhitcount >= 5) -> x
  
    print(paste0("Total number of pathways: ", nrow(x)))
    goResults <- x
    goResults %>% top_n(10, wt=-log(pval)) %>%
        dplyr::arrange(-log(pval)) -> goResults
    minp = min(goResults$pval) - 0.1*min(goResults$pval)
    maxp = max(goResults$pval) + 0.1*max(goResults$pval)
    #print(goResults$pval)
    sizemax = ceiling(max(goResults$pwhitcount)/10)*10  
    goResults %>% dplyr::mutate(pathwayname2 = stringr::str_replace_all(pathwayname, "_", " ")) -> goResults
    goResults %>% dplyr::mutate(pathwayname2 = stringr::str_wrap(pathwayname2,30)) -> goResults
    
    if (TRUE){
    
    goResults %>% dplyr::mutate(percorder = order(goResults$pval)) -> goResults
    goResults$pathwayname2 <- factor(goResults$pathwayname2, levels = goResults$pathwayname2[goResults$percorder])
    xmin = floor(min(goResults$pval))
    xmax = max(goResults$pval) 
    gplot <- goResults %>% 
               ggplot(aes(x=pval,
               y=pathwayname2, 
               colour=hitsPerc, 
               size=pwhitcount)) +
        geom_point() +
        theme(text = element_text(size=20), legend.position = "right", legend.key.height = unit(1, "cm"),
                axis.title.x = element_text(size = rel(1.2)),axis.title.y = element_text(size = rel(1.2))) +
        xlim(xmin,xmax) +
        expand_limits(colour = seq(minp, maxp, by = 10),
                size = seq(0, sizemax,by=10)) +
        labs(x="p value", y="GO term", colour="Hits (%)", size="Count") 
    print(gplot)
    }
    else{
    goResults %>% dplyr::mutate(percorder = order(goResults$hitsPerc)) -> goResults
    goResults$pathwayname2 <- factor(goResults$pathwayname2, levels = goResults$pathwayname2[goResults$percorder])
    xmin = floor(min(goResults$hitsPerc)-5)
    xmax = ceiling(max(goResults$hitsPerc)+5) 
    gplot <- goResults %>% 
               ggplot(aes(x=hitsPerc,
               y=pathwayname2, 
               colour=pval, 
               size=pwhitcount)) +
        geom_point() +
        theme_classic() +
        theme(text = element_text(size=20), legend.position = "right", legend.key.height = unit(1, "cm"),
                axis.title.x = element_text(size = rel(1.2)),axis.title.y = element_text(size = rel(1.2))) +
        xlim(xmin,xmax) +
        expand_limits(colour = seq(minp, maxp, by = 10),
                size = seq(0, sizemax,by=10)) +
        labs(x="Hits (%)", y="GO term", colour="p value", size="Count") 
    print(gplot)
    }
    return(x)
}

print("template_function_L2P_pathways.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
var_Correlation_analysis<-readRDS("var_Correlation_analysis.rds")
invisible(graphics.off())
var_L2P_pathways<-L2P_pathways(var_Correlation_analysis)
invisible(graphics.off())
saveRDS(var_L2P_pathways,"var_L2P_pathways.rds")
