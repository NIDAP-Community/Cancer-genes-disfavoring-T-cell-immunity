Normalized_counts_SKCM <- function(Filtered_Counts_SKCM, TCGA_SKCM_harmonized_clinical_data) {
    
    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(limma)
    library(tidyverse)
    library(edgeR)
    library(ggplot2)
    library(plotly)
    library(dplyr)
    library(RColorBrewer)
    library(colorspace)
    library(stringr)
    library(RCurl)
    library(reshape2)
    library(gridExtra)
    library(amap)
    library(lattice)
    library(gplots)
    library(gridGraphics)
    library(dendsort)

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##
    
    counts_matrix <- Filtered_Counts_SKCM
    sample_metadata <- TCGA_SKCM_harmonized_clinical_data
    samples_to_include = c("TCGA_3N_A9WB_06A_11R_A38C_07","TCGA_3N_A9WC_06A_11R_A38C_07","TCGA_3N_A9WD_06A_11R_A38C_07","TCGA_BF_A1PU_01A_11R_A18S_07","TCGA_BF_A1PV_01A_11R_A18U_07","TCGA_BF_A1PX_01A_12R_A18T_07","TCGA_BF_A1PZ_01A_11R_A18S_07","TCGA_BF_A1Q0_01A_21R_A18S_07","TCGA_BF_A3DJ_01A_11R_A20F_07","TCGA_BF_A3DL_01A_11R_A20F_07","TCGA_BF_A3DM_01A_11R_A20F_07","TCGA_BF_A3DN_01A_11R_A20F_07","TCGA_BF_A5EO_01A_12R_A27Q_07","TCGA_BF_A5EP_01A_12R_A27Q_07","TCGA_BF_A5EQ_01A_21R_A27Q_07","TCGA_BF_A5ER_01A_12R_A27Q_07","TCGA_BF_A5ES_01A_11R_A27Q_07","TCGA_BF_A9VF_01A_11R_A37K_07","TCGA_BF_AAOU_01A_12R_A39D_07","TCGA_BF_AAOX_01A_11R_A39D_07","TCGA_BF_AAP0_06A_11R_A39D_07","TCGA_BF_AAP1_01A_11R_A39D_07","TCGA_BF_AAP2_01A_11R_A40A_07","TCGA_BF_AAP4_01A_11R_A40A_07","TCGA_BF_AAP6_01A_11R_A40A_07","TCGA_BF_AAP7_01A_11R_A40A_07","TCGA_BF_AAP8_01A_11R_A40A_07","TCGA_D3_A1Q1_06A_21R_A18T_07","TCGA_D3_A1Q3_06A_11R_A18T_07","TCGA_D3_A1Q4_06A_11R_A18T_07","TCGA_D3_A1Q5_06A_11R_A18T_07","TCGA_D3_A1Q6_06A_11R_A18T_07","TCGA_D3_A1Q7_06A_11R_A18T_07","TCGA_D3_A1Q8_06A_11R_A18T_07","TCGA_D3_A1Q9_06A_11R_A18T_07","TCGA_D3_A1QA_07A_11R_A37K_07","TCGA_D3_A1QA_06A_11R_A18T_07","TCGA_D3_A1QB_06A_11R_A18T_07","TCGA_D3_A2J6_06A_11R_A18T_07","TCGA_D3_A2J7_06A_11R_A18T_07","TCGA_D3_A2J8_06A_11R_A18T_07","TCGA_D3_A2J9_06A_11R_A18T_07","TCGA_D3_A2JA_06A_11R_A18T_07","TCGA_D3_A2JB_06A_11R_A18T_07","TCGA_D3_A2JC_06A_11R_A18T_07","TCGA_D3_A2JD_06A_11R_A18T_07","TCGA_D3_A2JE_06A_11R_A37K_07","TCGA_D3_A2JF_06A_11R_A18S_07","TCGA_D3_A2JG_06A_11R_A18T_07","TCGA_D3_A2JH_06A_11R_A18T_07","TCGA_D3_A2JK_06A_11R_A18S_07","TCGA_D3_A2JL_06A_11R_A18S_07","TCGA_D3_A2JN_06A_11R_A18S_07","TCGA_D3_A2JO_06A_11R_A18S_07","TCGA_D3_A2JP_06A_11R_A18S_07","TCGA_D3_A3BZ_06A_12R_A18S_07","TCGA_D3_A3C1_06A_12R_A18S_07","TCGA_D3_A3C3_06A_12R_A18S_07","TCGA_D3_A3C6_06A_12R_A18U_07","TCGA_D3_A3C7_06A_11R_A18U_07","TCGA_D3_A3C8_06A_12R_A18S_07","TCGA_D3_A3CB_06A_11R_A18S_07","TCGA_D3_A3CC_06A_11R_A18S_07","TCGA_D3_A3CE_06A_11R_A18S_07","TCGA_D3_A3CF_06A_11R_A18T_07","TCGA_D3_A3ML_06A_11R_A21D_07","TCGA_D3_A3MO_06A_11R_A21D_07","TCGA_D3_A3MR_06A_11R_A21D_07","TCGA_D3_A3MU_06A_11R_A21D_07","TCGA_D3_A3MV_06A_11R_A21D_07","TCGA_D3_A51E_06A_11R_A266_07","TCGA_D3_A51F_06A_11R_A266_07","TCGA_D3_A51G_06A_11R_A266_07","TCGA_D3_A51H_06A_12R_A266_07","TCGA_D3_A51J_06A_11R_A266_07","TCGA_D3_A51K_06A_11R_A266_07","TCGA_D3_A51N_06A_11R_A266_07","TCGA_D3_A51R_06A_11R_A266_07","TCGA_D3_A51T_06A_11R_A266_07","TCGA_D3_A5GL_06A_11R_A27Q_07","TCGA_D3_A5GN_06A_11R_A27Q_07","TCGA_D3_A5GO_06A_12R_A27Q_07","TCGA_D3_A5GR_06A_11R_A27Q_07","TCGA_D3_A5GS_06A_11R_A27Q_07","TCGA_D3_A5GT_01A_12R_A311_07","TCGA_D3_A5GU_06A_11R_A27Q_07","TCGA_D3_A8GB_06A_11R_A37K_07","TCGA_D3_A8GC_06A_11R_A37K_07","TCGA_D3_A8GD_06A_11R_A37K_07","TCGA_D3_A8GE_06A_11R_A37K_07","TCGA_D3_A8GI_06A_11R_A37K_07","TCGA_D3_A8GJ_06A_11R_A37K_07","TCGA_D3_A8GK_06A_11R_A37K_07","TCGA_D3_A8GL_06A_11R_A37K_07","TCGA_D3_A8GM_06A_11R_A37K_07","TCGA_D3_A8GN_06A_11R_A37K_07","TCGA_D3_A8GO_06A_11R_A37K_07","TCGA_D3_A8GP_06A_11R_A37K_07","TCGA_D3_A8GQ_06A_11R_A37K_07","TCGA_D3_A8GR_06A_11R_A37K_07","TCGA_D3_A8GS_06A_12R_A37K_07","TCGA_D3_A8GV_06A_11R_A37K_07","TCGA_D9_A148_06A_11R_A18S_07","TCGA_D9_A149_06A_11R_A18S_07","TCGA_D9_A1JW_06A_11R_A18S_07","TCGA_D9_A1JX_06A_11R_A18S_07","TCGA_D9_A1X3_06A_11R_A18S_07","TCGA_D9_A3Z1_06A_11R_A239_07","TCGA_D9_A3Z3_06A_11R_A239_07","TCGA_D9_A3Z4_01A_11R_A239_07","TCGA_D9_A4Z2_01A_11R_A24X_07","TCGA_D9_A4Z3_01A_11R_A266_07","TCGA_D9_A4Z5_01A_11R_A266_07","TCGA_D9_A4Z6_06A_12R_A266_07","TCGA_D9_A6E9_06A_12R_A311_07","TCGA_D9_A6EA_06A_11R_A311_07","TCGA_D9_A6EC_06A_11R_A311_07","TCGA_D9_A6EG_06A_12R_A32P_07","TCGA_DA_A1HV_06A_21R_A18S_07","TCGA_DA_A1HW_06A_11R_A18U_07","TCGA_DA_A1HY_06A_11R_A18T_07","TCGA_DA_A1I0_06A_11R_A20F_07","TCGA_DA_A1I1_06A_12R_A18U_07","TCGA_DA_A1I2_06A_21R_A18U_07","TCGA_DA_A1I4_06A_11R_A18U_07","TCGA_DA_A1I5_06A_11R_A18T_07","TCGA_DA_A1I7_06A_22R_A18S_07","TCGA_DA_A1I8_06A_11R_A18T_07","TCGA_DA_A1IA_06A_11R_A18S_07","TCGA_DA_A1IB_06A_11R_A18S_07","TCGA_DA_A1IC_06A_11R_A18S_07","TCGA_DA_A3F2_06A_11R_A20F_07","TCGA_DA_A3F3_06A_11R_A20F_07","TCGA_DA_A3F5_06A_11R_A20F_07","TCGA_DA_A3F8_06A_11R_A20F_07","TCGA_DA_A95V_06A_11R_A37K_07","TCGA_DA_A95W_06A_11R_A37K_07","TCGA_DA_A95X_06A_11R_A37K_07","TCGA_DA_A95Y_06A_11R_A37K_07","TCGA_DA_A95Z_06A_11R_A37K_07","TCGA_DA_A960_01A_11R_A37K_07","TCGA_EB_A1NK_01A_11R_A18T_07","TCGA_EB_A24C_01A_11R_A18T_07","TCGA_EB_A24D_01A_11R_A18T_07","TCGA_EB_A299_01A_21R_A18U_07","TCGA_EB_A3HV_01A_11R_A21D_07","TCGA_EB_A3XB_01A_11R_A239_07","TCGA_EB_A3XC_01A_11R_A239_07","TCGA_EB_A3XD_01A_22R_A239_07","TCGA_EB_A3XE_01A_12R_A239_07","TCGA_EB_A3XF_01A_31R_A239_07","TCGA_EB_A3Y6_01A_21R_A239_07","TCGA_EB_A3Y7_01A_11R_A239_07","TCGA_EB_A41A_01A_11R_A24X_07","TCGA_EB_A41B_01A_11R_A24X_07","TCGA_EB_A42Y_01A_12R_A24X_07","TCGA_EB_A42Z_01A_12R_A24X_07","TCGA_EB_A430_01A_11R_A24X_07","TCGA_EB_A431_01A_11R_A266_07","TCGA_EB_A44N_01A_11R_A266_07","TCGA_EB_A44O_01A_11R_A266_07","TCGA_EB_A44P_01A_11R_A266_07","TCGA_EB_A44Q_06A_11R_A266_07","TCGA_EB_A44R_06A_41R_A266_07","TCGA_EB_A4IQ_01A_12R_A266_07","TCGA_EB_A4IS_01A_21R_A266_07","TCGA_EB_A4OY_01A_11R_A266_07","TCGA_EB_A4OZ_01A_12R_A266_07","TCGA_EB_A4P0_01A_41R_A266_07","TCGA_EB_A4XL_01A_11R_A27Q_07","TCGA_EB_A51B_01A_11R_A27Q_07","TCGA_EB_A550_01A_61R_A27Q_07","TCGA_EB_A551_01A_21R_A27Q_07","TCGA_EB_A553_01A_12R_A27Q_07","TCGA_EB_A57M_01A_51R_A311_07","TCGA_EB_A5FP_01A_11R_A27Q_07","TCGA_EB_A5KH_06A_11R_A27Q_07","TCGA_EB_A5SE_01A_11R_A311_07","TCGA_EB_A5SF_01A_11R_A311_07","TCGA_EB_A5SG_06A_11R_A311_07","TCGA_EB_A5SH_06A_11R_A311_07","TCGA_EB_A5UL_06A_11R_A311_07","TCGA_EB_A5UM_01A_11R_A311_07","TCGA_EB_A5UN_06A_11R_A311_07","TCGA_EB_A5VU_01A_21R_A32P_07","TCGA_EB_A5VV_06A_11R_A32P_07","TCGA_EB_A6L9_06A_11R_A32P_07","TCGA_EB_A6QY_01A_12R_A32P_07","TCGA_EB_A6QZ_01A_12R_A32P_07","TCGA_EB_A6R0_01A_12R_A32P_07","TCGA_EB_A82B_01A_11R_A352_07","TCGA_EB_A82C_01A_11R_A352_07","TCGA_EB_A85I_01A_11R_A352_07","TCGA_EB_A85J_01A_12R_A352_07","TCGA_EB_A97M_01A_11R_A38C_07","TCGA_EE_A17X_06A_11R_A18S_07","TCGA_EE_A17Y_06A_11R_A18T_07","TCGA_EE_A17Z_06A_11R_A18S_07","TCGA_EE_A180_06A_11R_A21D_07","TCGA_EE_A181_06A_11R_A18S_07","TCGA_EE_A182_06A_11R_A18T_07","TCGA_EE_A183_06A_11R_A18S_07","TCGA_EE_A184_06A_11R_A18S_07","TCGA_EE_A185_06A_11R_A18S_07","TCGA_EE_A20B_06A_11R_A18U_07","TCGA_EE_A20C_06A_11R_A18S_07","TCGA_EE_A20F_06A_21R_A18S_07","TCGA_EE_A20H_06A_11R_A18S_07","TCGA_EE_A20I_06A_11R_A18U_07","TCGA_EE_A29A_06A_12R_A18U_07","TCGA_EE_A29B_06A_11R_A18U_07","TCGA_EE_A29C_06A_21R_A18S_07","TCGA_EE_A29D_06A_11R_A18T_07","TCGA_EE_A29E_06A_11R_A18T_07","TCGA_EE_A29G_06A_12R_A18T_07","TCGA_EE_A29H_06A_12R_A18S_07","TCGA_EE_A29L_06A_12R_A18S_07","TCGA_EE_A29M_06A_11R_A18T_07","TCGA_EE_A29N_06A_12R_A18S_07","TCGA_EE_A29P_06A_11R_A18T_07","TCGA_EE_A29Q_06A_11R_A18T_07","TCGA_EE_A29R_06A_11R_A18T_07","TCGA_EE_A29S_06A_11R_A18T_07","TCGA_EE_A29T_06A_11R_A18T_07","TCGA_EE_A29V_06A_12R_A18S_07","TCGA_EE_A29W_06A_11R_A18U_07","TCGA_EE_A29X_06A_11R_A18T_07","TCGA_EE_A2A0_06A_11R_A18T_07","TCGA_EE_A2A1_06A_11R_A18T_07","TCGA_EE_A2A2_06A_11R_A18T_07","TCGA_EE_A2A5_06A_11R_A18T_07","TCGA_EE_A2A6_06A_11R_A18T_07","TCGA_EE_A2GB_06A_11R_A18T_07","TCGA_EE_A2GC_06A_11R_A18T_07","TCGA_EE_A2GD_06A_11R_A18T_07","TCGA_EE_A2GE_06A_11R_A18T_07","TCGA_EE_A2GH_06A_11R_A18T_07","TCGA_EE_A2GI_06A_11R_A18T_07","TCGA_EE_A2GJ_06A_11R_A18U_07","TCGA_EE_A2GK_06A_11R_A18S_07","TCGA_EE_A2GL_06A_11R_A18S_07","TCGA_EE_A2GM_06B_11R_A18S_07","TCGA_EE_A2GN_06A_11R_A18S_07","TCGA_EE_A2GO_06A_11R_A18S_07","TCGA_EE_A2GP_06A_11R_A18S_07","TCGA_EE_A2GR_06A_11R_A18S_07","TCGA_EE_A2GS_06A_12R_A18S_07","TCGA_EE_A2GT_06A_12R_A18S_07","TCGA_EE_A2GU_06A_11R_A18T_07","TCGA_EE_A2M5_06A_12R_A18S_07","TCGA_EE_A2M6_06A_12R_A18S_07","TCGA_EE_A2M7_06A_11R_A18U_07","TCGA_EE_A2M8_06A_12R_A18S_07","TCGA_EE_A2MC_06A_12R_A18S_07","TCGA_EE_A2MD_06A_11R_A18T_07","TCGA_EE_A2ME_06A_11R_A18T_07","TCGA_EE_A2MF_06A_11R_A21D_07","TCGA_EE_A2MG_06A_11R_A18T_07","TCGA_EE_A2MH_06A_11R_A18S_07","TCGA_EE_A2MI_06A_11R_A18U_07","TCGA_EE_A2MJ_06A_11R_A18S_07","TCGA_EE_A2MK_06A_11R_A18S_07","TCGA_EE_A2ML_06A_11R_A18S_07","TCGA_EE_A2MM_06A_11R_A18S_07","TCGA_EE_A2MN_06A_11R_A18S_07","TCGA_EE_A2MP_06A_11R_A18S_07","TCGA_EE_A2MQ_06A_11R_A18S_07","TCGA_EE_A2MR_06A_11R_A18S_07","TCGA_EE_A2MS_06A_11R_A18S_07","TCGA_EE_A2MT_06A_11R_A18S_07","TCGA_EE_A2MU_06A_21R_A18S_07","TCGA_EE_A3AA_06A_11R_A18S_07","TCGA_EE_A3AB_06A_11R_A18S_07","TCGA_EE_A3AC_06A_11R_A18S_07","TCGA_EE_A3AD_06A_11R_A18S_07","TCGA_EE_A3AF_06A_11R_A18S_07","TCGA_EE_A3AG_06A_31R_A18S_07","TCGA_EE_A3AH_06A_11R_A18S_07","TCGA_EE_A3J3_06A_11R_A20F_07","TCGA_EE_A3J4_06A_11R_A20F_07","TCGA_EE_A3J5_06A_11R_A20F_07","TCGA_EE_A3J7_06A_11R_A20F_07","TCGA_EE_A3J8_06A_11R_A20F_07","TCGA_EE_A3JA_06A_11R_A20F_07","TCGA_EE_A3JB_06A_11R_A21D_07","TCGA_EE_A3JD_06A_11R_A20F_07","TCGA_EE_A3JE_06A_11R_A20F_07","TCGA_EE_A3JH_06A_11R_A21D_07","TCGA_EE_A3JI_06A_11R_A21D_07","TCGA_ER_A193_06A_12R_A18S_07","TCGA_ER_A194_01A_11R_A18U_07","TCGA_ER_A195_06A_11R_A18U_07","TCGA_ER_A196_01A_11R_A18T_07","TCGA_ER_A197_06A_32R_A18S_07","TCGA_ER_A198_06A_11R_A18T_07","TCGA_ER_A199_06A_11R_A18T_07","TCGA_ER_A19A_06A_21R_A18U_07","TCGA_ER_A19B_06A_11R_A18S_07","TCGA_ER_A19C_06A_11R_A18S_07","TCGA_ER_A19D_06A_11R_A18S_07","TCGA_ER_A19E_06A_11R_A18S_07","TCGA_ER_A19F_06A_11R_A18S_07","TCGA_ER_A19G_06A_11R_A18U_07","TCGA_ER_A19H_06A_12R_A18S_07","TCGA_ER_A19J_06A_11R_A18S_07","TCGA_ER_A19K_01A_21R_A18T_07","TCGA_ER_A19L_06A_12R_A18S_07","TCGA_ER_A19M_06A_61R_A239_07","TCGA_ER_A19N_06A_11R_A18S_07","TCGA_ER_A19O_06A_11R_A18S_07","TCGA_ER_A19P_06A_11R_A18S_07","TCGA_ER_A19Q_06A_11R_A18U_07","TCGA_ER_A19S_06A_11R_A18U_07","TCGA_ER_A19T_06A_11R_A18U_07","TCGA_ER_A19T_01A_11R_A18T_07","TCGA_ER_A19W_06A_41R_A239_07","TCGA_ER_A1A1_06A_11R_A18U_07","TCGA_ER_A2NB_01A_12R_A18S_07","TCGA_ER_A2NC_06A_11R_A18T_07","TCGA_ER_A2ND_06A_11R_A18T_07","TCGA_ER_A2NE_06A_21R_A18T_07","TCGA_ER_A2NF_06A_11R_A18T_07","TCGA_ER_A2NF_01A_11R_A18T_07","TCGA_ER_A2NG_06A_11R_A18T_07","TCGA_ER_A2NH_06A_11R_A18S_07","TCGA_ER_A3ES_06A_11R_A20F_07","TCGA_ER_A3ET_06A_11R_A20F_07","TCGA_ER_A3EV_06A_11R_A20F_07","TCGA_ER_A3PL_06A_11R_A239_07","TCGA_ER_A42H_01A_11R_A24X_07","TCGA_ER_A42K_06A_11R_A24X_07","TCGA_ER_A42L_06A_11R_A24X_07","TCGA_FR_A2OS_01A_11R_A21D_07","TCGA_FR_A3R1_01A_11R_A239_07","TCGA_FR_A3YN_06A_11R_A239_07","TCGA_FR_A3YO_06A_11R_A239_07","TCGA_FR_A44A_06A_11R_A24X_07","TCGA_FR_A69P_06A_21R_A311_07","TCGA_FR_A726_01A_11R_A32P_07","TCGA_FR_A728_01A_11R_A32P_07","TCGA_FR_A729_06A_11R_A352_07","TCGA_FR_A7U8_06A_21R_A352_07","TCGA_FR_A7U9_06A_11R_A352_07","TCGA_FR_A7UA_06A_32R_A352_07","TCGA_FR_A8YC_06A_11R_A37K_07","TCGA_FR_A8YD_06A_11R_A37K_07","TCGA_FR_A8YE_06A_11R_A37K_07","TCGA_FS_A1YW_06A_11R_A18T_07","TCGA_FS_A1YX_06A_11R_A18T_07","TCGA_FS_A1YY_06A_11R_A18T_07","TCGA_FS_A1Z0_06A_11R_A18T_07","TCGA_FS_A1Z3_06A_11R_A18T_07","TCGA_FS_A1Z4_06A_11R_A18T_07","TCGA_FS_A1Z7_06A_11R_A18T_07","TCGA_FS_A1ZA_06A_11R_A18T_07","TCGA_FS_A1ZB_06A_12R_A18S_07","TCGA_FS_A1ZC_06A_11R_A18T_07","TCGA_FS_A1ZD_06A_11R_A18T_07","TCGA_FS_A1ZE_06A_11R_A18T_07","TCGA_FS_A1ZF_06A_12R_A18S_07","TCGA_FS_A1ZG_06A_11R_A18T_07","TCGA_FS_A1ZH_06A_11R_A18T_07","TCGA_FS_A1ZJ_06A_12R_A18S_07","TCGA_FS_A1ZK_06A_11R_A18T_07","TCGA_FS_A1ZM_06A_12R_A18S_07","TCGA_FS_A1ZN_01A_11R_A18T_07","TCGA_FS_A1ZP_06A_11R_A18T_07","TCGA_FS_A1ZQ_06A_11R_A18U_07","TCGA_FS_A1ZR_06A_21R_A18U_07","TCGA_FS_A1ZS_06A_12R_A18T_07","TCGA_FS_A1ZT_06A_11R_A18U_07","TCGA_FS_A1ZU_06A_12R_A18T_07","TCGA_FS_A1ZW_06A_12R_A18T_07","TCGA_FS_A1ZY_06A_11R_A18S_07","TCGA_FS_A1ZZ_06A_11R_A18S_07","TCGA_FS_A4F0_06A_11R_A24X_07","TCGA_FS_A4F2_06A_11R_A24X_07","TCGA_FS_A4F4_06A_12R_A266_07","TCGA_FS_A4F5_06A_11R_A266_07","TCGA_FS_A4F8_06A_11R_A266_07","TCGA_FS_A4F9_06A_11R_A24X_07","TCGA_FS_A4FB_06A_11R_A266_07","TCGA_FS_A4FC_06A_11R_A24X_07","TCGA_FS_A4FD_06A_11R_A266_07","TCGA_FW_A3I3_06A_11R_A21D_07","TCGA_FW_A3R5_06A_11R_A239_07","TCGA_FW_A3TU_06A_11R_A239_07","TCGA_FW_A3TV_06A_11R_A239_07","TCGA_FW_A5DX_01A_11R_A27Q_07","TCGA_FW_A5DY_06A_11R_A311_07","TCGA_GF_A2C7_01A_11R_A18T_07","TCGA_GF_A3OT_06A_23R_A239_07","TCGA_GF_A4EO_06A_12R_A24X_07","TCGA_GF_A6C8_06A_12R_A311_07","TCGA_GF_A6C9_06A_11R_A311_07","TCGA_GF_A769_01A_32R_A32P_07","TCGA_GN_A262_06A_11R_A18T_07","TCGA_GN_A263_01A_11R_A18T_07","TCGA_GN_A264_06A_11R_A18U_07","TCGA_GN_A265_06A_21R_A18T_07","TCGA_GN_A266_06A_11R_A18T_07","TCGA_GN_A267_06A_21R_A18T_07","TCGA_GN_A268_06A_11R_A18T_07","TCGA_GN_A26A_06A_11R_A18T_07","TCGA_GN_A26C_01A_11R_A18T_07","TCGA_GN_A26D_06A_11R_A18T_07","TCGA_GN_A4U3_06A_11R_A32P_07","TCGA_GN_A4U4_06A_11R_A32P_07","TCGA_GN_A4U5_01A_11R_A32P_07","TCGA_GN_A4U7_06A_21R_A32P_07","TCGA_GN_A4U8_11A_11R_A32P_07","TCGA_GN_A4U8_06A_11R_A32P_07","TCGA_GN_A4U9_06A_11R_A32P_07","TCGA_GN_A8LK_06A_11R_A37K_07","TCGA_GN_A8LL_06A_21R_A37K_07","TCGA_GN_A8LN_01A_11R_A37K_07","TCGA_GN_A9SD_06A_11R_A40A_07","TCGA_HR_A2OG_06A_21R_A18U_07","TCGA_HR_A2OH_06A_11R_A18U_07","TCGA_HR_A5NC_01A_11R_A27Q_07","TCGA_IH_A3EA_01A_11R_A20F_07","TCGA_LH_A9QB_06A_11R_A38C_07","TCGA_OD_A75X_06A_12R_A32P_07","TCGA_QB_A6FS_06A_11R_A311_07","TCGA_QB_AA9O_06A_11R_A39D_07","TCGA_RP_A690_06A_11R_A311_07","TCGA_RP_A693_06A_13R_A311_07","TCGA_RP_A694_06A_11R_A311_07","TCGA_RP_A695_06A_11R_A311_07","TCGA_RP_A6K9_06A_41R_A352_07","TCGA_W3_A824_06A_21R_A352_07","TCGA_W3_A825_06A_11R_A352_07","TCGA_W3_A828_06A_11R_A352_07","TCGA_W3_AA1O_06A_11R_A38C_07","TCGA_W3_AA1Q_06A_11R_A38C_07","TCGA_W3_AA1R_06A_11R_A39D_07","TCGA_W3_AA1V_06B_11R_A40A_07","TCGA_W3_AA1W_06A_11R_A38C_07","TCGA_W3_AA21_06A_11R_A38C_07","TCGA_WE_A8JZ_06A_11R_A37K_07","TCGA_WE_A8K1_06A_21R_A37K_07","TCGA_WE_A8K4_01A_12R_A37K_07","TCGA_WE_A8K5_06A_11R_A37K_07","TCGA_WE_A8K6_06A_11R_A37K_07","TCGA_WE_A8ZM_06A_11R_A37K_07","TCGA_WE_A8ZN_06A_11R_A37K_07","TCGA_WE_A8ZO_06A_11R_A37K_07","TCGA_WE_A8ZQ_06A_41R_A37K_07","TCGA_WE_A8ZR_06A_11R_A37K_07","TCGA_WE_A8ZT_06A_11R_A37K_07","TCGA_WE_A8ZX_06A_11R_A37K_07","TCGA_WE_A8ZY_06A_11R_A37K_07","TCGA_WE_AA9Y_06A_12R_A38C_07","TCGA_WE_AAA0_06A_11R_A38C_07","TCGA_WE_AAA3_06A_11R_A38C_07","TCGA_WE_AAA4_06A_12R_A38C_07","TCGA_XV_A9VZ_01A_11R_A38C_07","TCGA_XV_A9W2_01A_11R_A39D_07","TCGA_XV_A9W5_01A_11R_A38C_07","TCGA_XV_AAZV_01A_11R_A40A_07","TCGA_XV_AAZW_01A_12R_A40A_07","TCGA_XV_AAZY_01A_12R_A40A_07","TCGA_XV_AB01_06A_12R_A40A_07","TCGA_YD_A89C_06A_11R_A37K_07","TCGA_YD_A9TA_06A_11R_A39D_07","TCGA_YD_A9TB_06A_12R_A40A_07","TCGA_YG_AA3N_01A_11R_A38C_07","TCGA_YG_AA3O_06A_11R_A38C_07","TCGA_YG_AA3P_06A_11R_A38C_07","TCGA_Z2_A8RT_06A_11R_A37K_07","TCGA_Z2_AA3S_06A_11R_A39D_07","TCGA_Z2_AA3V_06A_11R_A39D_07")
    sample_names_column <- "Sample"
    gene_names_column <- "Gene"
    groups_column <- "shortLetterCode"
    
    input_in_log_counts <- FALSE
    normalization_method <- "none"
    
    principal_component_on_x_axis_for_pca <- 1
    principal_component_on_y_axis_for_pca <- 2
    colors_for_plots <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")
    point_size_for_pca <- 2
    add_labels_to_pca <- FALSE
    labels_column <- "Sample"
    label_offset_y_for_pca <- 2
    label_offset_x_for_pca <- 2
    label_font_size_for_pca <- 3
    legend_position_for_pca <- "top"
    samples_to_rename_manually_on_pca <- c("")
    
    
    color_histogram_by_group <- TRUE
    set_min_max_for_x_axis_for_histogram <- FALSE 
    minimum_for_x_axis_for_histogram <- -1
    maximum_for_x_axis_for_histogram <- 1
    legend_position_for_histogram <- "top"
    legend_font_size_for_histogram <- 10
    number_of_histogram_legend_columns <- 6
    
    make_plots_interactive <- FALSE
    plot_correlation_matrix_heatmap <- TRUE
    
    number_of_image_rows <- 2

    ##--------------- ##
    ## Error Messages ##
    ## -------------- ##

    
    ## --------- ##
    ## Functions ##
    ## --------- ##

    getourrandomcolors<-function(k){
        seed=10
        n <- 2e3
        ourColorSpace <- colorspace::RGB(runif(n), runif(n), runif(n))
        ourColorSpace <- as(ourColorSpace, "LAB")
        currentColorSpace <- ourColorSpace@coords
        # Set iter.max to 20 to avoid convergence warnings.
        set.seed(seed)
        km <- kmeans(currentColorSpace, k, iter.max=20)
        return( unname(hex(LAB(km$centers))))
    }

    make_heatmap <- function(counts_matrix, metadata,colorval) {
        mat <- as.matrix(counts_matrix) 
        tcounts=t(mat)
        d=Dist(tcounts,method="correlation",diag=TRUE)
        dend = rev(dendsort(as.dendrogram(hclust( d,method="average"))))
        m=as.matrix(d)
        sample_metadata <- metadata
        rownames(sample_metadata) = sample_metadata[[sample_names_column]]
        idx = as.factor(sample_metadata[rownames(m),groups_column])
        col = colorval
        cols <- col[idx]
        new.palette=colorRampPalette(c("blue","green","yellow"),space="rgb")
  
    mk<-function(){
        if(length(colnames(m))>20){
            par(mar=c(0,0,0,0))
            heatmap.2(m,labRow = NA, 
                        labCol = NA,
                        col=new.palette(20),
                        trace="none",
                        colRow = col[idx], 
                        colCol = col[idx],
                        rowDendrogram=dend,
                        colDendrogram=dend,
                        RowSideColors = col[idx],
                        ColSideColors = col[idx],
                        dendrogram = "row",
                        cexRow=3,
                        cexCol=3,
                        margins=c(0,0),   
                        lmat=rbind( c(0,0,2),c(4,1,3) ,c(0,5,6) ), 
                        lhei=c(.2,4,2), 
                        lwid=c(1, .2,4 ), 
                        key.par=list(mgp=c(1.75, .5, 0), 
                        mar=c(7, 2, 3.5, 0), 
                        cex.axis=.1, 
                        cex.lab=3, 
                        cex.main=1, 
                        cex.sub=1),
                        key.xlab = "Correlation",
                        key.ylab="Count",
                        key.title=" ")       
        } else {
            heatmap.2(m,col=new.palette(20),
                        trace="none",
                        colRow = col[idx], 
                        colCol = col[idx],
                        rowDendrogram=dend,
                        colDendrogram=dend,
                        RowSideColors = col[idx],
                        ColSideColors = col[idx],
                        dendrogram = "row",
                        cexRow=3,cexCol=3,margins=c(4,1),  
                        lmat=rbind( c(0,0,2),c(4,1,3) ,c(0,5,6) ), 
                        lhei=c( .2,4,2), 
                        lwid=c(1, .2,4),
                        key.par=list(mgp=c(1.75, .5, 0), mar=c(7, 2, 3.5, 0), cex.axis=.1, cex.lab=3, cex.main=1, cex.sub=1),
                        key.xlab = "Correlation",
                        key.ylab="Count",
                        key.title=" ")
            }
        }
  
            tg<-mk()
            grid.echo(mk)
            gh1<-grid.grab()
            mklegend<-function(){
            plot.new()
            legend(x="top", legend=levels(idx), col=col[as.factor(levels(idx))],pch=15,x.intersp=3,bty ="n",cex=2)
            }
        grid.echo(mklegend )
        gh2<-grid.grab()
        lay <- c(1,3)
        grid.newpage()
        grid.arrange(gh1,gh2,nrow=1,widths=c(unit(1000, "bigpts"),unit(300, "bigpts")))
        gh<-grid.grab()
        return(gh)
    }

    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##

    samples_to_include <- samples_to_include[samples_to_include != gene_names_column]
    samples_to_include <- samples_to_include[samples_to_include != "Gene"]
    samples_to_include <- samples_to_include[samples_to_include != "GeneName"]
    df.filt <- counts_matrix[,samples_to_include]
    gene_names <- NULL
    gene_names$GeneID <- counts_matrix[,1]
    
    sample_metadata <- sample_metadata[match(colnames(df.filt),sample_metadata[[sample_names_column]]),] #First match sample metadata to counts matrix
    sample_metadata <- sample_metadata[rowSums(is.na(sample_metadata)) != ncol(sample_metadata), ] # Remove empty rows
    sample_metadata <- sample_metadata[, colSums(is.na(sample_metadata)) == 0] #Remove empty columns
    rownames(sample_metadata) <- sample_metadata[[sample_names_column]]

    df.filt <- df.filt[,match(sample_metadata[[sample_names_column]],colnames(df.filt))] #Match counts matrix columns to sample metadata
    
    #If input is in log space, linearize
    if(input_in_log_counts == TRUE){
        x <- DGEList(counts=2^df.filt, genes=gene_names)
    } else {
        x <- DGEList(counts=df.filt, genes=gene_names)     
    }

    v <- voom(x,normalize=normalization_method)
    rownames(v$E) <- v$genes$GeneID
    as.data.frame(v$E) %>% rownames_to_column("Gene") -> df.voom
    print(paste0("Total number of genes included: ", nrow(df.voom)))

    #Start PCA Plot:
    
    edf <- v$E
    tedf <- t(edf)
    tedf <- tedf[, colSums(is.na(tedf)) != nrow(tedf)]
    tedf <- tedf[, apply(tedf, 2, var) != 0]
    pca <- prcomp(tedf, scale.=T)
    
    pcx <- paste0("PC",principal_component_on_x_axis_for_pca)
    pcy <- paste0("PC",principal_component_on_y_axis_for_pca)
    pca.df <- as.data.frame(pca$x) %>% dplyr::select(.data[[pcx]], .data[[pcy]])
    pca.df$group <- sample_metadata[[groups_column]]
    pca.df$sample <- sample_metadata[[labels_column]]
    perc.var <- (pca$sdev^2/sum(pca$sdev^2))*100
    perc.var <- formatC(perc.var,format = "g",digits=4)
    pc.x.lab <- paste0(pcx," ", perc.var[principal_component_on_x_axis_for_pca],"%")
    pc.y.lab <- paste0(pcy," ", perc.var[principal_component_on_y_axis_for_pca],"%")
    labelpos <- pca.df
    labelpos$mean_y <- pca.df[[pcy]]+label_offset_y_for_pca
    labelpos$mean_x <- pca.df[[pcx]]+label_offset_x_for_pca
    pca.df$xdata <- pca.df[[pcx]]
    pca.df$ydata <- pca.df[[pcy]]

    # Manual changes to sample names
    replacements = samples_to_rename_manually_on_pca

    if (!is.null(replacements)) {
        if (replacements != c("")) {
            for (x in replacements) {
                old <- strsplit(x, ": ?")[[1]][1]
                new <- strsplit(x, ": ?")[[1]][2]
                pca.df$sample <- ifelse(pca.df$sample==old, new, pca.df$sample)
            }
        }
    }

    colorlist <- c("#5954d6","#e1562c","#b80058","#00c6f8","#d163e6","#00a76c","#ff9287","#008cf9","#006e00","#796880","#FFA500","#878500")
    names(colorlist) <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")
    if(length(colors_for_plots) == 0){
        colors_for_plots <- c("indigo","carrot","lipstick","turquoise","lavender","jade","coral","azure","green","rum","orange","olive")
    }
    colorval <- colorlist[colors_for_plots]
    colorval <- unname(colorval) #remove names which affect ggplot

    if (length(unique(sample_metadata[[groups_column]])) > length(colorval)) {
        ## Original color-picking code.
        k=length(unique(sample_metadata[[groups_column]]))-length(colorval)
        more_cols<- getourrandomcolors(k) 
        colorval <- c(colorval , more_cols)
    }

    if (add_labels_to_pca){
    g <- ggplot(pca.df, aes(x=xdata, y=ydata)) +
        theme_bw() +
        theme(legend.title=element_blank()) +
        theme(legend.position=legend_position_for_pca) +
        geom_point(aes(color=group), size=point_size_for_pca) +
        geom_text(data=labelpos, aes(x=labelpos$mean_x, y=labelpos$mean_y, 
            label=sample, color=group, vjust="inward", hjust="inward"), size=label_font_size_for_pca, show.legend=FALSE) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
        scale_colour_manual(values = colorval) +
        xlab(pc.x.lab) + ylab(pc.y.lab)
    } else {
    g <- ggplot(pca.df, aes(x=xdata, y=ydata)) +
        theme_bw() +
        theme(legend.title=element_blank()) +
        theme(legend.position=legend_position_for_pca) +
        geom_point(aes(color=group,text=sample), size=point_size_for_pca) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
        scale_colour_manual(values = colorval) +
        xlab(pc.x.lab) + ylab(pc.y.lab)    
    }

    par(mfrow = c(2,1))

    #Histogram Plot:
    
    df.m <- melt(edf,id.vars=c(gene_names_column))
    df.m = dplyr::rename(df.m,sample=Var2)

    if(set_min_max_for_x_axis_for_histogram == TRUE){
        xmin = minimum_for_x_axis_for_histogram
        xmax = maximum_for_x_axis_for_histogram
    } else {
        xmin = min(df.m$value)
        xmax = max(df.m$value)
    }

    if(color_histogram_by_group == TRUE){
        df.m %>% mutate(colgroup = sample_metadata[sample,groups_column]) -> df.m
        df.m = df.m[complete.cases(df.m[, "colgroup"]),]
        df.m$colgroup = gsub("\\s","_",df.m$colgroup)
        df.m$colgroup = factor(df.m$colgroup, levels=unique(df.m$colgroup))
        print(unique(df.m$sample))

        # plot Density 
        g2 = ggplot(df.m, aes(x=value, group=sample)) + 
            geom_density(aes(colour = colgroup, linetype = colgroup)) +
            xlab("Filtered Counts") + ylab("Density") +
            theme_bw() +
            theme(legend.position=legend_position_for_histogram,legend.text = element_text(size = legend_font_size_for_histogram)) + 
            ggtitle("Frequency Histogram") +
            xlim(xmin,xmax) +
            scale_linetype_manual(values=rep(c('solid', 'dashed','dotted','twodash'),40)) +
            scale_colour_manual(values=colorval)
    } else {
        
        df.m$sample = sample_metadata[df.m$sample,labels_column]
        n=length(unique(df.m$sample))
        cols<- getourrandomcolors(n) 
        
        g2 = ggplot(df.m, aes(x=value, group=sample)) + 
            geom_density(aes(colour = sample, linetype = sample)) +
            xlab("Filtered Counts") + ylab("Density") +
            theme_bw() +
            theme(legend.position=legend_position_for_histogram,legend.text = element_text(size = legend_font_size_for_histogram)) +  
            ggtitle("Frequency Histogram") +
            xlim(xmin,xmax) +
            scale_linetype_manual(values=rep(c('solid', 'dashed','dotted','twodash'),n)) +
            scale_colour_manual(values=cols)#+
            guides(linetype = guide_legend(ncol = number_of_histogram_legend_columns))
    }

    dev.off()

    imageWidth = 3000
    imageHeight = 1500*2
    dpi = 300

    png(
      filename="Normalized_counts_SKCM.png",
      width=imageWidth,
      height=imageHeight,
      units="px",
      pointsize=4,
      bg="white",
      res=dpi,
      type="cairo")

    if(plot_correlation_matrix_heatmap == TRUE){
        if(make_plots_interactive == TRUE){
            p1=(g)%>%ggplotly(tooltip = c("sample","group"))
            p2=(g2+theme(legend.position = "none")) %>%ggplotly(tooltip = c("sample"))
            fig=subplot(p1,p2,which_layout = 'merge',margin=.05,shareX = F,shareY = F,titleY = T,titleX = T,widths=c(.5,.5),nrows = 1)
            fig=fig %>% layout(title = 'Interactive PCA and Histogram')
            print(fig)
        } else {
            require(gridExtra)
            gh<-make_heatmap(df.filt[,samples_to_include],sample_metadata,colorval)
            grid.arrange(g,g2,gh, nrow=number_of_image_rows)
            dev.off()
        }  
    } else {
        if(make_plots_interactive == TRUE){
            p1=(g)%>%ggplotly(tooltip = c("sample","group"))
            p2=(g2+theme(legend.position = "none")) %>%ggplotly(tooltip = "sample" )
            fig=subplot(p1,p2,which_layout = 'merge',margin=.05,shareX = F,shareY = F,titleY = T,titleX = T,widths=c(.5,.5),nrows = 1)
            fig=fig %>% layout(title = 'Interactive PCA and Histogram')
            print(fig)
        } else {
            grid.arrange(g,g2, nrow=number_of_image_rows)
            dev.off()
        }
        }    
    
    return(df.voom)
}

#################################################
## Global imports and functions included below ##
#################################################

print("template_function_Normalized_counts_SKCM.R #########################################################################")
library(plotly);library(ggplot2);library(jsonlite);
var_Filtered_Counts_SKCM<-readRDS("var_Filtered_Counts_SKCM.rds")
var_TCGA_SKCM_harmonized_clinical_data<-readRDS("var_TCGA_SKCM_harmonized_clinical_data.rds")
invisible(graphics.off())
var_Normalized_counts_SKCM<-Normalized_counts_SKCM(var_Filtered_Counts_SKCM,var_TCGA_SKCM_harmonized_clinical_data)
invisible(graphics.off())
saveRDS(var_Normalized_counts_SKCM,"var_Normalized_counts_SKCM.rds")
