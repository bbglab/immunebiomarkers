wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
source(paste0(wd,"/mission_control/helpers/figures/exhaustive_prep.R"))
source(paste0(wd,"/mission_control/helpers/figures/exhaustive_settings.R"))
library(tidyverse)
library(stringr)
library(RColorBrewer)

boom <- readRDS(paste0(TMP_DIR,"exhaustive-plots-base.Rds"))

main_filters <- (    
    boom 
        %>% filter( dataset == "all",
                    !feature %in% c("age", "pretreat",'purity','tmb','tcell','tgfb','prolif'),
                    feature != "isofox_gene_set_mariathan_CD_8_T_effector", 
                    feature != "isofox_gene_set_mariathan_Immune_Checkpoint", 
                    feature != "cibersort_TR4_mix_r",
                    !feature %in% c("isofox_gene_set_mariathan_tcga", "isofox_gene_set_KEGG_GLYCOSAMINOGLYCAN_BIOSYNTHESIS_HEPARAN_SULFATE", "isofox_gene_set_KEGG_GLYCOSPHINGOLIPID_BIOSYNTHESIS_GANGLIO_SERIES", "isofox_gene_set_HALLMARK_GLYCOLYSIS", "isofox_gene_set_KEGG_REGULATION_OF_ACTIN_CYTOSKELETON", "isofox_gene_set_mariathan_Mismatch_Repair","isofox_gene_set_KEGG_PANCREATIC_CANCER","isofox_gene_set_CELL_PROLIFERATION_GO_0008283", "isofox_gene_set_KEGG_SMALL_CELL_LUNG_CANCER","isofox_gene_set_KEGG_PROSTATE_CANCER") ,
                    !grepl("_10", feature),
                    !grepl("_05", feature),
                    !grepl("vhio", tolower(feature)), 
                    !grepl("rand", tolower(feature)), 
                    !grepl("per ", tolower(feature)), 
                    !grepl("high", tolower(feature)), 
                    !grepl("battle", tolower(feature)),
                    (covariates %in% c("age_biopsy_purity_tissue") & model %in% c("bor", "os")) | ### removed PFS
                    (covariates %in% c("residuals") & model %in% c( "os" ))
                  )
        %>% mutate( feature_group = ifelse( feature ==  "hla_lilac_aneuploidy_score", "cnv", feature_group))
        %>% mutate( model = ifelse( covariates == "residuals", "os_res", model))
        %>% mutate( feature_group = ifelse(grepl("gene_set", feature), "gene_set", feature_group))
)

finer_grouper <- function( discovery_group = "X", feature_group = "X" ){
    if( discovery_group == "RNA: T-cell"){
        if( feature_group == "isofox"){
            "T-cell genes"
        } else if (feature_group == "gene_set"){
            "T-cell gene sets"
        } else {
            feature_group
        }
    } else if( discovery_group == "RNA: TGFB"){
        if( feature_group == "isofox"){
            "TGFB genes"
        } else if (feature_group == "gene_set"){
            "TGFB gene sets"
        } else {
            feature_group
        }
    } else if( discovery_group == "RNA: Proliferation"){
        if( feature_group == "isofox"){
            "Proliferation genes"
        } else if (feature_group == "gene_set"){
            "Proliferation gene sets"
        } else {
            feature_group
        }
    } else if ( discovery_group == "RNA: Remaining"){
        if( feature_group == "isofox"){
            "Remaining genes"
        } else if (feature_group == "gene_set"){
            "Remaining gene sets"
        } else {
            feature_group
        }
    } else {
        feature_group
    }   
}

finer_groups <- (
    main_filters 
        %>% mutate(discovery_group = as.character(discovery_group))
        %>% rowwise() 
        %>% mutate( go_group = finer_grouper(discovery_group, feature_group))
        %>% group_by( go_group, model ) 
        %>% mutate( rk = row_number(p_val))
        %>% ungroup()
)

top_features <- (
    finer_groups 
        %>% filter(
            feature %in% c("hla_HLA_all_LOH", "hla_HLA_all_tumor_heterozygous") |
            model %in% c("os_res") & discovery_group %in% c("RNA: Proliferation", "RNA: TGFB") & rk < 13 |
            model %in% c("bor") & discovery_group %in% c("RNA: T-cell") & rk < 13 |
            model %in% c("bor") & discovery_group %in% c("RNA: Remaining") & rk < 8 |
            model %in% c("bor") & discovery_group %in% c("Clinical") & rk < 12 |
            model %in% c("bor") & !discovery_group %in% c("Clinical","RNA: Proliferation", "RNA: TGFB", "RNA: T-cell", "RNA: Remaining") & rk < 9 | 
            model %in% c("bor") & feature_group %in% c("somatic") & rk < 15 ) %>% pull(feature)
)

base <- finer_groups %>% filter( feature %in% top_features)

feature_levels <- (
     base 
        %>% filter(
             model == "bor" & !discovery_group %in% c("RNA: Proliferation", "RNA: TGFB") | 
             model == "os_res" & discovery_group %in% c("RNA: Proliferation", "RNA: TGFB"))
        %>% arrange(go_group, desc(rk))
        %>% pull(feature)
)

base_camp <- base %>% mutate(feature = factor(feature, levels = feature_levels))

z_alpha <- -qnorm( base_camp$by_05_fdr[1]/2)
    
camp1 <- (
    base_camp
      %>% mutate( est = ifelse( (model == "os_res") & ((cor_tmb > .4 & cor_tmb != 1)| (cor_tcell > .4 & cor_tcell != 1) | (cor_pretreat > .4 & cor_pretreat != 1)), NA, est))
      %>% mutate( ci_low = est - z_alpha*se, ci_high = est + z_alpha*se)
      %>% select(feature, est, ci_low, ci_high, p_val, go_group, feature_group, model, contains("cor"))
      %>% arrange(go_group, est)
)

cor_map <- list("cor_pretreat" = "Pretr", 
                "cor_tcell" = "T-cell",
                "cor_tmb" = "TMB",
                "cor_tgfb" = "TGFB",
                "cor_prolif" = "Prolif")
mod_map <- list("os" = "OS", 
                "pfs" = "PFS", 
                "bor" = "BOR", 
                "os_res" = "OS Res")
gp_map <- list("cibersort" = "Cibersort",
               "clinical" = "Clinical",
               "cnv" = "CNV Summary",
               "cnv.region" = "CNV Region",
               "driver_cnv" = "CNV Drivers",
               "driver_somatic" = "Mutation Driver",
               "hla" = "HLA",
               "sig" = "SBS Sigs",
               "somatic" = "TMB Variants",
               "somatic.gene" = "Gene Mutations",
               "sv" = "SV Summary")

feature_map <- list(
    "cibersort_LM22_Dendritic.cells.resting" = "Dendritic Resting",	
	 "cibersort_LM22_Macrophages.M2" = "Macrophage M2",
	 "cibersort_TR4_immune" = "Immune cells",
	 "cibersort_LM22_T.cells.CD4.memory.activated" = "T-cells CD4 activated",
	 "cibersort_LM22_T.cells.CD8" = "T-cells CD8",
	 "cibersort_LM22_T.cells.follicular.helper" = "T-cells Follicular",
	 "cibersort_LM22_T.cells.gamma.delta" = "T-cells gamma delta",
	 "cibersort_LM22_Macrophages.M1" = "Macrophage M1",
	 "clinical_number_pretreatment" = "Number pretreatment",
	 "clinical_meta_hasSystemicPreTreatment2" = "Prior Systemic therapy",
	 "clinical_pre_treated" = "Prior therapy",
	 "clinical_pre_contains_Chemotherapy" = "Prior Chemotherapy",
	 "clinical_meta_hasRadiotherapyPreTreatment" = "Prior Radiotherapy",
	 "clinical_pre_contains_Immunotherapy" = "Prior Immunotherapy",
	 "clinical_pre_to_post_treatment_time" = "Pre-to-Post treatment time",
	 "clinical_post_contains_Targeted" = "Targeted therapy given",
	 "cnv_scna" = "SCNA",
	 "cnv_summary_wholeGenomeDuplication" = "WGD",
	 "cnv_summary_ploidy" = "Ploidy",
	 "cnv_summary_diploidProportion" = "Diploid %",
	 "cnv_copy_loss_burden" = "Copy Loss Burden",
	 "hla_lilac_aneuploidy_score" = "Anueploidy",
	 "cnv_summary_polyclonalProportion" = "Polyclonal %",
	 "cnv.region_loh_chr2.p25.3" = "LOH Chr2.p25.3",
	 "cnv.region_loh_chr2.p25.1.p25.2" = "LOH Chr2.p25.1",
	 "cnv.region_loh_chr2.p25.2" = "LOH Chr2.p25.2",
	 "cnv.region_cn_chr9.p24.1" = "CN Chr9.p24.1",
	 "cnv.region_cn_chr9.p23.p24.1" = "CN Chr9.p23.p24.1",
	 "cnv.region_loh_chr10.q26.13.q26.2" = "LOH Chr10.q26.13",
	 "cnv.region_loh_chr10.q26.3" = "LOH Chr10.q26.3",
	 "cnv.region_loh_chr15.q26.2" = "LOH Chr15.q26.2",
	 "driver_LRP1B_DEL" = "LRP1B Del",
	 "driver_FGF3_AMP" = "FGF3 Amp",
	 "driver_MYC_AMP" = "MYC Amp",
	 "driver_CCND1_AMP" = "CCND1 Amp",
	 "driver_TERT_AMP" = "TERT Amp",
	 "driver_PTEN_DEL" = "PTEN Del",
	 "driver_CDKN2A_DEL" = "CDKN2A Del",
	 "driver_PTPRD_DEL" = "PTPRD Del",
	 "driver_KRAS_MUTATION" = "KRAS",
	 "driver_PIK3CA_MUTATION" = "PIK3CA",
	 "driver_NRAS_MUTATION" = "NRAS",
	 "driver_PBRM1_MUTATION" = "PBRM1",
	 "driver_PTEN_MUTATION" = "PTEN",
	 "driver_VHL_MUTATION" = "VHL",
	 "driver_TP53_MUTATION" = "TP53",
	 "driver_TERT_MUTATION" = "TERT",
	 "hla_HLA_contains_B27" = "HLA contains B27",
	 "hla_HLA_contains_B08" = "HLA contains B08",
	 "hla_lilac_targeted_escape" = "HLA targeted escape (lilac)",
	 "hla_HLA_contains_B62" = "HLA contains B62",
	 "hla_lilac_imbalance" = "HLA imbalance (lilac)",
	 "hla_HLA_contains_A24" = "HLA contains A24",
	 "hla_HLA_contains_A02" = "HLA contains A02",
	 "hla_HLA_contains_B44" = "HLA contains B44",
	 "isofox_gene_set_KEGG_PROGESTERONE_MEDIATED_OOCYTE_MATURATION" = "KEGG OOCYTE Maturation",
	 "isofox_gene_set_HALLMARK_MITOTIC_SPINDLE" = "Hallmark Mitotic Spindle",
     "isofox_gene_set_HALLMARK_NOTCH_SIGNALING" = "Hallmark Notch Signaling",
	 "isofox_gene_set_KEGG_OOCYTE_MEIOSIS" = "KEGG Meiosis",
	 "isofox_gene_set_mariathan_Cell_cycle" = "Cell Cycle",
	 "isofox_gene_set_KEGG_CELL_CYCLE" = "KEGG Cell Cycle",
	 "isofox_gene_set_prolif_cluster" = "Proliferation cluster",
	 "isofox_gene_set_HALLMARK_G2M_CHECKPOINT" = "G2M Checkpoint",
	 "isofox_gene_set_prolif" = "Proliferation gene set",
	 "isofox_ANLN" = "ANLN",
	 "isofox_MKI67" = "MKI67",
	 "isofox_KIF11" = "KIF11",
	 "isofox_CENPF" = "CENPF",
	 "isofox_CENPA" = "CENPA",
	 "isofox_CCNA2" = "CCNA2",
	 "isofox_TOP2A" = "TOP2A",
	 "isofox_SKA1" = "SKA1",
	 "isofox_gene_set_mariathan_EMT3" = "EMT3",
	 "isofox_gene_set_KEGG_GLYCOSPHINGOLIPID_BIOSYNTHESIS_GANGLIO_SERIES" = "Glycosphingolipid Biosynthesis",
	 "isofox_gene_set_KEGG_MELANOGENESIS" = "KEGG Melanogenesis",
	 "isofox_gene_set_KEGG_GLYCOSAMINOGLYCAN_BIOSYNTHESIS_HEPARAN_SULFATE" = "Glycosaminoglycan Biosynthesis",
	 "isofox_gene_set_HALLMARK_KRAS_SIGNALING_DN" = "Hallmark KRAS Signaling",
	 "isofox_gene_set_KEGG_CYSTEINE_AND_METHIONINE_METABOLISM" = "KEGG Cysteine Metabolism",
	 "isofox_gene_set_KEGG_NITROGEN_METABOLISM" = "KEGG Nitrogen Metabolism",
	 "isofox_gene_set_KEGG_PROTEASOME" = "KEGG Proteasome",
	 "isofox_LATS2" = "LATS2",
	 "isofox_CTTN" = "CTTN",
	 "isofox_ELL3" = "ELL3",
	 "isofox_BAIAP2" = "BAIAP2",
	 "isofox_CD276" = "CD276",
	 "isofox_WTIP" = "WTIP",
	 "isofox_ADORA2A" = "ADORA2A",
	 "isofox_TLR3" = "TLR3",
	 "sig_SBS8" = "SBS8",
	 "sig_SBS3" = "SBS3",
	 "sig_SBS1" = "SBS1",
	 "sig_SBS5" = "SBS5",
	 "sig_SBS4" = "SBS4",
	 "sig_SBS40" = "SBS40",
	 "sig_SBS38" = "SBS38",
	 "sig_SBS93" = "SBS93",
	 "somatic_summary_tmbStatus" = "TMB high",
	 "somatic_summary_tmlStatus" = "TML high",
	 "somatic_summary_tmbPerMb" = "TMB per MB",
	 "somatic_TMB_clonal" = "TMB clonal",
	 "somatic_TMB_damage_pathway" = "TMB damage pathway",
	 "somatic_summary_tml" = "Tumor Mutational Load",
	 "somatic_TMB_exome" = "TMB Exome",
	 "somatic_TMB" = "TMB",
	 "somatic_TMB_sbs" = "TMB SBS",
	 "somatic_TMB_indel" = "TMB Indel",
	 "somatic_TMB_mbs" = "TMB MBS",
	 "somatic_TMB_SBS10b" = "TMB SBS10b",
	 "somatic_TMB_frameshift" = "TMB frameshift",
	 "somatic_TMB_dbs" = "TMB DBS",
	 "somatic.gene_TMEM117.mb" = "TMEM117",
	 "somatic.gene_KIRREL3.mb" = "KIRREL3",
	 "somatic.gene_EML5.mb" = "EML5",
	 "somatic.gene_EIPR1.mb" = "EIPR1",
	 "somatic.gene_DEPDC1B.mb" = "DEPDC1B.",
	 "somatic.gene_NEDD4L.mb" = "NEDD4L",
	 "somatic.gene_SOD2.mb" = "SOD2",
	 "somatic.gene_ITIH6.mb" = "ITIH6",
	 "sv_clusters" = "Clusters",
	 "sv_svs" = "SVs",
	 "sv_breakend" = "Breakends",
	 "sv_summary_svTumorMutationalBurden" = "SV TMB",
	 "sv_links" = "Links",
	 "sv_fusion" = "Fusions",
	 "isofox_gene_set_t_cell_effector" = "T-cell effector",
	 "isofox_gene_set_immune_checkpoint_genes" = "Immune checkpoint genes",
	 "isofox_gene_set_tcell_cluster" = "T-cell cluster",
	 "isofox_gene_set_t_cell_gep_18" = "T-cell GEP 18",
	 "isofox_gene_set_cd8_t_effector" = "CD8 T-cell effector",
	 "isofox_gene_set_infiltrate" = "Infiltrate",
	 "isofox_gene_set_t_cell_gep_6" = "T-cell GEP 6",
	 "isofox_gene_set_cyt" = "CYT",
	 "isofox_TIGIT" = "TIGIT",
	 "isofox_PLA2G2D" = "PLA2G2D",
	 "isofox_GBP5" = "GBP5",
	 "isofox_CD8A" = "CD8A",
	 "isofox_CXCL9" = "CXCL9",
	 "isofox_UBD" = "UBD",
	 "isofox_STAT1" = "STAT1",
	 "isofox_CALHM6" = "CALHM6",
	 "isofox_gene_set_KEGG_ECM_RECEPTOR_INTERACTION" = "KEGG ECM receptor interaction",
	 "isofox_gene_set_tgfb_cluster" = "TGFB cluster",
	 "isofox_gene_set_KEGG_PATHWAYS_IN_CANCER" = "KEGG Pathways in cancer",
	 "isofox_gene_set_KEGG_SMALL_CELL_LUNG_CANCER" = "",
	 "isofox_gene_set_KEGG_PROSTATE_CANCER" = "",
	 "isofox_gene_set_HALLMARK_ANGIOGENESIS" = "Hallmark Angiogenesis",
	 "isofox_gene_set_HALLMARK_APICAL_SURFACE" = "Hallmark Apical Surface",
	 "isofox_gene_set_CELL_PROLIFERATION_GO_0008283" = "GO Cell proliferation",
	 "isofox_C3orf36" = "C3orf36",
	 "isofox_CCDC3" = "CCDC3",
	 "isofox_COL4A2" = "COL4A2",
	 "isofox_COL4A1" = "COL4A1",
	 "isofox_NID1" = "NID1",
	 "isofox_CD248" = "CD248",
	 "isofox_PALLD" = "PALLD",
	 "isofox_S1PR3" = "S1PR3",
     "isofox_gene_set_HALLMARK_COAGULATION" = "Hallmark Coagulation",
     "isofox_gene_set_KEGG_GLIOMA" = "KEGG Glioma",
     "isofox_gene_set_Pan_TBRS" = "Pan-TBRS",
     "hla_HLA_all_LOH" = "HLA LOH",
     "isofox_NKG7" = "NKG7",
     "isofox_IRF1" = "IRF1",
     "isofox_GBP4" = "GBP4",
     "isofox_TCF4" = "TCF4",
     "isofox_SPARC" = "SPARC",
     "isofox_GGT5" = "GGT5",
     "isofox_KNL1" = "KNL1",
     "isofox_IQGAP3" = "IQGAP3",
     "isofox_GTSE1" = "GTSE1",
     "isofox_gene_set_F_TBRS" = "F-TBRS",
     "isofox_gene_set_KEGG_FOCAL_ADHESION" = "KEGG Focal Adhesion",
     "isofox_gene_set_HALLMARK_TGF_BETA_SIGNALING" = "Hallmark TGFB",
     "isofox_gene_set_mariathan_APM" = "APM",
     "isofox_gene_set_12_chemokine" = "Chemokine 12",
     "isofox_gene_set_KEGG_PRIMARY_IMMUNODEFICIENCY" = "KEGG primary immunodeficiency",
     "isofox_gene_set_HALLMARK_E2F_TARGETS" = "Hallmark E2F Targets",
     "isofox_gene_set_KEGG_MISMATCH_REPAIR" = "KEGG Mismatch Repair",
     "isofox_gene_set_KEGG_PYRIMIDINE_METABOLISM" = "KEGG Pyrimidine Metabolism",
     "isofox_gene_set_KEGG_CIRCADIAN_RHYTHM_MAMMAL" = "KEGG Circadian Rhythm",
     "isofox_gene_set_mariathan_tcga" = "TCGA",
     "clinical_pre_contains_Targeted"  = "Prior Targeted Therapy",
     "clinical_meta_tumorPurity" = "Tumor Purity",
     "clinical_age_at_treatment_start" = "Age",
     "isofox_NDC80" = "NDC80", 
     "isofox_CD3E" = "CD3E",
     "isofox_VSTM4" = "VSTM4",
     "hla_HLA_all_tumor_heterozygous" = "HLA Heterzygosity",
     "isofox_gene_set_KEGG_HOMOLOGOUS_RECOMBINATION" = "KEGG Homologous Recombination",
     "isofox_gene_set_KEGG_ALLOGRAFT_REJECTION" = "KEGG Allograft Rejection",
     "isofox_gene_set_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "Hallmark EMT"
)

name_mapper <- function(i, map){
    if( i %in% names(map)){
        map[[i]]
    } else {
        i
    }
}

camp2 <- (
  camp1 
    %>% rowwise() 
    %>% mutate(model = factor(name_mapper(as.character(model), mod_map), 
                          levels = c("BOR", "PFS", "OS", "OS Res")),
               go_group = factor(name_mapper(as.character(go_group), gp_map), 
                          levels = c("Clinical", 
                                     "HLA", 
                                     "TMB Variants", 
                                     "Gene Mutations", 
                                     "Mutation Driver", 
                                     "SBS Sigs",
                                     "Cibersort", 
                                     "T-cell gene sets", 
                                     "T-cell genes", 
                                     "TGFB gene sets", 
                                     "TGFB genes", 
                                     "Proliferation gene sets", 
                                     "Proliferation genes", 
                                     "Remaining gene sets", 
                                     "Remaining genes",
                                     "CNV Summary", 
                                     "CNV Region", 
                                     "CNV Drivers", 
                                     "SV Summary")),
               feature = name_mapper(as.character(feature), feature_map))
    %>% ungroup()
)

better_or_worse <- function( model, ci_low, ci_high){
    if( is.na(ci_low) | is.na(ci_high) | is.na(model)){
        "none"
    } else if( model == "BOR"){
        if(ci_low > 0){
            "better"
        } else if (ci_high < 0){
            "worse"
        } else {
            "none"
        }
    } else {
        if(ci_low > 0){
            "worse"
        } else if (ci_high < 0){
            "better"
        } else {
            "none"
        }
    }
}

camp2 <- camp2 %>% rowwise() %>% mutate( dir = better_or_worse( model, ci_low, ci_high) ) %>% ungroup()

res <- c("TGFB gene sets", "TGFB genes","Proliferation gene sets","Proliferation genes")
feature_order <- c( camp2 %>% filter(go_group %in% res,model == "BOR") %>% arrange(go_group, est) %>% pull(feature),
   camp2 %>% filter(!go_group %in% res,model == "OS Res") %>% arrange(go_group, est) %>% pull(feature))

#camp2 <- camp2 %>% mutate(feature = factor(feature, levels = feature_order))

camp3 <- (
  camp2 
    %>% filter(model == "OS") 
    %>% select(feature, go_group, contains("cor"), -cor_purity)
    %>% gather("latent", "cor", -feature, -go_group)
    %>% rowwise()
    %>% mutate( latent = factor(name_mapper(latent, cor_map), levels = c("TMB", "T-cell", "TGFB", "Prolif", "Pretr")))
)

cor_grouper <- function(cor){
    if( abs(cor) > .5 ){
        "high"
    } else if ( abs(cor) > .3){
        "medium"
    } else {
        "low"
    }
}

camp3 <- camp3 %>% rowwise() %>% mutate(alpha = cor_grouper(cor)) %>% ungroup()

gps <- list(
    "cln_hla" = c("Clinical", "HLA"), 
    "cnv_sv" = c("CNV Summary", "CNV Region", "SV Summary", "CNV Drivers"),
    "somatic" = c("Mutation Driver", "Gene Mutations", "TMB Variants", "SBS Sigs"),
    "rna" = c("Cibersort", 
              "T-cell genes", "T-cell gene sets", 
              "TGFB genes", "TGFB gene sets", 
              "Proliferation genes", "Proliferation gene sets",
              "Remaining genes", "Remaining gene sets")
)

climb <- function( i = "cln_hla"){
    list( "ests" = camp2 %>% filter(go_group %in% gps[[i]]), "cors" = camp3 %>% filter(go_group %in% gps[[i]]))
}

ready <- list()

for( i in names(gps)){
    ready[[i]] <- climb(i)
}

saveRDS( ready, paste0(TMP_DIR, "exhaustive-plots-highlights2.Rds"))
