wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
library(tidyverse)

hmf <- readRDS(paste0(TMP_DIR, "validation-hmf-go.Rds")) %>% mutate(extra = "bbg_cluster", extra2 = NA)
vhio <- readRDS(paste0(TMP_DIR, "validation-vhio-go.Rds")) %>% mutate(extra = NA, extra2 = NA)
inspire <- readRDS(paste0(TMP_DIR, "validation-inspire-go.Rds")) %>% mutate(extra = NA, extra2 = NA)
mariathan <- readRDS(paste0(TMP_DIR, "validation-mariathan-go.Rds"))
ravi <- readRDS(paste0(TMP_DIR, "validation-ravi-go.Rds"))
parker <- readRDS(paste0(TMP_DIR, "validation-parker-go.Rds"))
lyon <- readRDS(paste0(TMP_DIR, "validation-lyon-go.Rds"))

tissue_update <- function(i) {
        if( i == "D: Melanoma"){
            "skin"
        } else if (i == "A: Head and Neck") {
            "hncc"
        } else if (i == "E: Mixed") {
            "other"
        } else if (i == "B: Breast"){
            "breast"
        } else if (i == "C: Ovary"){
            "ovary"
        }
}
inspire$tissue_full <- unlist(lapply(inspire$tissue, tissue_update))

inspire <- (
    inspire 
        %>% mutate(
            tissue = ifelse( tissue_full == "skin", "skin", "other"),
            model_apply = ifelse( tissue_full == "skin", "skin", "all"),
            cohort = "all"
        )
)

vhio <- (
    vhio 
      %>% mutate(
        model_apply =  ifelse( tissue == "Colon", "all", tolower(tissue)),
        cohort = "all",
        purity = NA,
        tcell_set = tcell, 
        prolif_set = prolif,
        tgfb_set = tgfb,
        tcell_cluster5 = tcell, 
        prolif_cluster5 = prolif,
        tgfb_cluster5 = tgfb,  
    )
)

hmf <- hmf %>% mutate( model_apply = tissue, cohort = "all")

mariathan <- mariathan %>% mutate( model_apply = "bladder", cohort = "all")

ravi <- ravi %>% mutate( model_apply = "lung", cohort = "all")

parker <- parker %>% mutate( model_apply = "skin") %>% select(-sample.id)

lyon <- lyon %>% mutate( model_apply = ifelse(tissue == "lung", "lung", "all"))

validation_ready <- (
    rbind(vhio, inspire, hmf, mariathan, ravi , parker, lyon)
        %>% mutate(
            tissue = tolower(tissue),
            tmb_bin = ifelse(exp(tmb)-1 > 10, 1, 0),
            complete = !is.na(tmb) & !is.na(pretreat) & !is.na(tcell),
            complete_dna = !is.na(tmb),
            complete_rna = !is.na(tcell),
            complete_pretreat = !is.na(pretreat), 
            qc = NA,
            tcell = as.numeric(tcell),
            prolif = as.numeric(prolif),
            tgfb = as.numeric(tgfb),
            tcell_cluster5 = as.numeric(tcell_cluster5),
            prolif_cluster5 = as.numeric(prolif_cluster5),
            tgfb_cluster5 = as.numeric(tgfb_cluster5),
            tcell_set = as.numeric(tcell_set),
            prolif_set = as.numeric(prolif_set),
            tgfb_set = as.numeric(tgfb_set),
            purity = as.numeric(purity)
        )
)

validation_go <- (
  validation_ready 
    %>% group_by(Study) 
    %>% mutate_at(vars(all_of(c("tmb", "tcell", "tgfb", "prolif", "tcell_set", "prolif_set", "tgfb_set","pretreat"))), scale)
    %>% ungroup()
    %>% mutate(study_tissue_cohort = paste0(Study,"-", tissue, "-", cohort))
    %>% group_by(Study, tissue)
    %>% mutate(
               pretreat=ifelse(is.na(pretreat),mean(pretreat,na.rm=TRUE),pretreat),
               tmb=ifelse(is.na(tmb),median(tmb,na.rm=TRUE),tmb), 
               tcell=ifelse(is.na(tcell),median(tcell,na.rm=TRUE),tcell),
               tgfb=ifelse(is.na(tgfb),median(tgfb,na.rm=TRUE),tgfb),
               prolif=ifelse(is.na(prolif),median(prolif,na.rm=TRUE),prolif)
              )
    %>% ungroup()
    %>% mutate( 
            pretreat = ifelse(pretreat == "NaN", mean(pretreat), pretreat),
            pretreat = ifelse(is.na(pretreat),mean(pretreat,na.rm=TRUE),pretreat),
            tmb= ifelse(is.na(tmb),median(tmb,na.rm=TRUE),tmb)
    )
    %>% mutate(Study_cohort = ifelse(Study == "HMF-CPCT", paste0(Study,"-", tissue), Study))
    %>% mutate(Study_cohort = ifelse(Study %in% c("Lyon", "PARKER"), paste0(Study,"-", cohort), Study_cohort))
)

validation_go <- (
 validation_go %>% bind_rows(
    validation_go %>% filter( !grepl("HMF-CPCT", Study)) %>% mutate(Study_cohort = "External Studies"),
    validation_go %>% filter( grepl("PARKER", Study_cohort)) %>% mutate(Study_cohort = "PARKER"),
    validation_go %>% filter( grepl("Lyon", Study_cohort)) %>% mutate(Study_cohort = "Lyon"),
    validation_go %>% filter( Study == "HMF-CPCT") %>% mutate(Study_cohort = "HMF-CPCT")
) %>% filter(os_days > 0)
  %>% unique() 
  %>% mutate(purity = as.numeric(purity))
)

table(validation_go$Study_cohort)

saveRDS(validation_ready, paste0(TMP_DIR, "validation-ready.Rds"))
saveRDS(validation_go, paste0(TMP_DIR, "validation-go.Rds"))
