wd <- dirname(dirname(getwd()))
source(paste0(wd,"/mission_control/treasure_map.R"))
library(tidyverse)

ready <- readRDS(paste0(TMP_DIR, "validation-ready.Rds"))

summary_table_1 <- (
  ready
    %>% group_by(Study, cohort, tissue)
    %>% summarise( tot = n(), 
                   age_mn = round(mean( as.numeric(age), na.rm = TRUE)),
                   age_me = round(1.96*round(sd(as.numeric(age), na.rm = TRUE))/sqrt(tot), 1),
                 bor_ct = sum(!is.na(bor)),
                 responders_ct = sum(bor, na.rm = TRUE),
                 non_responders_ct = bor_ct - responders_ct,
                 os_ct = sum(!is.na(os) & os!= 0),
                 tmb_ct = sum(!is.na(tmb)), 
                 rna_ct = sum(!is.na(tcell)), 
                 pretrt_ct = sum(!is.na(pretreat)),
                 complete_genomic_data_ct = sum(!is.na(tmb) & !is.na(tcell) & !is.na(pretreat)),
                 complete_data_ct = sum(!is.na(tmb) & !is.na(tcell) & !is.na(pretreat) & !is.na(bor)  & !is.na(os) & os != 0)
))

1491 - 339

summary_table_1 %>% group_by(Study == "HMF-CPCT") %>% summarise(sum(tot),
                                                                sum(bor_ct),
                                                                sum(os_ct),
                                                                sum(complete_genomic_data_ct),
                                                                sum(complete_data_ct)
                                                               )

summary_table_2 <- (
  ready
    %>% group_by(Study, cohort)
    %>% summarise( tot = n(), 
                   age_mn = round(mean( as.numeric(age), na.rm = TRUE)),
                   age_me = round(1.96*round(sd(as.numeric(age), na.rm = TRUE))/sqrt(tot), 1),
                   bor_ct = sum(!is.na(bor)),
                   responders_ct = sum(bor, na.rm = TRUE),
                   non_responders_ct = bor_ct - responders_ct,
                   os_ct = sum(!is.na(os) & os!= 0),
                   tmb_ct = sum(!is.na(tmb)), 
                   rna_ct = sum(!is.na(tcell)), 
                   pretrt_ct = sum(!is.na(pretreat)),
                   complete_data_ct = sum(!is.na(tmb) & !is.na(tcell) & !is.na(pretreat) & !is.na(bor)  & !is.na(os) & os != 0)
))

write.csv(summary_table_1, paste0(FIG_FINAL_DIR, "9_study_summaries_table1.csv"))
write.csv(summary_table_2, paste0(FIG_FINAL_DIR, "9_study_summaries_table2.csv"))
