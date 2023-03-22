
               
                   

install.packages(c("ggridges", "ggplot2", "reshape2", "plyr", "dplyr", "tidyverse","plotly","GGaly","rio","ggthemes","foreach","parallel","doParallel","foreach","EnvStats","rstatix","psych"))
install.packages("tidyverse")


#install.packages("rlang")


library(readxl)
library(ggplot2)
library(reshape2) 
library(plyr)
library(dplyr)
library(ggthemes) # Charger
library(tidyverse)
library(plotly)
library(GGally)
library(rio)
library(data.table)
library(ggridges)
library(doParallel)


install.packages("remotes")

remotes::install_github('rpkgs/gg.layers')
library(gg.layers)

install.packages("devtools")

devtools::install_github("associatedpress/aptheme")

library(aptheme)

library(psych)
library(rstatix)
library(EnvStats)



# Set the working directory of R studio to Outputs_Fish or change absolute path

results_total <- fread(file="result_fish_2_23_561982_Total.csv")



sub_cm_eutro= results_total %>% select(starts_with("genericAH"))

sub_GWP = results_total %>% select(starts_with("GWP100AH"))
sub_HTPinf = results_total %>% select(starts_with("HTPinfAH"))
sub_FETPinf = results_total %>% select(starts_with("FETPinfAH"))
sub_eutrophication = results_total %>% select(starts_with("eutrophicationAH"))
sub_TETPinf = results_total %>% select(starts_with("TETPinfAH"))
sub_TAP100 = results_total %>% select(starts_with("TAP100AH"))
sub_PMFP = results_total %>% select(starts_with("PMFPAH"))
sub_FE = results_total %>% select(starts_with("FEPAH"))
sub_ME = results_total %>% select(starts_with("MEPAH"))

sub_ODPinf = results_total %>% select(starts_with("ODPinfAH"))

sub_ECO_FCR = results_total %>% select(starts_with("Economic FCR"))
sub_BIO_FCR = results_total %>% select(starts_with("Biological FCR"))




sub_BIO_FCR$paired = 1:nrow(sub_BIO_FCR)
sub_ECO_FCR$paired = 1:nrow(sub_ECO_FCR)

sub_GWP$paired = 1:nrow(sub_GWP)
sub_HTPinf$paired =1:nrow(sub_HTPinf)
sub_FETPinf$paired =1:nrow(sub_FETPinf)
sub_eutrophication$paired =1:nrow(sub_eutrophication)
sub_TETPinf$paired =1:nrow(sub_TETPinf)
sub_TAP100$paired =1:nrow(sub_TAP100)
sub_PMFP$paired =1:nrow(sub_PMFP)
sub_FE$paired =1:nrow(sub_FE)
sub_ME$paired =1:nrow(sub_ME)
sub_ODPinf$paired =1:nrow(sub_ODPinf)
sub_cm_eutro$paired =1:nrow(sub_cm_eutro)




sub_GW_melted=melt(sub_GWP,measure.vars = colnames(sub_GWP)[-length(colnames(sub_GWP))],variable.name = 'Scenario',value.name='Score')
sub_HTPinf_melted=melt(sub_HTPinf,measure.vars = colnames(sub_HTPinf)[-length(colnames(sub_HTPinf))],variable.name = 'Scenario',value.name='Score')
sub_FETPinf_melted=melt(sub_FETPinf,measure.vars = colnames(sub_FETPinf)[-length(colnames(sub_FETPinf))],variable.name = 'Scenario',value.name='Score')
sub_eutrophication_melted=melt(sub_eutrophication,measure.vars = colnames(sub_eutrophication)[-length(colnames(sub_eutrophication))],variable.name = 'Scenario',value.name='Score')
sub_TETPinf_melted=melt(sub_TETPinf,measure.vars = colnames(sub_TETPinf)[-length(colnames(sub_TETPinf))],variable.name = 'Scenario',value.name='Score')

sub_TAP100_melted=melt(sub_TAP100,measure.vars = colnames(sub_TAP100)[-length(colnames(sub_TAP100))],variable.name = 'Scenario',value.name='Score')
sub_PMFP_melted=melt(sub_PMFP,measure.vars = colnames(sub_PMFP)[-length(colnames(sub_PMFP))],variable.name = 'Scenario',value.name='Score')


sub_ODPinf_melted=melt(sub_ODPinf,measure.vars = colnames(sub_ODPinf)[-length(colnames(sub_ODPinf))],variable.name = 'Scenario',value.name='Score')


sub_FE_melted=melt(sub_FE,measure.vars = colnames(sub_FE)[-length(colnames(sub_FE))],variable.name = 'Scenario',value.name='Score')
sub_ME_melted=melt(sub_ME,measure.vars = colnames(sub_ME)[-length(colnames(sub_ME))],variable.name = 'Scenario',value.name='Score')



sub_BIO_FCR_melted = melt(sub_BIO_FCR,measure.vars = colnames(sub_BIO_FCR)[-length(colnames(sub_BIO_FCR))],variable.name = 'Scenario',value.name='Score')
sub_ECO_FCR_melted = melt(sub_ECO_FCR,measure.vars = colnames(sub_ECO_FCR)[-length(colnames(sub_ECO_FCR))],variable.name = 'Scenario',value.name='Score')

sub_cm_eutro_melted = melt(sub_cm_eutro,measure.vars = colnames(sub_cm_eutro)[-length(colnames(sub_cm_eutro))],variable.name = 'Scenario',value.name='Score')




sub_ECO_FCR_melted$Scenario <- recode_factor(sub_ECO_FCR_melted$Scenario,
                                        `Economic FCR` = "Current",
                                        `Economic FCRdict_update_sce_A` = "No mortality",
                                        `Economic FCRdict_update_sce_B` = "No mortality and minimum FCR",
                                        `Economic FCRdict_update_sce_C` = "Current mortality and minimum FCR",
                                        `Economic FCRdict_update_sce_stage_hatch` = "Mort. hatchery",
                                        `Economic FCRdict_update_sce_stage_fry` = "Mort. fry-fingerling",
                                        `Economic FCRdict_update_sce_stage_goit` = "Mort. On-growing 1",
                                        `Economic FCRdict_update_sce_stage_godk` = "Mort. On-growing DK",
                                        `Economic FCRdict_update_sce_stage_sfdk1` = "Mort. Seafarm 1",
                                        `Economic FCRdict_update_sce_stage_sfdk2` = "Mort. Seafarm 2",
                                        `Economic FCRdeter` = "Deterministic")






sub_BIO_FCR_melted$Scenario <- recode_factor(sub_BIO_FCR_melted$Scenario,
                                             `Biological FCR` = "Current",
                                             `Biological FCRdict_update_sce_A` = "No mortality",
                                             `Biological FCRdict_update_sce_B` = "No mortality and minimum FCR",
                                             `Biological FCRdict_update_sce_C` = "Current mortality and minimum FCR",
                                             `Biological FCRdict_update_sce_stage_hatch` = "Mort. hatchery",
                                             `Biological FCRdict_update_sce_stage_fry` = "Mort. fry-fingerling",
                                             `Biological FCRdict_update_sce_stage_goit` = "Mort. On-growing 1",
                                             `Biological FCRdict_update_sce_stage_godk` = "Mort. On-growing DK",
                                             `Biological FCRdict_update_sce_stage_sfdk1` = "Mort. Seafarm 1",
                                             `Biological FCRdict_update_sce_stage_sfdk2` = "Mort. Seafarm 2",
                                             `Biological FCRdeter` = "Deterministic")




sub_GW_melted$Scenario <- recode_factor(sub_GW_melted$Scenario,
                                        `GWP100AH` = "Current",
                                        `GWP100AHdict_update_sce_A` = "No mortality",
                                        `GWP100AHdict_update_sce_B` = "No mortality and minimum FCR",
                                        `GWP100AHdict_update_sce_C` = "Current mortality and minimum FCR",
                                        `GWP100AHdict_update_sce_stage_hatch` = "Mort. hatchery",
                                        `GWP100AHdict_update_sce_stage_fry` = "Mort. fry-fingerling",
                                        `GWP100AHdict_update_sce_stage_goit` = "Mort. On-growing 1",
                                        `GWP100AHdict_update_sce_stage_godk` = "Mort. On-growing DK",
                                        `GWP100AHdict_update_sce_stage_sfdk1` = "Mort. Seafarm 1",
                                        `GWP100AHdict_update_sce_stage_sfdk2` = "Mort. Seafarm 2",
                                        `GWP100AHdeter` = "Deterministic")




sub_HTPinf_melted$Scenario <- recode_factor(sub_HTPinf_melted$Scenario,
                                        `HTPinfAH` = "Current",
                                        `HTPinfAHdict_update_sce_A` = "No mortality",
                                        `HTPinfAHdict_update_sce_B` = "No mortality and minimum FCR",
                                        `HTPinfAHdict_update_sce_C` = "Current mortality and minimum FCR",
                                        `HTPinfAHdict_update_sce_stage_hatch` = "Mort. hatchery",
                                        `HTPinfAHdict_update_sce_stage_fry` = "Mort. fry-fingerling",
                                        `HTPinfAHdict_update_sce_stage_goit` = "Mort. On-growing 1",
                                        `HTPinfAHdict_update_sce_stage_godk` = "Mort. On-growing DK",
                                        `HTPinfAHdict_update_sce_stage_sfdk1` = "Mort. Seafarm 1",
                                        `HTPinfAHdict_update_sce_stage_sfdk2` = "Mort. Seafarm 2",
                                        `HTPinfAHdeter` = "Deterministic")






sub_FETPinf_melted$Scenario <- recode_factor(sub_FETPinf_melted$Scenario,
                                             `FETPinfAH` = "Current",
                                             `FETPinfAHdict_update_sce_A` = "No mortality",
                                             `FETPinfAHdict_update_sce_B` = "No mortality and minimum FCR",
                                             `FETPinfAHdict_update_sce_C` = "Current mortality and minimum FCR",
                                             `FETPinfAHdict_update_sce_stage_hatch` = "Mort. hatchery",
                                             `FETPinfAHdict_update_sce_stage_fry` = "Mort. fry-fingerling",
                                             `FETPinfAHdict_update_sce_stage_goit` = "Mort. On-growing 1",
                                             `FETPinfAHdict_update_sce_stage_godk` = "Mort. On-growing DK",
                                             `FETPinfAHdict_update_sce_stage_sfdk1` = "Mort. Seafarm 1",
                                             `FETPinfAHdict_update_sce_stage_sfdk2` = "Mort. Seafarm 2",
                                             `FETPinfAHdeter` = "Deterministic")




sub_eutrophication_melted$Scenario <- recode_factor(sub_eutrophication_melted$Scenario,
                                            `eutrophicationAH` = "Current",
                                            `eutrophicationAHdict_update_sce_A` = "No mortality",
                                            `eutrophicationAHdict_update_sce_B` = "No mortality and minimum FCR",
                                            `eutrophicationAHdict_update_sce_C` = "Current mortality and minimum FCR",
                                            `eutrophicationAHdict_update_sce_stage_hatch` = "Mort. hatchery",
                                            `eutrophicationAHdict_update_sce_stage_fry` = "Mort. fry-fingerling",
                                            `eutrophicationAHdict_update_sce_stage_goit` = "Mort. On-growing 1",
                                            `eutrophicationAHdict_update_sce_stage_godk` = "Mort. On-growing DK",
                                            `eutrophicationAHdict_update_sce_stage_sfdk1` = "Mort. Seafarm 1",
                                            `eutrophicationAHdict_update_sce_stage_sfdk2` = "Mort. Seafarm 2",
                                            `eutrophicationAHdeter` = "Deterministic")




sub_TETPinf_melted$Scenario <- recode_factor(sub_TETPinf_melted$Scenario,
                                                    `TETPinfAH` = "Current",
                                                    `TETPinfAHdict_update_sce_A` = "No mortality",
                                                    `TETPinfAHdict_update_sce_B` = "No mortality and minimum FCR",
                                                    `TETPinfAHdict_update_sce_C` = "Current mortality and minimum FCR",
                                                    `TETPinfAHdict_update_sce_stage_hatch` = "Mort. hatchery",
                                                    `TETPinfAHdict_update_sce_stage_fry` = "Mort. fry-fingerling",
                                                    `TETPinfAHdict_update_sce_stage_goit` = "Mort. On-growing 1",
                                                    `TETPinfAHdict_update_sce_stage_godk` = "Mort. On-growing DK",
                                                    `TETPinfAHdict_update_sce_stage_sfdk1` = "Mort. Seafarm 1",
                                                    `TETPinfAHdict_update_sce_stage_sfdk2` = "Mort. Seafarm 2",
                                                    `TETPinfAHdeter` = "Deterministic")




sub_TAP100_melted$Scenario <- recode_factor(sub_TAP100_melted$Scenario,
                                             `TAP100AH` = "Current",
                                             `TAP100AHdict_update_sce_A` = "No mortality",
                                             `TAP100AHdict_update_sce_B` = "No mortality and minimum FCR",
                                             `TAP100AHdict_update_sce_C` = "Current mortality and minimum FCR",
                                             `TAP100AHdict_update_sce_stage_hatch` = "Mort. hatchery",
                                             `TAP100AHdict_update_sce_stage_fry` = "Mort. fry-fingerling",
                                             `TAP100AHdict_update_sce_stage_goit` = "Mort. On-growing 1",
                                             `TAP100AHdict_update_sce_stage_godk` = "Mort. On-growing DK",
                                             `TAP100AHdict_update_sce_stage_sfdk1` = "Mort. Seafarm 1",
                                             `TAP100AHdict_update_sce_stage_sfdk2` = "Mort. Seafarm 2",
                                             `TAP100AHdeter` = "Deterministic")


sub_PMFP_melted$Scenario <- recode_factor(sub_PMFP_melted$Scenario,
                                            `PMFPAH` = "Current",
                                            `PMFPAHdict_update_sce_A` = "No mortality",
                                            `PMFPAHdict_update_sce_B` = "No mortality and minimum FCR",
                                            `PMFPAHdict_update_sce_C` = "Current mortality and minimum FCR",
                                            `PMFPAHdict_update_sce_stage_hatch` = "Mort. hatchery",
                                            `PMFPAHdict_update_sce_stage_fry` = "Mort. fry-fingerling",
                                            `PMFPAHdict_update_sce_stage_goit` = "Mort. On-growing 1",
                                            `PMFPAHdict_update_sce_stage_godk` = "Mort. On-growing DK",
                                            `PMFPAHdict_update_sce_stage_sfdk1` = "Mort. Seafarm 1",
                                            `PMFPAHdict_update_sce_stage_sfdk2` = "Mort. Seafarm 2",
                                            `PMFPAHdeter` = "Deterministic")





sub_ODPinf_melted$Scenario <- recode_factor(sub_ODPinf_melted$Scenario,
                                         `ODPinfAH` = "Current",
                                         `ODPinfAHdict_update_sce_A` = "No mortality",
                                         `ODPinfAHdict_update_sce_B` = "No mortality and minimum FCR",
                                         `ODPinfAHdict_update_sce_C` = "Current mortality and minimum FCR",
                                         `ODPinfAHdict_update_sce_stage_hatch` = "Mort. hatchery",
                                         `ODPinfAHdict_update_sce_stage_fry` = "Mort. fry-fingerling",
                                         `ODPinfAHdict_update_sce_stage_goit` = "Mort. On-growing 1",
                                         `ODPinfAHdict_update_sce_stage_godk` = "Mort. On-growing DK",
                                         `ODPinfAHdict_update_sce_stage_sfdk1` = "Mort. Seafarm 1",
                                         `ODPinfAHdict_update_sce_stage_sfdk2` = "Mort. Seafarm 2",
                                         `ODPinfAHdeter` = "Deterministic")





sub_ME_melted$Scenario <- recode_factor(sub_ME_melted$Scenario,
                                            `MEPAH` = "Current",
                                            `MEPAHdict_update_sce_A` = "No mortality",
                                            `MEPAHdict_update_sce_B` = "No mortality and minimum FCR",
                                            `MEPAHdict_update_sce_C` = "Current mortality and minimum FCR",
                                            `MEPAHdict_update_sce_stage_hatch` = "Mort. hatchery",
                                            `MEPAHdict_update_sce_stage_fry` = "Mort. fry-fingerling",
                                            `MEPAHdict_update_sce_stage_goit` = "Mort. On-growing 1",
                                            `MEPAHdict_update_sce_stage_godk` = "Mort. On-growing DK",
                                            `MEPAHdict_update_sce_stage_sfdk1` = "Mort. Seafarm 1",
                                            `MEPAHdict_update_sce_stage_sfdk2` = "Mort. Seafarm 2",
                                            `MEPAHdeter` = "Deterministic")


sub_FE_melted$Scenario <- recode_factor(sub_FE_melted$Scenario,
                                        `FEPAH` = "Current",
                                        `FEPAHdict_update_sce_A` = "No mortality",
                                        `FEPAHdict_update_sce_B` = "No mortality and minimum FCR",
                                        `FEPAHdict_update_sce_C` = "Current mortality and minimum FCR",
                                        `FEPAHdict_update_sce_stage_hatch` = "Mort. hatchery",
                                        `FEPAHdict_update_sce_stage_fry` = "Mort. fry-fingerling",
                                        `FEPAHdict_update_sce_stage_goit` = "Mort. On-growing 1",
                                        `FEPAHdict_update_sce_stage_godk` = "Mort. On-growing DK",
                                        `FEPAHdict_update_sce_stage_sfdk1` = "Mort. Seafarm 1",
                                        `FEPAHdict_update_sce_stage_sfdk2` = "Mort. Seafarm 2",
                                        `FEPAHdeter` = "Deterministic")





sub_cm_eutro_melted$Scenario <- recode_factor(sub_cm_eutro_melted$Scenario,
                                        `genericAH` = "Current",
                                        `genericAHdict_update_sce_A` = "No mortality",
                                        `genericAHdict_update_sce_B` = "No mortality and minimum FCR",
                                        `genericAHdict_update_sce_C` = "Current mortality and minimum FCR",
                                        `genericAHdict_update_sce_stage_hatch` = "Mort. hatchery",
                                        `genericAHdict_update_sce_stage_fry` = "Mort. fry-fingerling",
                                        `genericAHdict_update_sce_stage_goit` = "Mort. On-growing 1",
                                        `genericAHdict_update_sce_stage_godk` = "Mort. On-growing DK",
                                        `genericAHdict_update_sce_stage_sfdk1` = "Mort. Seafarm 1",
                                        `genericAHdict_update_sce_stage_sfdk2` = "Mort. Seafarm 2",
                                        `genericAHdeter` = "Deterministic")
















sub_TETPinf_melted$Category = "TETinf~(kg~1.4-DC-eq.)"

sub_TAP100_melted$Category = "TA100~(kg~SO[2]-eq.)"

sub_PMFP_melted$Category = "PM~(kg~PM10-eq.)"

sub_ME_melted$Category = "ME~(kg~N-eq.)"

sub_FE_melted$Category = "FE~(kg~P-eq.)"

sub_FETPinf_melted$Category = "FETinf~(kg~1.4-DC-eq.)"

sub_eutrophication_melted$Category = "Eutro.~(kg~N-eq.)"

sub_ODPinf_melted$Category = "ODinf.~(kg~CFC-11-eq.)"  


sub_HTPinf_melted$Category = "HTinf~(kg~1.4-DC-eq.)"

sub_GW_melted$Category = "GW100~(kg~CO[2]-eq)"


sub_cm_eutro_melted$Category  = "Eutro.~(PO[4]-eq.)"


sub_ECO_FCR_melted$Category="Economic"

sub_BIO_FCR_melted$Category="Biological"









total <- rbind(
                 sub_TETPinf_melted,
                 sub_TAP100_melted,
                 sub_PMFP_melted,
                 sub_ME_melted,
                 sub_ODPinf_melted,
                 sub_FE_melted,
                 sub_FETPinf_melted,
                 sub_HTPinf_melted,
                 sub_GW_melted,
                 sub_cm_eutro_melted)



total_with_FCR <- rbind(
  sub_TETPinf_melted,
  sub_TAP100_melted,
  sub_PMFP_melted,
  sub_ME_melted,
  sub_ODPinf_melted,
  sub_FE_melted,
  sub_FETPinf_melted,
  sub_HTPinf_melted,
  sub_GW_melted,
  sub_BIO_FCR_melted,
  sub_ECO_FCR_melted,
  sub_cm_eutro_melted)



total_FCR <- rbind(
 
  sub_BIO_FCR_melted,
  sub_ECO_FCR_melted)




# 
# total_with_FCR$Category <- recode_factor(total_with_FCR$Category,
#                                              `eutrophication` = "Eutro.")
# 
# 
# total$Category <- recode_factor(total$Category,
#                                          `eutrophication` = "Eutro.")











#total_best_incworst = total[ which(total$Scenario %in% c('Incumbent', "No mortality and minimum FCR", "No mortality", "Incumbent mortality and minimum FCR")),]
total_best_incworst = total[ which(total$Scenario %in% c('Current', "No mortality and minimum FCR", "No mortality", "Current mortality and minimum FCR")),]

total_best_incworst_FCR = total_with_FCR[which(total_with_FCR$Scenario %in% c('Current', "No mortality and minimum FCR", "No mortality", "Current mortality and minimum FCR")),]


total_best_incworst_FCR_only = total_FCR[which(total_FCR$Scenario %in% c('Current', "No mortality and minimum FCR", "No mortality", "Current mortality and minimum FCR")),]




total_best_incworst_mean_FCR_only= total_best_incworst_FCR_only %>% group_by(Scenario,Category)%>% 
  mutate(Mean= mean(`Score`)) %>% # calculate mean for plotting as well
  ungroup()

total_best_incworst_mean_FCR_only = unique(total_best_incworst_mean_FCR_only[,c('Scenario','Category','Mean')])



# Print bio and ECO FCR alone

total_best_incworst_mean_FCR_only$FCR = total_best_incworst_mean_FCR_only$Mean



# OOOK



total_best_incworst_mean_FCR_only %>%
  ggplot(aes(x=Scenario, y=FCR,shape=Category,color=Category))+ 
  geom_point(stat="identity",size=4,alpha=0.8) +theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(size=15))+
  ylim(0, 1.5)
ggsave('FCRs.jpeg', width = 7, height = 2,dpi=600)






total_best_incworst_mean_median= total_best_incworst_FCR %>% group_by(Scenario,Category)%>% 
  mutate(Mean= mean(`Score`),Median=median(`Score`)) %>% # calculate mean for plotting as well
  ungroup()


total_best_incworst_mean = unique(total_best_incworst_mean_median[,c('Scenario','Category','Mean')])
total_best_incworst_median = unique(total_best_incworst_mean_median[,c('Scenario','Category','Median')])


total_best_incworst_tosave = unique(total_best_incworst_mean_median[,c('Scenario','Category','Mean','Median')])

# expot this as csv
write.csv(total_best_incworst_tosave,"Mean_and_median_scores_opportunity_costs_2.csv",row.names = TRUE)







# Best, incumbnt and worst vertical



total_best_incworst %>%
  ggplot(aes(Scenario,Score, fill=Scenario)) + 
  facet_grid(rows=vars(Category),scales="free",labeller = label_parsed)+
  geom_boxplot2(width.errorbar = 0.001,) +theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size=15),
        strip.text.y = element_text(size = 10,angle=0),
        axis.text.y = element_text(size=9))
ggsave('total_dd22_3.jpeg', width = 7, height = 9,dpi=600)











### Sensitivity OAT mortality



total_sensitivity_mortality = total[ which(total$Scenario %in% c('Current', "Mort. hatchery",
                                                                 "Mort. fry-fingerling",
                                                                 "Mort. On-growing 1",
                                                                 "Mort. On-growing DK",
                                                                 "Mort. Seafarm 1",
                                                                 "Mort. Seafarm 2")),]





total_sensitivity_mean_median= total_sensitivity_mortality %>% group_by(Scenario,Category)%>% 
  mutate(Mean= mean(`Score`),Median=median(`Score`)) %>% # calculate mean for plotting as well
  ungroup()


total_sensi_mean = unique(total_sensitivity_mean_median[,c('Scenario','Category','Mean')])
total_sensi_median = unique(total_sensitivity_mean_median[,c('Scenario','Category','Median')])


total_sensi_tosave = unique(total_sensitivity_mean_median[,c('Scenario','Category','Mean','Median')])












total_sensitivity_mortality %>%
  ggplot(aes(Scenario,Score)) + 
  facet_grid(rows=vars(Category),scales="free",labeller = label_parsed)+
  geom_boxplot2(width.errorbar = 0.001, fill="blue",alpha=0.5) +
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=9),
        
        text = element_text(size=15),
        strip.text.y = element_text(size = 9))
ggsave('total_sensiti_2.jpeg', width = 4, height = 9,dpi=600)






#Subselection

total_sensitivity_mortality_sub = subset(total_sensitivity_mortality, !(Category %in% c("TETinf~(kg~1.4-DC-eq.)",
                                                                                        "TA100~(kg~SO[2]-eq.)",
                                                                                        "FETinf~(kg~1.4-DC-eq.)",
                                                                                        "ME~(kg~N-eq.)",
                                                                                        "FE~(kg~P-eq.)",
                                                                                        "PM~(kg~PM10-eq.)",
                                                                                        "Eutro.~(kg~N-eq.)",
                                                                                        "ODinf.~(kg~CFC-11-eq.)")))




total_sensitivity_mean_median_sub= total_sensitivity_mortality_sub %>% group_by(Scenario,Category)%>% 
  mutate(Mean= mean(`Score`),Median=median(`Score`)) %>% # calculate mean for plotting as well
  ungroup()


total_sensi_mean_sub = unique(total_sensitivity_mean_median_sub[,c('Scenario','Category','Mean')])
total_sensi_median_sub = unique(total_sensitivity_mean_median_sub[,c('Scenario','Category','Median')])


total_sensi_tosave_sub = unique(total_sensitivity_mean_median_sub[,c('Scenario','Category','Mean','Median')])








total_sensitivity_mortality_sub %>%
  ggplot(aes(Scenario,Score)) + 
  facet_grid(rows=vars(Category),scales="free",labeller = label_parsed)+
  geom_boxplot2(width.errorbar = 0.001, fill="blue",alpha=0.5) +
  theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size=9),
        
        text = element_text(size=15),
        strip.text.y = element_text(size = 9))
ggsave('total_sensiti_sub_2.jpeg', width = 4, height = 9,dpi=600)





# Median in %




total_sensi_tosave_percent= total_sensi_tosave %>% group_by(Category)%>% 
  mutate(Min_median =min(Median)) %>%
ungroup()




  
total_sensi_tosave_percent$Median_inc =1
for (category in unique(total_sensi_tosave_percent$Category)) {
    median_inc = total_sensi_tosave_percent[which(total_sensi_tosave_percent$Category == category & total_sensi_tosave_percent$Scenario == 'Current'),"Median"]
    #median_incline = total_sensi_tosave[which(total_sensi_tosave$Category == category | total_sensi_tosave$Scenario == 'Inc.'),]
    total_sensi_tosave_percent[which(total_sensi_tosave_percent$Category == category),"Median_inc"] =median_inc
    #print(median_inc)
    }

  
  




total_sensi_tosave_percent= total_sensi_tosave_percent %>% group_by(Category)%>% 
  mutate(`% median increase` =100*(abs(Median_inc)/Median_inc)*(Median-Median_inc)/Median_inc) %>%
  ungroup()




total_sensi_tosave_percent_without_inc = total_sensi_tosave_percent[ which(total_sensi_tosave_percent$Scenario %in% c( "Mort. hatchery",
                                                        "Mort. fry-fingerling",
                                                        "Mort. On-growing 1",
                                                        "Mort. On-growing DK",
                                                        "Mort. Seafarm 1",
                                                        "Mort. Seafarm 2")),]


# expot this as csv
write.csv(total_sensi_tosave_percent,"Mean_and_median_scores_sensitivity_2.csv",row.names = TRUE)









total_sensi_tosave_percent_without_inc %>%
  ggplot(aes(x=Scenario,y=`% median increase`,group=1)) + 
  facet_grid(rows=vars(Category),scales="free",labeller = label_parsed)+
  geom_line(color="red",linewidth=1.5,alpha=0.5) + geom_point() + theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
        text = element_text(size=15),
        strip.text.y = element_text(size = 8,angle=0))+
  ylim(-2, 17)
ggsave('total_sensiti_median_percentage_2.jpeg', width = 3.5, height = 9,dpi=600)








#Subselection




total_sensi_tosave_percent_sub= total_sensi_tosave_sub %>% group_by(Category)%>% 
  mutate(Min_median =min(Median)) %>%
  ungroup()





total_sensi_tosave_percent_sub$Median_inc =1
for (category in unique(total_sensi_tosave_percent_sub$Category)) {
  median_inc = total_sensi_tosave_percent_sub[which(total_sensi_tosave_percent_sub$Category == category & total_sensi_tosave_percent_sub$Scenario == 'Current'),"Median"]
  #median_incline = total_sensi_tosave[which(total_sensi_tosave$Category == category | total_sensi_tosave$Scenario == 'Inc.'),]
  total_sensi_tosave_percent_sub[which(total_sensi_tosave_percent_sub$Category == category),"Median_inc"] =median_inc
  #print(median_inc)
}






total_sensi_tosave_percent_sub= total_sensi_tosave_percent_sub %>% group_by(Category)%>% 
  mutate(`% median increase` =100*(abs(Median_inc)/Median_inc)*(Median-Median_inc)/Median_inc) %>%
  ungroup()




total_sensi_tosave_percent_without_inc_sub = total_sensi_tosave_percent_sub[ which(total_sensi_tosave_percent_sub$Scenario %in% c( "Mort. hatchery",
                                                                                                                       "Mort. fry-fingerling",
                                                                                                                       "Mort. On-growing 1",
                                                                                                                       "Mort. On-growing DK",
                                                                                                                       "Mort. Seafarm 1",
                                                                                                                       "Mort. Seafarm 2")),]
      







total_sensi_tosave_percent_without_inc_sub %>%
  ggplot(aes(x=Scenario,y=`% median increase`,group=1)) + 
  facet_grid(rows=vars(Category),scales="free",labeller = label_parsed)+
  geom_line(color="red",linewidth=1.5,alpha=0.5) + geom_point() + theme_bw()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
        text = element_text(size=15),
        strip.text.y = element_text(size = 11))+
  ylim(-2, 17)
ggsave('total_sensiti_median_percentage_sub_2.jpeg', width = 3.5, height = 9,dpi=600)







### PAIRED RESULTS FOR SI



sub_GW_melted %>%
  ggplot(aes(Scenario,Score, fill=Scenario)) +
  geom_boxplot(outlier.shape=NA) +
  geom_point(size=0.001)+ 
  geom_line(aes(group=paired),size=0.001,alpha=0.2) +
  theme(legend.position = "none")
ggsave('GW_paired.jpeg', width = 10, height = 7,dpi=600)








#Subselection

total_best_incworst_paired_sub = subset(total_best_incworst, !(Category %in% c("TETinf~(kg~1.4-DC-eq.)",
                                                                                        "TA100~(kg~SO[2]-eq.)",
                                                                                        "FETinf~(kg~1.4-DC-eq.)",
                                                                                        "ME~(kg~N-eq.)",
                                                                                        "FE~(kg~P-eq.)",
                                                                                        "PM~(kg~PM10-eq.)",
                                                                                        "Eutro.~(kg~N-eq.)",
                                                                                        "ODinf.~(kg~CFC-11-eq.)")))



total_best_incworst %>%
  ggplot(aes(Scenario,Score, fill=Scenario)) + 
  facet_grid(rows=vars(Category),scales="free",labeller = label_parsed)+
  geom_boxplot2(width.errorbar = 0.001,) +theme_bw()+
  geom_line(aes(group=paired),size=0.001,alpha=0.2) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size=15),
        strip.text.y = element_text(size = 10,angle=0),
        axis.text.y = element_text(size=9))
ggsave('total_dd2paired_2.jpeg', width = 7, height = 9,dpi=600)









total_best_incworst_paired_sub %>%
  ggplot(aes(Scenario,Score, fill=Scenario)) + 
  facet_grid(rows=vars(Category),scales="free",labeller = label_parsed)+
  geom_boxplot2(width.errorbar = 0.001,) +theme_bw()+
  geom_line(aes(group=paired),size=0.05,alpha=0.3) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        text = element_text(size=15),
        strip.text.y = element_text(size = 10,angle=0),
        axis.text.y = element_text(size=9))
ggsave('total_dd2paired.jpeg', width = 7, height = 12,dpi=600)



