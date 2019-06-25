library(dplyr)

summary_dat <- Clinic_lncRNA_Exprs[, colnames(Clinic_lncRNA_Exprs) %in% c("id", cox_vip_list)]

row.names(summary_dat) = summary_dat$id

meregedata = merge(risk_df, summary_dat, by = "id")
meregedata = meregedata %>% arrange(desc(log10.riskScore))
meregedata =meregedata %>% dplyr::select(-id,-status) %>% mutate(
group = factor(group, levels = c("Testing", "Training"), labels = c("Testing", "Training")),
survival.status = factor(survival.status, levels = c("Alive", "Dead", "P-value"),
 labels = c("Alive", "Dead", "P-value")))
#
str(meregedata)
colnames(meregedata)[8:15] =lookup(colnames(meregedata)[8:15], Annotation23[, -2])

#CREATE A TABLE FOR SUMMARY INFORMATIONS AND STATISTICAL TEST .
library(table1)
table1::label(meregedata$OS) <- "Overall surrivall"

table1::label(meregedata$survival.status) <- "Status"
table1::label(meregedata$group) <- "Group"

table1::units(meregedata$OS) <- "days"


source("C:/Users/woodhaha/Desktop/CRC_data_mining/Scrptis/rndr3.R")


table1(~OS + riskScore + risk +log10.riskScore+ group + `RP11-440D17.3` + `EIF3J-AS1` + `TNRC6C-AS1` + MIR210HG + `AC113189.5` + LINC00261 + `CTB-25B13.12` + AC114730.3 | survival.status, group, topclass = "Rtable1-zebra",data = meregedata, droplevels = F, render = rndr3, overall = F, render.strat = rndr.strat) #
