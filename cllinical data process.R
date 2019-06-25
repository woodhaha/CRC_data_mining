#RETRIEVAL AND SUMMARY OF TCGA-COAD CLINICALAL FOLLOWUP DATA.
setwd("C:/Users/woodhaha/Desktop/CRC_data_mining")
require(FirebrowseR)
library(dplyr)
cohorts = Metadata.Cohorts(format = "csv")
cancer.Type = cohorts[grep("Colon", cohorts$description, ignore.case = T), 1]
Colon.followup = Samples.Clinical(cohort = cancer.Type, format = "tsv")
dim(Colon.followup)
all.Received = F;
page.Counter = 1;
page.size = 150
Colon.followup = list()
while (all.Received == F) {
	Colon.followup[[page.Counter]] = Samples.Clinical(format = "csv",
	cohort = cancer.Type,page_size = page.size, page = page.Counter)
	
	if (page.Counter > 1)
		colnames(Colon.followup[[page.Counter]]) = colnames(Colon.followup[[page.Counter - 1]])

	if (nrow(Colon.followup[[page.Counter]]) < page.size) {
		all.Received = T
	} else {
		page.Counter = page.Counter + 1
	}
}

Colon.followup = do.call(rbind, Colon.followup)

write.csv(Colon.followup,"TCGA-COAD.followup-fullData.csv")
clinical = Colon.followup %>% select("tcga_participant_barcode", 
				     "age_at_initial_pathologic_diagnosis", "anatomic_neoplasm_subdivision",
				     "days_to_death", "days_to_last_followup", "gender",
				    "pathologic_stage", "venous_invasion", "vital_status",
				    "preoperative_pretreatment_cea_level","perineural_invasion_present", 
				     "number_of_lymphnodes_positive_by_he","lymphatic_invasion", "person_neoplasm_cancer_status")

colnames(clinical) = c("id", "age", "location", "days_to_death", 
                      "days_to_last_followup", "gender",
			"stage", "venous_invasion", "vital_status",
			"preoperative_CEA","perineural_invasion", 
		          "lymphnodes_positive_by_HE","lymphatic_invasion", 
					 "person_neoplasm_status")

clinical$status <- ifelse(clinical$vital_status == "alive", 0, 1)
clinical$OS <- ifelse(is.na(clinical$days_to_death), clinical$days_to_last_followup, clinical$days_to_death)
clinical$DFS <- ifelse(clinical$person_neoplasm_status == "tumor free", clinical$OS, NA)
clinical <- clinical[clinical$OS > 0,] #follow up days >=0
dim(clinical)
#Disease-free survival (DFS) is a number that tells the chances of staying free of a disease or cancer after a particular treatment. 
clinical$stage= lava::trim(gsub("?stage|a|b|c", "", clinical$stage),all=T)
#
clinical$id <- gsub("-", ".", clinical$id)
#
clinical <-clinical%>% dplyr::select(-days_to_death,-days_to_last_followup) #,
table(is.na(clinical$stage));
table(clinical$person_neoplasm_status)
colnames(clinical)= gsub("_", ".", colnames(clinical))

dim(clinical)
#
library(VIM)

complete = clinical[complete.cases(clinical),]
incomplete = clinical[!complete.cases(clinical),]
## Visualising the missing patterns
library(mice)
md.pattern(clinical)
 par(mfrow = c(1,1), omi = c(0.5,0.5,0.5,0.5), mgp=c(3, 1, 0),mar=c(5,4.1,4.1,2.1), las = 1 ,family="Lato Light")

aggr(clinical, prop = F, numbers = T,combined =F)
write.csv(clinical, "tcga_coad_clinical_followup_data.csv",row.names=F)
