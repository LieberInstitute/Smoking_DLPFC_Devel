#Table One
setwd('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq')
#load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rawCounts_szControlDLPFC.rda') #load this file for the mapObjects 
load('/dcl01/lieber/ajaffe/Brain/DLPFC_PolyA/szControl/data/rpkmCounts_szControlDLPFC.rda') #load rpkm data to filter junctions before DESeq2 pipeline--speed up computation!
library(LIBDpheno)
library(ggplot2)
library(reshape2)
#####################
#Add in LIMS information to data on cluster which isn't otherwise available
LIMS = toxicant[[1]]
pd.model = pd[pd$Dx =="Control",]
id <- match( brnum(pd.model$BrNum), LIMS$brnumerical  )
id <- id[!is.na(id)]
pd.model <- cbind( LIMS[id,], pd.model)
pd.model$source <- droplevels(pd.model$source)
pd.model <- pd.model[,colSums(sapply(pd.model, is.na))!=nrow(pd.model)]
pd.model <- pd.model[ ,!colnames(pd.model) %in% 
						colnames(pd.model[,sapply(pd.model, is.logical)])[!sapply( pd.model[,sapply(pd.model, is.logical)], any, na.rm=TRUE)] ] #drop logicals without any trues
#seeing smoking and toxicant status

################################################################
###### Nicotine Brain Numbers
pd.model$nicotine_brain_number <- NA
pd.model$nicotine_brain_number <- ifelse( grepl('ng/g|pg/mg', pd.model$nicotine_comments), pd.model$nicotine_comments, NA)
pd.model$nicotine_brain_number <- gsub(';.*', '', pd.model$nicotine_brain_number, ignore.case = TRUE) #removing all text after ;
pd.model$nicotine_brain_number <- gsub("[^0-9|.]", "", pd.model$nicotine_brain_number)
pd.model$nicotine_brain_number <- as.numeric(pd.model$nicotine_brain_number)
#pd.model$nicotine_brain_number[which(!pd.model$nicotine & !df$cotinine & df$primarydx=='Control'& grepl('CB', df$nicotine_comments))] <- 0
######Cotinine Brain Number
pd.model$cotinine_brain_number <- NA
pd.model$cotinine_brain_number <- ifelse( grepl('ng/g|pg/mg', pd.model$nicotine_comments), pd.model$nicotine_comments, NA)
pd.model$cotinine_brain_number <- gsub('.*cotinine', '', pd.model$cotinine_brain_number, ignore.case = TRUE) #removing all text after ;
pd.model$cotinine_brain_number <- gsub('nicotine.*', '', pd.model$cotinine_brain_number, ignore.case = TRUE) #removing all text after ;
pd.model$cotinine_brain_number <- gsub("[^0-9|.]", "", pd.model$cotinine_brain_number)
pd.model$cotinine_brain_number <- as.numeric(pd.model$cotinine_brain_number)
#pd.model$nicotine_brain_number[which(!pd.model$nicotine & !df$cotinine & df$primarydx=='Control'& grepl('CB', df$nicotine_comments))] <- 0
###### Nicotine Blood Numbers
pd.model$nicotine_blood_number <- NA
pd.model$nicotine_blood_number <- ifelse(grepl('ng/mL', pd.model$nicotine_comments), pd.model$nicotine_comments, NA)
pd.model$nicotine_blood_number <- gsub('cotinine.*', '', pd.model$nicotine_blood_number, ignore.case = TRUE) #removing all text after ng/mL
pd.model$nicotine_blood_number <- gsub('.*Nicotine', '', pd.model$nicotine_blood_number, ignore.case = TRUE) #removing all text before Nicotine
pd.model$nicotine_blood_number <- gsub("[^0-9|.]", "", pd.model$nicotine_blood_number)
pd.model$nicotine_blood_number <- as.numeric(pd.model$nicotine_blood_number)
###### Cotinine Blood Concentrations
pd.model$cotinine_blood_number <- NA
pd.model$cotinine_blood_number <- ifelse(grepl('ng/mL', pd.model$nicotine_comments), pd.model$nicotine_comments, NA)
pd.model$cotinine_blood_number <- gsub('.*cotinine', '', pd.model$cotinine_blood_number, ignore.case = TRUE) #removing all text before ; (usually nicotine said first)
pd.model$cotinine_blood_number <- gsub("[^0-9|.]", "", pd.model$cotinine_blood_number)
pd.model$cotinine_blood_number <- as.numeric(pd.model$cotinine_blood_number)
################################################################

#pd.model = pd.model[pd.model$Age <0,]
pd.model[pd.model$BrNum %in% c('Br1813','Br1826','Br2047', "Br2049"),
		 c("nicotine_comments","other_drugs_comments","final_dx_comments") ]   
pd.model$modelGroup <- pd.model$cotinine| pd.model$BrNum %in% c('Br1813','Br1826','Br2047', "Br2049")
pd.model$modelGroup <- as.factor(pd.model$modelGroup)
pd.model$modelGroup <- plyr::revalue(pd.model$modelGroup, c("TRUE" = "Smoker", "FALSE" = "Non-Smoker"))
pd.model$modelGroup[which(pd.model$nicotine)] <- "Smoker"
pd.model <- pd.model[!is.na(pd.model$modelGroup),]
pd.model <- pd.model[pd.model$agedeath<0|pd.model$agedeath>16,]
pd.model$AgeGroup = factor(ifelse(pd.model$agedeath<0, "Fetal","Adult"))
pd.model<- pd.model[-which(pd.model$AgeGroup=="Adult"&pd.model$smoking & pd.model$modelGroup =="Non-Smoker"),]

##########Adult Table
#setwd('/users/ssemick/Tox_Expr/Adult_Smoking/RNA_Seq/rdas')
#save(pd.model, file="PhenotypeData_Fetal_Adult.rda")
##########
setwd('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/csvs')
library(tableone)
pd.adult <- subset(pd.model,AgeGroup=="Adult")
pd.adult$source <- droplevels(pd.adult$source)
pd.adult$cotinine_blood_number[which(!pd.adult$cotinine)] <- 0
pd.adult$nicotine_blood_number[which(!pd.adult$nicotine)] <- 0

pd.adult<- pd.adult[!pd.adult$BrNum %in%c("Br1179","Br1105","Br2267"),]#outliers on SNP MDS plot droppped
adult_smoking <- CreateTableOne(vars = c("Age","Race","Sex","source", "RIN","ph","PMI","mitoRate", "smoking", "cotinine_blood_number", "nicotine_blood_number"), strata = "modelGroup", data = pd.adult)
table1 <- print(adult_smoking, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, pDigits=100)
## Save to a CSV file
write.csv(table1, file = "Adult_Smoking_TableI_OutliersDropped.csv")

##########Fetal Table
setwd('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/csvs')
library(tableone)
pd.fetal <- subset(pd.model,AgeGroup=="Fetal")
pd.fetal$source <- droplevels(pd.fetal$source)
pd.fetal<- pd.fetal[!pd.fetal$BrNum %in%c("Br1779","Br1794"),]#outliers on SNP MDS plot droppped
pd.fetal$Age = pd.fetal$Age*52 + 40 
fetal_smoking <- CreateTableOne(vars = c("Age","Race","Sex","source", "RIN","ph","PMI","mitoRate", "smoking", "cotinine_blood_number", "nicotine_blood_number"), strata = "modelGroup", data = pd.fetal)
table1 <- print(fetal_smoking, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, pDigits=100)
## Save to a CSV file
write.csv(table1, file = "Fetal_Smoking_TableI_OutliersDropped.csv")

########## Schizo Table
setwd('/dcl01/lieber/ajaffe/lab/Nicotine/DLPFC/RNAseq/csvs')
LIMS = toxicant[[1]]
pd.model = pd[pd$Dx =="Schizo",]
pd.model = pd.model[pd.model$Age >16,]
id <- match( brnum(pd.model$BrNum), LIMS$brnumerical  )
id <- id[!is.na(id)]
pd.model <- cbind( LIMS[id,], pd.model)
pd.model$source <- droplevels(pd.model$source)
pd.model <- pd.model[,colSums(sapply(pd.model, is.na))!=nrow(pd.model)]
pd.model <- pd.model[ ,!colnames(pd.model) %in% 
						colnames(pd.model[,sapply(pd.model, is.logical)])[!sapply( pd.model[,sapply(pd.model, is.logical)], any, na.rm=TRUE)] ] #drop logicals without any trues
#seeing smoking and toxicant status
table(Smoking = pd.model$smoking, Toxicant = (pd.model$nicotine | pd.model$cotinine), useNA = 'ifany' )


pd.model$modelGroup <- NA
pd.model$modelGroup[which(pd.model$smoking&!(pd.model$cotinine| pd.model$nicotine))] <- "Smoking History Only"
pd.model$modelGroup[which(pd.model$smoking&(pd.model$cotinine| pd.model$nicotine))] <- "Smoking History and Tox"
pd.model$modelGroup[which(!pd.model$smoking&!(pd.model$cotinine| pd.model$nicotine))] <- "Non-Smoker"
pd.model <- pd.model[!is.na(pd.model$modelGroup),]
pd.model$modelGroup <- factor(pd.model$modelGroup, levels=c("Non-Smoker", "Smoking History Only", "Smoking History and Tox"))
pd.model <-pd.model[!pd.model$BrNum%in%"Br1178",] 


pd.model <-pd.model[which(pd.model$modelGroup != "Smoking History Only"),]
pd.model$source <- droplevels(pd.model$source)
pd.model$modelGroup <- droplevels(pd.model$modelGroup)


schizo_smoking <- CreateTableOne(vars = c("Age","Race","Sex","source", "RIN","ph","PMI","mitoRate", "smoking"), strata = "modelGroup", data = pd.model)
table1 <- print(schizo_smoking, exact = "stage", quote = FALSE, noSpaces = TRUE, printToggle = FALSE, pDigits=100)
## Save to a CSV file
write.csv(table1, file = "Schizo_Smoking_TableI_OutliersDropped.csv")

