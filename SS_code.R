####################################################################################################
####   block 0: set directory
####################################################################################################
rm(list=ls())
setwd("~/Desktop/SS_code_for_share/")
require("survey")
require("openxlsx")
require("mediation")
require("Hmisc")
require("haven")
require("ggplot2")
require("gridExtra")

extract_mediation_summary <- function (x) {
  
  clp <- 100 * x$conf.level
  isLinear.y <- ((class(x$model.y)[1] %in% c("lm", "rq")) ||
                   (inherits(x$model.y, "glm") && x$model.y$family$family ==
                      "gaussian" && x$model.y$family$link == "identity") ||
                   (inherits(x$model.y, "survreg") && x$model.y$dist ==
                      "gaussian"))
  
  printone <- !x$INT && isLinear.y
  
  if (printone) {
    
    smat <- c(x$d1, x$d1.ci, x$d1.p)
    smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
    smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
    smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
    
    rownames(smat) <- c("ACME", "ADE", "Total Effect", "Prop. Mediated")
    
  } else {
    smat <- c(x$d0, x$d0.ci, x$d0.p)
    smat <- rbind(smat, c(x$d1, x$d1.ci, x$d1.p))
    smat <- rbind(smat, c(x$z0, x$z0.ci, x$z0.p))
    smat <- rbind(smat, c(x$z1, x$z1.ci, x$z1.p))
    smat <- rbind(smat, c(x$tau.coef, x$tau.ci, x$tau.p))
    smat <- rbind(smat, c(x$n0, x$n0.ci, x$n0.p))
    smat <- rbind(smat, c(x$n1, x$n1.ci, x$n1.p))
    smat <- rbind(smat, c(x$d.avg, x$d.avg.ci, x$d.avg.p))
    smat <- rbind(smat, c(x$z.avg, x$z.avg.ci, x$z.avg.p))
    smat <- rbind(smat, c(x$n.avg, x$n.avg.ci, x$n.avg.p))
    
    rownames(smat) <- c("ACME (control)", "ACME (treated)",
                        "ADE (control)", "ADE (treated)", "Total Effect",
                        "Prop. Mediated (control)", "Prop. Mediated (treated)",
                        "ACME (average)", "ADE (average)", "Prop. Mediated (average)")
    
  }
  
  colnames(smat) <- c("Estimate", paste(clp, "% CI Lower", sep = ""),
                      paste(clp, "% CI Upper", sep = ""), "p-value")
  smat
  
}

clock_list <- c("HannumAge","HorvathAge","WeidnerAge","LinAge","VidalBraloAge",
                "SkinBloodAge","ZhangAge","YangCell","PhenoAge","GrimAgeMort","HorvathTelo","GrimAge2Mort",
                "DunedinPoAm")
# clock_list <- c("GDF15Mort","B2MMort","CystatinCMort","TIMP1Mort","ADMMort","PAI1Mort","LeptinMort",
#                 "PACKYRSMort","CRPMort","logA1CMort")

#supplementary bio markers
if ("CystatinCMort" %in% clock_list){namelist <- "supplementary biomarkers"} else {namelist <- "main"}


tt_tic <- proc.time()
####################################################################################################


####################################################################################################
####    block 1: loading datasets
####################################################################################################
nhBioA <- read.csv("Raw_clean_SS_11_24_2024.csv")

# Creat death age
nhBioA$deathage <- nhBioA$permth_int/12+ nhBioA$age

# Create binary variables for subgroups
nhBioA$blackwomen <- 0
nhBioA$blackwomen[nhBioA$black==1&nhBioA$women==1] <- 1
table(nhBioA$blackwomen,useNA = "a")
nhBioA$blackmen <- 0
nhBioA$blackmen[nhBioA$black==1&nhBioA$men==1] <- 1
table(nhBioA$blackmen,useNA = "a")
nhBioA$whitewomen <- 0
nhBioA$whitewomen[nhBioA$white==1&nhBioA$women==1] <- 1
table(nhBioA$whitewomen,useNA = "a")
nhBioA$whitemen <- 0
nhBioA$whitemen[nhBioA$white==1&nhBioA$men==1] <- 1
table(nhBioA$whitemen,useNA = "a")
table(nhBioA$mexican)
nhBioA$mexicanwomen <- 0
nhBioA$mexicanwomen[nhBioA$mexican==1&nhBioA$women==1] <- 1
table(nhBioA$mexicanwomen,useNA = "a")
nhBioA$mexicanmen <- 0
nhBioA$mexicanmen[nhBioA$mexican==1&nhBioA$men==1] <- 1
table(nhBioA$mexicanmen,useNA = "a")

# Interaction smoking
nhBioA$formerSmoker_packyrslt30[!is.na(nhBioA$packyrslt30)] <- 0
nhBioA$formerSmoker_packyrslt30[nhBioA$formerSmoker==1&nhBioA$packyrslt30==1&!is.na(nhBioA$packyrslt30)] <- 1
nhBioA$formerSmoker_packyrs30[!is.na(nhBioA$packyrslt30)] <- 0
nhBioA$formerSmoker_packyrs30[nhBioA$formerSmoker==1&(nhBioA$packyrs60==1|nhBioA$packyrs30_59==1)&!is.na(nhBioA$packyrslt30)] <- 1
nhBioA$currentSmoker_packyrslt30[!is.na(nhBioA$packyrslt30)] <- 0
nhBioA$currentSmoker_packyrslt30[nhBioA$currentSmoker==1&nhBioA$packyrslt30==1&!is.na(nhBioA$packyrslt30)] <- 1
nhBioA$currentSmoker_packyrs30[!is.na(nhBioA$packyrslt30)] <- 0
nhBioA$currentSmoker_packyrs30[nhBioA$currentSmoker==1&(nhBioA$packyrs60==1|nhBioA$packyrs30_59==1)&!is.na(nhBioA$packyrslt30)] <- 1

table(paste0(nhBioA$nonSmoker,nhBioA$formerSmoker_packyrslt30,nhBioA$formerSmoker_packyrs30,
             nhBioA$currentSmoker_packyrslt30,nhBioA$currentSmoker_packyrs30))

# Social support variables
for (i in c("anyoneSS","moresupptSS","anyonefinSS",
            "ssq060","nofriendsSS","friends1_4SS","friends5SS")){
  print(i)
  print(table(nhBioA[,i],useNA = "a"))
  nhBioA[nhBioA[,i]==9999&!is.na(nhBioA[,i]) |  nhBioA[,i]==7777&!is.na(nhBioA[,i]),i] <- NA
  print(table(nhBioA[,i],useNA = "a"))
}

table(nhBioA$anyoneSS,useNA = "a")
table(nhBioA$anyoneSS,nhBioA$moresupptSS,useNA = "a")
nhBioA$noneedmoreSS[!is.na(nhBioA$anyoneSS)] <- 0
nhBioA$noneedmoreSS[!is.na(nhBioA$anyoneSS)&nhBioA$moresupptSS==0] <- 1
nhBioA$needmoreSS[!is.na(nhBioA$anyoneSS)] <- 0
nhBioA$needmoreSS[!is.na(nhBioA$anyoneSS)&nhBioA$moresupptSS==1] <- 1
nhBioA$needmoreSS[!is.na(nhBioA$anyoneSS)&is.na(nhBioA$moresupptSS)] <- 1
table(paste0(nhBioA$anyoneSS,nhBioA$noneedmoreSS,nhBioA$needmoreSS),useNA = "a")
nhBioA$needmoreSS_nooneSS[!is.na(nhBioA$anyoneSS)] <- 0
nhBioA$needmoreSS_nooneSS[!is.na(nhBioA$anyoneSS)&nhBioA$noneedmoreSS==0] <- 1
table(paste0(nhBioA$noneedmoreSS,nhBioA$needmoreSS_nooneSS),useNA = "a")


table(nhBioA$anyonefinSS,useNA = "a")
table(paste0(nhBioA$nofriendsSS,nhBioA$friends1_4SS,nhBioA$friends5SS))
range(nhBioA$ssq060,na.rm = T)
nhBioA$HavefriendSS[nhBioA$nofriendsSS==1&!is.na(nhBioA$nofriendsSS)] <- 0
nhBioA$HavefriendSS[nhBioA$nofriendsSS==0&!is.na(nhBioA$nofriendsSS)] <- 1
table(nhBioA$HavefriendSS,nhBioA$nofriendsSS,useNA = "a")
nhBioA$Have5friendSS[nhBioA$ssq060<5|nhBioA$ssd061<5] <- 0
nhBioA$Have5friendSS[nhBioA$ssq060>=5|nhBioA$ssd061>=5] <- 1
table(nhBioA$HavefriendSS,nhBioA$Have5friendSS,useNA = "a")
table(nhBioA$HavefriendSS[nhBioA$blackwomen==1],nhBioA$anyoneSS[nhBioA$blackwomen==1])
table(nhBioA$HavefriendSS[nhBioA$whitemen==1],nhBioA$anyoneSS[nhBioA$whitemen==1])
table(nhBioA$HavefriendSS[nhBioA$pir_above5==1],nhBioA$anyoneSS[nhBioA$pir_above5==1])

nhBioA$drinker[!is.na(nhBioA$abstainer)] <- 0
nhBioA$drinker[nhBioA$abstainer==0] <- 1
table(nhBioA$abstainer,nhBioA$drinker,useNA = "a")
nhBioA$NooneSS[!is.na(nhBioA$anyoneSS)] <- 0
nhBioA$NooneSS[nhBioA$anyoneSS==0] <- 1
table(nhBioA$anyoneSS,nhBioA$NooneSS,useNA = "a")
nhBioA$NoonefinSS[!is.na(nhBioA$anyonefinSS)] <- 0
nhBioA$NoonefinSS[nhBioA$anyonefinSS==0] <- 1
table(nhBioA$anyonefinSS,nhBioA$NoonefinSS,useNA = "a")
nhBioA$Havelt5friendSS[!is.na(nhBioA$Have5friendSS)] <- 0
nhBioA$Havelt5friendSS[nhBioA$Have5friendSS==0] <- 1
table(nhBioA$Have5friendSS,nhBioA$Havelt5friendSS,useNA = "a")

####################################################################################################
# Load DNA methylation data
DNAm_Meta <- read_sas("dnmepi.sas7bdat")
names(DNAm_Meta)
analysis <- merge(nhBioA,DNAm_Meta,by.x="seqn",by.y="SEQN",all = T)

# Create indicator variable for clock data
analysis$overallSample<- 0
analysis$overallSample[!is.na(analysis$HorvathAge)] <- 1
table(analysis$overallSample,useNA = "a")
tapply(analysis$WTDN4YR,analysis$overallSample,median,na.rm=T)
####################################################################################################


####################################################################################################
####    block 2: multiple imputation data
####################################################################################################
multiimpu_analysis <- read.csv("Impute_clean_SS_11_24_2024.csv")
multiimpu_analysis$unmarried <- 0
multiimpu_analysis$unmarried[multiimpu_analysis$married==0] <- 1
multiimpu_analysis$forborn <- 0
multiimpu_analysis$forborn[multiimpu_analysis$nativity==0] <- 1
multiimpu_analysis$sedentary <- NA
multiimpu_analysis$sedentary <- ifelse(multiimpu_analysis$active==0,1,0)
multiimpu_analysis$men <- NA
multiimpu_analysis$men <- ifelse(multiimpu_analysis$women==0,1,0)
multiimpu_analysis <- merge(multiimpu_analysis,analysis[,c("seqn","WTDN4YR","sdmvpsu","sdmvstra","overallSample",clock_list,
                                                           "dead","ucod_113","deathage",
                                                           "blackmen","blackwomen","whitemen","whitewomen","mexicanmen","mexicanwomen")],
                            all=T,by="seqn")
# interaction smoking
table(paste0(multiimpu_analysis$hei_quantile1,multiimpu_analysis$hei_quantile2,multiimpu_analysis$hei_quantile3,
             multiimpu_analysis$hei_quantile4,multiimpu_analysis$hei_quantile5),
      paste0(multiimpu_analysis$pir_below1,multiimpu_analysis$pir_1_2,multiimpu_analysis$pir_2_5,
             multiimpu_analysis$pir_above5))
table(paste0(multiimpu_analysis$hei_quantile1,multiimpu_analysis$hei_quantile2,multiimpu_analysis$hei_quantile3,
             multiimpu_analysis$hei_quantile4,multiimpu_analysis$hei_quantile5),
      multiimpu_analysis$overallSample,useNA = "a")
multiimpu_analysis$formerSmoker_packyrslt30[!is.na(multiimpu_analysis$packyrslt30)] <- 0
multiimpu_analysis$formerSmoker_packyrslt30[multiimpu_analysis$formerSmoker==1&multiimpu_analysis$packyrslt30==1&!is.na(multiimpu_analysis$packyrslt30)] <- 1
multiimpu_analysis$formerSmoker_packyrs30[!is.na(multiimpu_analysis$packyrslt30)] <- 0
multiimpu_analysis$formerSmoker_packyrs30[multiimpu_analysis$formerSmoker==1&(multiimpu_analysis$packyrs60==1|multiimpu_analysis$packyrs30_59==1)&!is.na(multiimpu_analysis$packyrslt30)] <- 1
multiimpu_analysis$currentSmoker_packyrslt30[!is.na(multiimpu_analysis$packyrslt30)] <- 0
multiimpu_analysis$currentSmoker_packyrslt30[multiimpu_analysis$currentSmoker==1&multiimpu_analysis$packyrslt30==1&!is.na(multiimpu_analysis$packyrslt30)] <- 1
multiimpu_analysis$currentSmoker_packyrs30[!is.na(multiimpu_analysis$packyrslt30)] <- 0
multiimpu_analysis$currentSmoker_packyrs30[multiimpu_analysis$currentSmoker==1&(multiimpu_analysis$packyrs60==1|multiimpu_analysis$packyrs30_59==1)&!is.na(multiimpu_analysis$packyrslt30)] <- 1
table(paste0(multiimpu_analysis$nonSmoker,multiimpu_analysis$formerSmoker_packyrslt30,multiimpu_analysis$formerSmoker_packyrs30_59,multiimpu_analysis$formerSmoker_packyrs60,
             multiimpu_analysis$currentSmoker_packyrslt30,multiimpu_analysis$currentSmoker_packyrs30_59,multiimpu_analysis$currentSmoker_packyrs60))

multiimpu_analysis$drinker[!is.na(multiimpu_analysis$abstainer)] <- 0
multiimpu_analysis$drinker[multiimpu_analysis$abstainer==0] <- 1
table(multiimpu_analysis$abstainer,multiimpu_analysis$drinker,useNA = "a")

analysis$nonSmoker_interref <- analysis$nonSmoker
analysis$nonSmoker_interref[is.na(analysis$formerSmoker_packyrslt30)] <- NA
table(paste0(analysis$nonSmoker_interref,analysis$formerSmoker_packyrslt30,analysis$formerSmoker_packyrs30_59,
             analysis$formerSmoker_packyrs60,analysis$currentSmoker_packyrslt30,
             analysis$currentSmoker_packyrs30_59,analysis$currentSmoker_packyrs60))
multiimpu_analysis$nonSmoker_interref <- multiimpu_analysis$nonSmoker
###################################################################################################


####################################################################################################
####    block 3: creating raw data
####################################################################################################
rawanalysis <- analysis[analysis$overallSample==1&!is.na(analysis$overallSample),]
rawanalysis <- rawanalysis[rawanalysis$age>=60&rawanalysis$age<85,]
rawanalysis <- rawanalysis[!rawanalysis$WTDN4YR==0,]
rawanalysis <- rawanalysis[rawanalysis$otherRE==0,]
dim(rawanalysis)

rawanalysis_multiimpu <- multiimpu_analysis[multiimpu_analysis$overallSample==1&!is.na(multiimpu_analysis$overallSample),]
rawanalysis_multiimpu <- rawanalysis_multiimpu[rawanalysis_multiimpu$age>=60&rawanalysis_multiimpu$age<85,]
rawanalysis_multiimpu <- rawanalysis_multiimpu[!rawanalysis_multiimpu$WTDN4YR==0,]
rawanalysis_multiimpu <- rawanalysis_multiimpu[rawanalysis_multiimpu$race %in% c(1:3),]
dim(rawanalysis_multiimpu)
####################################################################################################


####################################################################################################
####    block 4: survey data for final estimates
####################################################################################################
analysis <- analysis[!is.na(analysis$WTDN4YR),]
svyNHE <- svydesign(id = ~sdmvpsu , strata = ~sdmvstra , nest = TRUE ,
                            weights = ~WTDN4YR, data = analysis)
svyNHEanalysis <- subset(svyNHE,overallSample==1&!is.na(overallSample))
svyNHEanalysis <- subset(svyNHEanalysis,age<85)
svyNHEanalysis <- subset(svyNHEanalysis,age>=60)
svyNHEanalysis <- subset(svyNHEanalysis,otherRE==0)
dim(svyNHEanalysis$variables)


multiimpu_analysis <- multiimpu_analysis[!is.na(multiimpu_analysis$WTDN4YR),]
svyNHE_multiimpu <- svydesign(id = ~sdmvpsu , strata = ~sdmvstra , nest = TRUE ,
                                      weights = ~WTDN4YR, data = multiimpu_analysis)
svyNHEanalysis_multiimpu <- subset(svyNHE_multiimpu,overallSample==1&!is.na(overallSample))
svyNHEanalysis_multiimpu <- subset(svyNHEanalysis_multiimpu,age<85)
svyNHEanalysis_multiimpu <- subset(svyNHEanalysis_multiimpu,age>=60)
svyNHEanalysis_multiimpu <- subset(svyNHEanalysis_multiimpu,race %in% c(1:3))
dim(svyNHEanalysis_multiimpu$variables)
####################################################################################################


####################################################################################################
####    block 5: sample descriptive statistics
####################################################################################################

####################################################################################################
## block 5.1: sample descriptive for categorical variables without missing
####################################################################################################
table(rawanalysis$nonSmoker_interref ,rawanalysis$nonSmoker,useNA = "a")
variable_list_sample_categorical <- c("white","black","mexican","",
                                      "women","men","", "nativity","forborn","",
                                      "married","unmarried","",
                                      "lths","hs","somecoll","coll","",
                                      "pir_below1","pir_1_2","pir_2_5","pir_above5","",
                                      "hiwhite","lowwhite","hiblue","lowblue","nowork","",
                                      "nonSmoker","formerSmoker","currentSmoker","",
                                      "packyrs0","packyrslt30","packyrs30_59","packyrs60","",
                                      "nonSmoker_interref","formerSmoker_packyrslt30",
                                      "formerSmoker_packyrs30","currentSmoker_packyrslt30",
                                      "currentSmoker_packyrs30","",
                                      "abstainer","moddrinker","hvydrinker","",
                                      "sedentary","active","",
                                      "hei_quantile1","hei_quantile2","hei_quantile3","hei_quantile4","hei_quantile5")

# descriptive tables
i=1; descriptive_categorical <- NULL
for (var in variable_list_sample_categorical){
  explore_dataset <- get("rawanalysis") ## get completed cases dataset
  explore_surveydata <- get("svyNHEanalysis") ## get completed cases survey object

  if (var==""){
    temp_2 <- as.data.frame(matrix(nrow=1,ncol = 8))
  } else{
    print(var)
    
    temp_1 <- as.data.frame(table(explore_dataset[,var],useNA = "a"))
    temp_1$Var1 <- as.character(temp_1$Var1)
    temp_1$Var1[temp_1$Var1==1] <- "Yes"
    temp_1$Var1[temp_1$Var1==0] <- "No"
    
    temp_1$Per <- as.numeric(temp_1$Freq)/nrow(explore_dataset[!is.na(explore_dataset[,var]),])
    
    temp_1$Var1 <- factor(temp_1$Var1,levels = c("Yes"))
    
    temp_1 <- temp_1[order(temp_1$Var1),]
    temp_1 <- temp_1[!is.na(temp_1$Var1),]
    
    svytt <- as.data.frame(svytotal(design=explore_surveydata,make.formula(var),na.rm=T))
    names(svytt)[2] <- "total_se"
    svyper <- as.data.frame(svymean(design=explore_surveydata,make.formula(var),na.rm = T))
    names(svyper)[2] <- "mean_se"
    
    temp_1 <- cbind(temp_1,svytt,svyper)
    temp_2 <- cbind(c(var,rep("",nrow(temp_1)-1)),temp_1)
  }
  
  if (!is.null(names(descriptive_categorical))){
    names(temp_2) <- names(descriptive_categorical)}
  
  descriptive_categorical <- rbind(descriptive_categorical,temp_2)
  
  i=i+1
}
descriptive_categorical$Var1 <- NULL
names(descriptive_categorical) <- c("variable",paste0("sample_",c("count","percentage")),
                                    paste0("survey_",c("total","total_se","percentage","percentage_se")))
row.names(descriptive_categorical) <- NULL
wb <- createWorkbook()
addWorksheet(wb, "Table 1_nonmissing_raw")
writeData(wb,x=descriptive_categorical,sheet = "Table 1_nonmissing_raw",
          rowNames = F)
###########################################
# multiple imputed data
i=1; descriptive_categorical <- NULL
for (var in variable_list_sample_categorical){
  explore_dataset <- get("rawanalysis_multiimpu") ## get imputed dataset
  explore_surveydata <- get("svyNHEanalysis_multiimpu") ## get imputed survey object
  
  if (var==""){
    temp_2 <- as.data.frame(matrix(nrow=1,ncol = 8))
  } else if (!var %in% names(explore_dataset)){
    next
  } else{
    print(var)
    
    temp_1 <- as.data.frame(table(explore_dataset[,var],useNA = "a"))
    temp_1$Var1 <- as.character(temp_1$Var1)
    temp_1$Var1[temp_1$Var1==1] <- "Yes"
    temp_1$Var1[temp_1$Var1==0] <- "No"
    
    temp_1$Per <- as.numeric(temp_1$Freq)/nrow(explore_dataset[!is.na(explore_dataset[,var]),])
    
    temp_1$Var1 <- factor(temp_1$Var1,levels = c("Yes"))
    
    temp_1 <- temp_1[order(temp_1$Var1),]
    temp_1 <- temp_1[!is.na(temp_1$Var1),]
    
    svytt <- as.data.frame(svytotal(design=explore_surveydata,make.formula(var),na.rm=T))
    names(svytt)[2] <- "total_se"
    svyper <- as.data.frame(svymean(design=explore_surveydata,make.formula(var),na.rm = T))
    names(svyper)[2] <- "mean_se"
    
    temp_1 <- cbind(temp_1,svytt,svyper)
    temp_2 <- cbind(c(var,rep("",nrow(temp_1)-1)),temp_1)
  }
  
  if (!is.null(names(descriptive_categorical))){
    names(temp_2) <- names(descriptive_categorical)}
  
  descriptive_categorical <- rbind(descriptive_categorical,temp_2)
  
  i=i+1
}
descriptive_categorical$Var1 <- NULL
names(descriptive_categorical) <- c("variable",paste0("sample_",c("count","percentage")),
                                    paste0("survey_",c("total","total_se","percentage","percentage_se")))
row.names(descriptive_categorical) <- NULL
addWorksheet(wb, "Table 1_nonmissing_multiimpu")
writeData(wb,x=descriptive_categorical,sheet = "Table 1_nonmissing_multiimpu",
          rowNames = F)
####################################################################################################

####################################################################################################
## block 5.2: sample descriptive for continuous variables
####################################################################################################
variable_list_sample_continuous <- c("age","agesq","lbxlypct","lbxmopct","lbxnepct","lbxeopct","lbxbapct")
# descriptive tables
i=1; descriptive_continuous <- NULL
for (var in variable_list_sample_continuous){
  explore_dataset <- get("rawanalysis") ## get completed cases dataset
  explore_surveydata <- get("svyNHEanalysis") ## get completed cases survey object
  
  print(var)
  
  if (length(table(is.na(explore_dataset[,var])))>1) {
    temp_0 <- table(!is.na(explore_dataset[,var]))[1]/nrow(explore_dataset)
  } else {temp_0 <- 0}
  
  temp_1 <- mean(explore_dataset[,var],na.rm=TRUE)
  temp_2 <- sd(explore_dataset[,var],na.rm=TRUE)
  temp_3 <- quantile(explore_dataset[,var],na.rm=TRUE,
                     c(.25,.5,.75))
  temp_4 <- cbind(var,as.data.frame(temp_0),
                  as.data.frame(temp_1),as.data.frame(temp_2),
                  t(as.data.frame(temp_3)),table(is.na(explore_dataset[,var]))[1])
  
  svytt <- as.data.frame(svytotal(design=explore_surveydata,make.formula(var),na.rm=T))
  names(svytt)[2] <- "total_se"
  svyper <- as.data.frame(svymean(design=explore_surveydata,make.formula(var),na.rm = T))
  svyvariance <- as.data.frame(svyvar(design=explore_surveydata,make.formula(var),na.rm = T))
  svyiqr <- as.data.frame(svyquantile(design=explore_surveydata,make.formula(var),
                                              quantiles = c(0.25,0.5,0.75), na.rm = T)[var])
  names(svyper)[2] <- "mean_se"
  
  temp_4 <- cbind(temp_4,svytt,svyper,t(svyiqr[,1]))
  
  descriptive_continuous <- rbind(descriptive_continuous,temp_4)
  i=i+1
}
names(descriptive_continuous) <- c("variable",paste0("sample_",c("missing_rate","mean","se",
                                                                 "lower_quartile","median","upper_quartile","count")),
                                   paste0("survey_",c("total","total_se","mean","mean_se","lower_quartile",
                                                      "median","upper_quartile")))

addWorksheet(wb, "Table 1_continuous_raw")
writeData(wb,x=descriptive_continuous,sheet = "Table 1_continuous_raw",
          rowNames = F)
###########################################
# multiple imputed data
i=1; descriptive_continuous <- NULL
for (var in variable_list_sample_continuous){
  explore_dataset <- get("rawanalysis_multiimpu") ## get imputed dataset
  explore_surveydata <- get("svyNHEanalysis_multiimpu") ## get imputed survey object
  
  print(var)
  
  if (length(table(is.na(explore_dataset[,var])))>1) {
    temp_0 <- table(!is.na(explore_dataset[,var]))[1]/nrow(explore_dataset)
  } else {temp_0 <- 0}
  
  temp_1 <- mean(explore_dataset[,var],na.rm=TRUE)
  temp_2 <- sd(explore_dataset[,var],na.rm=TRUE)
  temp_3 <- quantile(explore_dataset[,var],na.rm=TRUE,
                     c(.25,.5,.75))
  temp_4 <- cbind(var,as.data.frame(temp_0),
                  as.data.frame(temp_1),as.data.frame(temp_2),
                  t(as.data.frame(temp_3)),table(is.na(explore_dataset[,var]))[1])
  
  svytt <- as.data.frame(svytotal(design=explore_surveydata,make.formula(var),na.rm=T))
  names(svytt)[2] <- "total_se"
  svyper <- as.data.frame(svymean(design=explore_surveydata,make.formula(var),na.rm = T))
  svyvariance <- as.data.frame(svyvar(design=explore_surveydata,make.formula(var),na.rm = T))
  svyiqr <- as.data.frame(svyquantile(design=explore_surveydata,make.formula(var),
                                              quantiles = c(0.25,0.5,0.75), na.rm = T)[var])
  names(svyper)[2] <- "mean_se"
  
  temp_4 <- cbind(temp_4,svytt,svyper,t(svyiqr[,1]))
  
  descriptive_continuous <- rbind(descriptive_continuous,temp_4)
  i=i+1
}
names(descriptive_continuous) <- c("variable",paste0("sample_",c("missing_rate","mean","se",
                                                                 "lower_quartile","median","upper_quartile","count")),
                                   paste0("survey_",c("total","total_se","mean","mean_se","lower_quartile",
                                                      "median","upper_quartile")))
addWorksheet(wb, "Table 1_continuous_multiimpu")
writeData(wb,x=descriptive_continuous,sheet = "Table 1_continuous_multiimpu",
          rowNames = F)
####################################################################################################
saveWorkbook(wb, file=paste0("Results/Table 1",format(Sys.Date(),"_%m_%d_%Y"), ".xlsx"), overwrite = T)
####################################################################################################


####################################################################################################
####    block 6: SS analysis models
####################################################################################################
table(rawanalysis$anyoneSS,rawanalysis$age>=60,useNA = "a")

####################################################################################################
##  block 6.0: Subgroup analysis
####################################################################################################
subgroup_SS <- c("blackmen","blackwomen","whitemen","whitewomen","mexicanmen","mexicanwomen")
SSvar_list <- c("anyoneSS","noneedmoreSS","HavefriendSS","Have5friendSS","anyonefinSS","married")

for (S_pp in subgroup_SS){
  #S_pp <- "whitewomen"；S_pp <- "mexicanmen"
  print(S_pp)
 
  explore_surveydata <- get("svyNHEanalysis_multiimpu")  ## get imputed survey object
  # explore_surveydata <- get("svyNHEanalysis") ## get completed cases survey object
  
  explore_surveydata <- subset(explore_surveydata,age>=60)
  explore_surveydata <- subset(explore_surveydata,!is.na(anyoneSS))
  dim(explore_surveydata$variables)

  explore_surveydata$variables$test <- explore_surveydata$variables[,S_pp]
  explore_surveydata <- subset(explore_surveydata,test==1)
  dim(explore_surveydata$variables)

  ####################################################################################################
  #### 6.1 association between clocks and SS variables
  for (SSvar in SSvar_list){
    #SSvar <- "HavefriendSS"
    print(SSvar)

    for (clock in clock_list){
      # clock <- "HorvathAge";clock <- "DNAmADN";clock <- "GrimAgeMort"
      # print(clock)

      #### model 0
      model_0_SS <- svyglm(design=explore_surveydata,paste0(clock,"~",SSvar),family="gaussian")
      summary(model_0_SS, df.resid=Inf)
      model_0_coef_SS <-  as.data.frame(rbind(c(summary(model_0_SS, df.resid=Inf)$df.null+1,NA,NA,NA),
                                              cbind(coef(summary(model_0_SS, df.resid=Inf))[,1],
                                                    coef(summary(model_0_SS, df.resid=Inf))[,1]-1.96*coef(summary(model_0_SS, df.resid=Inf))[,2],
                                                    coef(summary(model_0_SS, df.resid=Inf))[,1]+1.96*coef(summary(model_0_SS, df.resid=Inf))[,2],
                                                    coef(summary(model_0_SS, df.resid=Inf))[,4])))
      row.names(model_0_coef_SS)[1] <- "Raw_sample_size"
      model_0_coef_SS$label <- rownames(model_0_coef_SS)

      #### model 1
      model_1_SS <- svyglm(design=explore_surveydata,paste0(clock,"~",SSvar,"+women+age+agesq +forborn"),family="gaussian")
      summary(model_1_SS, df.resid=Inf)
      model_1_coef_SS <-  as.data.frame(rbind(c(summary(model_1_SS, df.resid=Inf)$df.null+1,NA,NA,NA),
                                              cbind(coef(summary(model_1_SS, df.resid=Inf))[,1],
                                                    coef(summary(model_1_SS, df.resid=Inf))[,1]-1.96*coef(summary(model_1_SS, df.resid=Inf))[,2],
                                                    coef(summary(model_1_SS, df.resid=Inf))[,1]+1.96*coef(summary(model_1_SS, df.resid=Inf))[,2],
                                                    coef(summary(model_1_SS, df.resid=Inf))[,4])))
      row.names(model_1_coef_SS)[1] <- "Raw_sample_size"
      model_1_coef_SS$label <- rownames(model_1_coef_SS)

      #### model 2
      model_2_SS <- svyglm(design=explore_surveydata,paste0(clock,"~",SSvar,"+women+age+agesq +forborn+
                            lths+hs+somecoll+  pir_below1+pir_1_2+pir_2_5+
                            lowwhite+hiblue+lowblue+nowork"),family="gaussian")
      summary(model_2_SS, df.resid=Inf)
      model_2_coef_SS <-  as.data.frame(rbind(c(summary(model_2_SS, df.resid=Inf)$df.null+1,NA,NA,NA),
                                              cbind(coef(summary(model_2_SS, df.resid=Inf))[,1],
                                                    coef(summary(model_2_SS, df.resid=Inf))[,1]-1.96*coef(summary(model_2_SS, df.resid=Inf))[,2],
                                                    coef(summary(model_2_SS, df.resid=Inf))[,1]+1.96*coef(summary(model_2_SS, df.resid=Inf))[,2],
                                                    coef(summary(model_2_SS, df.resid=Inf))[,4])))
      row.names(model_2_coef_SS)[1] <- "Raw_sample_size"
      model_2_coef_SS$label <- rownames(model_2_coef_SS)

      #### model 3
      model_3_SS <- svyglm(design=explore_surveydata,paste0(clock,"~",SSvar,"+women+age+agesq +forborn+
                            lths+hs+somecoll+  pir_below1+pir_1_2+pir_2_5+
                            lowwhite+hiblue+lowblue+nowork+
                            abstainer+
                            active+
                            hei_quantile1+hei_quantile2+hei_quantile3+hei_quantile4+
                            formerSmoker_packyrslt30+formerSmoker_packyrs30+
                            currentSmoker_packyrslt30+currentSmoker_packyrs30"),family="gaussian")
      summary(model_3_SS, df.resid=Inf)
      model_3_coef_SS <-  as.data.frame(rbind(c(summary(model_3_SS, df.resid=Inf)$df.null+1,NA,NA,NA),
                                              cbind(coef(summary(model_3_SS, df.resid=Inf))[,1],
                                                    coef(summary(model_3_SS, df.resid=Inf))[,1]-1.96*coef(summary(model_3_SS, df.resid=Inf))[,2],
                                                    coef(summary(model_3_SS, df.resid=Inf))[,1]+1.96*coef(summary(model_3_SS, df.resid=Inf))[,2],
                                                    coef(summary(model_3_SS, df.resid=Inf))[,4])))
      row.names(model_3_coef_SS)[1] <- "Raw_sample_size"
      model_3_coef_SS$label <- rownames(model_3_coef_SS)

      #### model 4
      model_4_SS <- svyglm(design=explore_surveydata,paste0(clock,"~",SSvar,"+women+age+agesq +forborn+
                            lths+hs+somecoll+  pir_below1+pir_1_2+pir_2_5+
                            lowwhite+hiblue+lowblue+nowork+
                            abstainer+
                            active+
                            hei_quantile1+hei_quantile2+hei_quantile3+hei_quantile4+
                            formerSmoker_packyrslt30+formerSmoker_packyrs30+
                            currentSmoker_packyrslt30+currentSmoker_packyrs30+
                            lbxlypct+lbxmopct+lbxnepct+lbxeopct+lbxbapct
                            "),family="gaussian")
      summary(model_4_SS, df.resid=Inf)
      model_4_coef_SS <-  as.data.frame(rbind(c(summary(model_4_SS, df.resid=Inf)$df.null+1,NA,NA,NA),
                                              cbind(coef(summary(model_4_SS, df.resid=Inf))[,1],
                                                    coef(summary(model_4_SS, df.resid=Inf))[,1]-1.96*coef(summary(model_4_SS, df.resid=Inf))[,2],
                                                    coef(summary(model_4_SS, df.resid=Inf))[,1]+1.96*coef(summary(model_4_SS, df.resid=Inf))[,2],
                                                    coef(summary(model_4_SS, df.resid=Inf))[,4])))
      row.names(model_4_coef_SS)[1] <- "Raw_sample_size"
      model_4_coef_SS$label <- rownames(model_4_coef_SS)

      ##### combine and print table
      model_list_SS <- grep("_coef_SS",ls(),value=T)
      modelline_SS <- nrow(model_3_coef_SS)+3
      model_printtable <-  as.data.frame(matrix(c(model_3_coef_SS$label),nrow = modelline_SS,ncol = 1))
      names(model_printtable) <- "label"
      model_printtable$order <- 1:modelline_SS

      for (model in model_list_SS){
        temp_model <- get(model)
        names(temp_model) <- c(paste0(substr(model,1,(nchar(model)-8)),c("_Beta","_lowCI","_highCI","_pvalue")),"label")

        model_printtable <- merge(model_printtable,temp_model,by="label",all=T)
      }
      model_printtable <- model_printtable[order(model_printtable$order),]

              if (clock %in% c("HannumAge", "GDF15Mort")){
          wb <- createWorkbook()
        }
        addWorksheet(wb, paste0(clock))
        writeData(wb,x=model_printtable,sheet = paste0(clock),
                  rowNames = F)
      }
      saveWorkbook(wb, file=paste0("Results/SS (",SSvar,") associated with clocks among ",S_pp," participants ",namelist,
                                   format(Sys.Date(),"_%m_%d_%Y"), ".xlsx"), overwrite = T)
    }
  ####################################################################################################
  ####  6.2: GGplots for bar plots for clocks and SS variables
  library(ggplot2)
  pdf(paste0("Results/bar plots for SS associated with clocks among ",S_pp,
             " participants ",namelist,format(Sys.Date(),"_%m_%d_%Y"),".pdf"),
      width = 12,height = 5)

  for (tt in SSvar_list){
    #tt <- "married"

    for (j in c(4)){1:5
      #j=4

      dataforplots <-NULL
      for (i in 1:length(clock_list)){
        #i=1;i=13;i=15
        # print(i)
        # print(clock_list[i])
        work <-  readxl::read_excel(paste0("Results/SS (",tt,") associated with clocks among ",S_pp," participants ",namelist,
                                           format(Sys.Date(),"_%m_%d_%Y"), ".xlsx"),
                                    sheet = i)
        SS_list <- strsplit(tt,"\\+")
        work <- work[work$label %in% unlist(SS_list),]

        #j=1
        temp <- work[,c(1,2,3+(j-1)*4,4+(j-1)*4,5+(j-1)*4,6+(j-1)*4)]
        model_title <- substr(names(temp)[3],1,nchar(names(temp)[3])-5)
        names(temp) <- c("variable","order","model_coef","coef_lowCI","coef_highCI","coef_pvalue")
        temp$clock <- clock_list[i]

        dataforplots <- rbind(dataforplots,temp)
        #names(dataforplots)
      }
      if(namelist== "main"){plot_clock_list <- c("HannumAge","HorvathAge","WeidnerAge","LinAge","VidalBraloAge","SkinBloodAge","ZhangAge",
                                                  " ","YangCell","PhenoAge","GrimAgeMort","HorvathTelo","GrimAge2Mort",
                                                  "   ","DunedinPoAm")
      dataforplots <- rbind(dataforplots,c(as.character(dataforplots[1,1]),NA,NA,NA,NA,NA," "))
      dataforplots <- rbind(dataforplots,c(as.character(dataforplots[1,1]),NA,NA,NA,NA,NA,"   "))
      } else {plot_clock_list <- clock_list}

      dataforplots <- as.data.frame(dataforplots)
      for (temp in c("model_coef","coef_lowCI","coef_highCI")){
        dataforplots[,temp] <- as.numeric(dataforplots[,temp])
        dataforplots[dataforplots$clock %in% c("HorvathTelo"),temp] <- dataforplots[dataforplots$clock %in% c("HorvathTelo"),temp] *10
        dataforplots[dataforplots$clock %in% c("YangCell","DunedinPoAm"),temp] <- dataforplots[dataforplots$clock %in% c("YangCell","DunedinPoAm"),temp] *100
      }

      dataforplots$clock <- factor(dataforplots$clock,levels=plot_clock_list)
      dataforplots$variable <- factor(dataforplots$variable,levels=unlist(SS_list))
      dataforplots$sign <- ""
      dataforplots$coef_pvalue <- as.numeric(dataforplots$coef_pvalue)
      dataforplots$sign[dataforplots$coef_pvalue<(0.05/(6*6))&!is.na(dataforplots$coef_pvalue)] <- "**"

      print(ggplot(data=dataforplots,aes(x=clock,y=model_coef))+
              #ylim(limit_y)+
              theme_classic()+ xlab(paste0("Coefficient for clocks from SS model (",model_title,") among ",S_pp," participants"))+
              geom_errorbar(aes(x=clock,y=model_coef,ymin=coef_lowCI,ymax=coef_highCI,color=variable),
                            position = position_dodge(0.3), width= 0.2)+
              geom_point(aes(x=clock,y=model_coef,color=variable),position = position_dodge(0.3))+
              scale_color_manual(values = c("black","orange","blue","red"))+ geom_hline(yintercept = 0)+
              geom_vline(xintercept = 8,color="purple",linetype="dashed")+geom_vline(xintercept = 14,color="purple",linetype="dashed")+
              labs(y="Model coefficients") + theme(axis.text.x = element_text(angle=30,hjust = 1))+
              geom_text(aes(label=sign), position=position_dodge(width=0.9),  hjust=-0.5,color="red")  )
    }
  }
  dev.off()
  ####################################################################################################
  ####  6.3: Survival mediation for SS variables and overall mortality

  for (var in SSvar_list){
    # var <- "anyoneSS";var <- "HavefriendSS"
    explore_surveydata$variables$test <- explore_surveydata$variables[,var]
    explore_surveydata_sub <- subset(explore_surveydata,!is.na(test)&WTDN4YR>0)
    print(var)

    ##############   MODELS
    model_printtable <-  NULL

    for (clock in clock_list){
      #clock <- "HannumAge"; clock <- "HorvathTelo"; clock <- "LinAge"

      cov_list <- "+age+agesq +forborn+lths+hs+somecoll+  pir_below1+pir_1_2+pir_2_5+
      lowwhite+hiblue+lowblue+nowork+drinker+active+
      hei_quantile1+hei_quantile2+hei_quantile3+hei_quantile4+
      formerSmoker_packyrslt30+formerSmoker_packyrs30+
      currentSmoker_packyrslt30+currentSmoker_packyrs30"
      #+lbxlypct+lbxmopct+lbxnepct+lbxeopct+lbxbapct
      
      print(clock)
      #### mediator model
      med.fit <- svyglm(design=explore_surveydata_sub,paste0(clock,"~",var,cov_list),family="gaussian")
      summary(med.fit)
      
      #### mediator_outcome model
      med_out.fit <- survreg(as.formula(paste0("Surv(deathage,dead)~",clock,cov_list)), dist="weibull",data=rawanalysis)
      summary(med_out.fit)
      
      #### outcome model without and with mediators
      model_unweight_nomed <- survreg(as.formula(paste0("Surv(deathage,dead)~",var,cov_list)), dist="weibull",data=rawanalysis)

      model_unweight <- survreg(as.formula(paste0("Surv(deathage,dead)~",var,"+",clock,cov_list)), dist="weibull",data=rawanalysis)
      summary(model_unweight)


      out.fit_nomed <- svysurvreg(design=explore_surveydata_sub,
                                          as.formula(paste0("Surv(deathage,dead)~",var,cov_list)), dist="weibull")

      out.fit <- svysurvreg(design=explore_surveydata_sub,
                                    as.formula(paste0("Surv(deathage,dead)~",var,"+",clock,cov_list)), dist="weibull")###*/+ for testing interaction between exposures and mediators
      summary(out.fit)

      #### causal mediation models
      tryCatch({
        med.out <- mediate(med.fit,out.fit,treat=var,mediator=clock)
        summary(med.out,robustSE = TRUE, sims = 100)
        # print(test.TMint(med.out, conf.level = .95))
        # plot(med.out)

        ###################################################################################################
        #### Create table
        temp <- as.data.frame(extract_mediation_summary(summary(med.out,robustSE = TRUE, sims = 100)))[c(5,8:10),]
        temp$label <- rownames(temp)
        temp <- rbind(names(temp),temp)
        temp <- cbind(temp,matrix(NA,nrow = nrow(temp),ncol=1))
        temp$sample <- med.out$nobs
        names(temp) <- 1:7
      },error=function(e){
        cat("Error",conditionMessage(e),"\n");temp <- as.data.frame(matrix(data = NA,nrow = 7,ncol = 4 ))
      })

      ###################################################################################################
      #### Create table
      temp <- as.data.frame(extract_mediation_summary(summary(med.out)))[c(5,8:10),]
      temp$label <- rownames(temp)
      temp <- rbind(names(temp),temp)
      temp <- cbind(temp,matrix(NA,nrow = nrow(temp),ncol=2))
      names(temp) <- 1:7
      temp[6,1:2] <- c(med.fit$df.null+1,dim(out.fit$model)[1])

      #################################
      weib_results_1 <- as.data.frame(cbind(summary(out.fit_nomed, df.resid=Inf)$table[,1],(summary(out.fit_nomed, df.resid=Inf)$table[,1]-summary(out.fit_nomed, df.resid=Inf)$table[,2]*1.96),
                                            (summary(out.fit_nomed, df.resid=Inf)$table[,1]+summary(out.fit_nomed, df.resid=Inf)$table[,2]*1.96),
                                            exp(summary(out.fit_nomed, df.resid=Inf)$table[,1]*(-1)/out.fit_nomed$scale),
                                            exp((summary(out.fit_nomed, df.resid=Inf)$table[,1]+summary(out.fit_nomed, df.resid=Inf)$table[,2]*1.96)*(-1)/out.fit_nomed$scale),
                                            exp((summary(out.fit_nomed, df.resid=Inf)$table[,1]-summary(out.fit_nomed, df.resid=Inf)$table[,2]*1.96)*(-1)/out.fit_nomed$scale)))
      weib_results_1$label <- rownames(weib_results_1)
      names(weib_results_1) <- c("estimate","estimate lower95%CI","estimate higher95%CI",
                                 "Hazard Ratio","HR Lower95%CI","HR Higher 95%CI")
      weib_results_1 <- rbind(names(weib_results_1),weib_results_1)
      names(weib_results_1) <- 1:7
      #################################
      weib_results_2 <- as.data.frame(cbind(summary(out.fit, df.resid=Inf)$table[,1],(summary(out.fit, df.resid=Inf)$table[,1]-summary(out.fit, df.resid=Inf)$table[,2]*1.96),
                                            (summary(out.fit, df.resid=Inf)$table[,1]+summary(out.fit, df.resid=Inf)$table[,2]*1.96),
                                            exp(summary(out.fit, df.resid=Inf)$table[,1]*(-1)/out.fit$scale),
                                            exp((summary(out.fit, df.resid=Inf)$table[,1]+summary(out.fit, df.resid=Inf)$table[,2]*1.96)*(-1)/out.fit$scale),
                                            exp((summary(out.fit, df.resid=Inf)$table[,1]-summary(out.fit, df.resid=Inf)$table[,2]*1.96)*(-1)/out.fit$scale)))
      weib_results_2$label <- rownames(weib_results_2)
      names(weib_results_2) <- c("estimate","estimate lower95%CI","estimate higher95%CI",
                                 "Hazard Ratio","HR Lower95%CI","HR Higher 95%CI")
      weib_results_2 <- rbind(names(weib_results_2),weib_results_2)
      names(weib_results_2) <- 1:7
      #################################
      weib_results_3 <- as.data.frame(cbind(summary(med_out.fit)$table[,1],(summary(med_out.fit)$table[,1]-summary(med_out.fit)$table[,2]*1.96),
                                            (summary(med_out.fit)$table[,1]+summary(med_out.fit)$table[,2]*1.96),
                                            exp(summary(med_out.fit)$table[,1]*(-1)/med_out.fit$scale),
                                            exp((summary(med_out.fit)$table[,1]+summary(med_out.fit)$table[,2]*1.96)*(-1)/med_out.fit$scale),
                                            exp((summary(med_out.fit)$table[,1]-summary(med_out.fit)$table[,2]*1.96)*(-1)/med_out.fit$scale)))
      weib_results_3$label <- rownames(weib_results_3)
      names(weib_results_3) <- c("estimate","estimate lower95%CI","estimate higher95%CI",
                                 "Hazard Ratio","HR Lower95%CI","HR Higher 95%CI")
      weib_results_3 <- rbind(names(weib_results_3),weib_results_3)
      names(weib_results_3) <- 1:7
      #################################
      exp_m_results <-  as.data.frame(rbind(c(summary(med.fit, df.resid=Inf)$df.null+1,NA,NA,NA),
                                   cbind(coef(summary(med.fit, df.resid=Inf))[,1],
                                         coef(summary(med.fit, df.resid=Inf))[,1]-1.96*coef(summary(med.fit, df.resid=Inf))[,2],
                                         coef(summary(med.fit, df.resid=Inf))[,1]+1.96*coef(summary(med.fit, df.resid=Inf))[,2],
                                         coef(summary(med.fit, df.resid=Inf))[,4])))
      row.names(exp_m_results)[1] <- "Raw_sample_size"
      names(exp_m_results) <- c("Beta","lowCI","highCI","pvalue")
      exp_m_results$label <- rownames(exp_m_results)
      exp_m_results <- cbind(exp_m_results,matrix(NA,nrow = nrow(exp_m_results),ncol=2))
      names(exp_m_results) <- 1:7
      #################################
      unweight_results_1 <- as.data.frame(cbind(summary(model_unweight_nomed)$table[,1],(summary(model_unweight_nomed)$table[,1]-summary(model_unweight_nomed)$table[,2]*1.96),
                                                (summary(model_unweight_nomed)$table[,1]+summary(model_unweight_nomed)$table[,2]*1.96),
                                                exp(summary(model_unweight_nomed)$table[,1]*(-1)/model_unweight_nomed$scale),
                                                exp((summary(model_unweight_nomed)$table[,1]+summary(model_unweight_nomed)$table[,2]*1.96)*(-1)/model_unweight_nomed$scale),
                                                exp((summary(model_unweight_nomed)$table[,1]-summary(model_unweight_nomed)$table[,2]*1.96)*(-1)/model_unweight_nomed$scale)))
      unweight_results_1$label <- rownames(unweight_results_1)
      names(unweight_results_1) <- c("estimate","estimate lower95%CI","estimate higher95%CI",
                                     "Hazard Ratio","HR Lower95%CI","HR Higher 95%CI")
      unweight_results_1 <- rbind(names(unweight_results_1),unweight_results_1)
      names(unweight_results_1) <- 1:7
      #################################
      unweight_results_2 <- as.data.frame(cbind(summary(model_unweight)$table[,1],(summary(model_unweight)$table[,1]-summary(model_unweight)$table[,2]*1.96),
                                                (summary(model_unweight)$table[,1]+summary(model_unweight)$table[,2]*1.96),
                                                exp(summary(model_unweight)$table[,1]*(-1)/model_unweight$scale),
                                                exp((summary(model_unweight)$table[,1]+summary(model_unweight)$table[,2]*1.96)*(-1)/model_unweight$scale),
                                                exp((summary(model_unweight)$table[,1]-summary(model_unweight)$table[,2]*1.96)*(-1)/model_unweight$scale)))
      unweight_results_2$label <- rownames(unweight_results_2)
      names(unweight_results_2) <- c("estimate","estimate lower95%CI","estimate higher95%CI",
                                     "Hazard Ratio","HR Lower95%CI","HR Higher 95%CI")
      unweight_results_2 <- rbind(names(unweight_results_2),unweight_results_2)
      names(unweight_results_2) <- 1:7

      ####################################################################################################
      printout_table <-  as.data.frame(rbind(c("Causal mediation model:",rep(NA,6)),temp,rep(NA,7),rep(NA,7),rep(NA,7),
                                             c("Weibull model:",rep(NA,6)),weib_results_2,rep(NA,7),
                                             rep(NA,7),rep(NA,7),
                                             c("Weibull model without mediator:",rep(NA,6)),weib_results_1,rep(NA,7),
                                             c("Mediator Weibull model without exposure:",rep(NA,6)),weib_results_3,rep(NA,7),
                                             c("Exposure mediator model:",rep(NA,6)),exp_m_results,rep(NA,7),
                                             rep(NA,7),rep(NA,7),
                                             c("Unweighted Weibull model with mediator:",rep(NA,6)),unweight_results_2,rep(NA,7),
                                             c("Unweighted Weibull model without mediator:",rep(NA,6)),unweight_results_1,rep(NA,7)))

      if (clock %in% c("HannumAge", "GDF15Mort")){
        wb <- createWorkbook()
      }
      addWorksheet(wb, paste0(clock))
      writeData(wb,x=printout_table,sheet = paste0(clock),colNames = F,
                rowNames = F)
    }
    saveWorkbook(wb, file=paste0("Results/survival model for SS ",var," among ",S_pp," participants ",namelist,
                                 format(Sys.Date(),"_%m_%d_%Y"), ".xlsx"), overwrite = T)
  }
  
  ####################################################################################################
  ####  6.4: GGplots for bar plots for Survival mediation
  library(ggplot2)
  pdf(paste0("Results/bar plots for mediation analysis among ",S_pp," participants ",namelist,format(Sys.Date(),"_%m_%d_%Y"),".pdf"),
      width = 12,height = 5)
  #S_pp <- "blackmen";namelist <- "main"

  for (tt in SSvar_list){
    #tt="anyoneSS"
    print(tt)

    dataforplots <-NULL
    for (i in 1:length(clock_list)){
      #i=1;i=13
      print(i)
      print(clock_list[i])
      work <-  readxl::read_excel(paste0("Results/survival model for SS ",tt," among ",S_pp," participants ",namelist,
                                         format(Sys.Date(),"_%m_%d_%Y"), ".xlsx"),sheet = i)
      work <- as.data.frame(work[5,])
      names(work) <- c("Estimate","coef_lowCI","coef_highCI","p_value")
      work$set <- clock_list[i]

      dataforplots <- rbind(dataforplots,work)
    }

    if(namelist== "main"){plot_clock_list <- c("HannumAge","HorvathAge","WeidnerAge","LinAge","VidalBraloAge","SkinBloodAge","ZhangAge",
                                                " ","YangCell","PhenoAge","GrimAgeMort","HorvathTelo","GrimAge2Mort",
                                                "   ","DunedinPoAm")
    dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA," "))
    dataforplots <- rbind(dataforplots,c(NA,NA,NA,NA,NA,NA,NA,"   "))
    } else {plot_clock_list <- clock_list}

    dataforplots <- as.data.frame(dataforplots)

    dataforplots$set <- factor(dataforplots$set,levels=plot_clock_list)

    for (temp in c("Estimate","coef_lowCI","coef_highCI")){
      dataforplots[,temp] <- as.numeric(dataforplots[,temp])*100
    }
    print(range(c(dataforplots$coef_lowCI,dataforplots$coef_highCI),na.rm = T))


    print(ggplot(data=dataforplots,aes(x=set,y=Estimate))+
            #ylim(limit_y)+
            coord_cartesian(ylim = c(-120,120))+
            theme_classic()+ xlab(paste0("Proportion of mediated effects of DNA methylation clocks and biomarkers on ",tt," among ",S_pp," participants"))+
            geom_errorbar(aes(x=set,y=Estimate,ymin=coef_lowCI,ymax=coef_highCI),
                          position = position_dodge(0.3), width= 0.2)+
            geom_point(aes(x=set,y=Estimate),position = position_dodge(0.3))+
            geom_hline(yintercept = 0)+geom_hline(yintercept = 50,linetype=3)+geom_hline(yintercept = 100,linetype=4)+
            geom_hline(yintercept = -50,linetype=3)+geom_hline(yintercept = -100,linetype=4)+
            geom_vline(xintercept = 8,color="purple",linetype="dashed")+geom_vline(xintercept = 14,color="purple",linetype="dashed")+
            labs(y="Proportion of mediated effects")+theme(axis.text.x = element_text(angle=30,hjust = 1))   )
  }
  dev.off()
  ####################################################################################################
}


####################################################################################################
##  block 6.5: GGplots for bar plots for clocks and SS variables by multiple page
####################################################################################################
for (j in c(3:5)){#1:5
  #j=4
  pdf(paste0("Results/bar plots for SS associated with clocks ","(model ",j-1,") ",
             namelist,format(Sys.Date(),"_%m_%d_%Y"),".pdf"), width = 12,height = 10)

  for (tt in SSvar_list){
    #tt <- "anyoneSS";tt <- "HavefriendSS"
    print(tt)

    for (S_pp in subgroup_SS){
      #S_pp <- "Whitewomen"
      print(S_pp)

      dataforplots <-NULL
      for (i in 1:length(clock_list)){
        #i=1;i=13;i=15
        # print(i)
        # print(clock_list[i])
        work <-  readxl::read_excel(paste0("Results/SS (",tt,") associated with clocks among ",S_pp," participants ",namelist,
                                           format(Sys.Date(),"_%m_%d_%Y"), ".xlsx"),
                                    sheet = i)
        SS_list <- strsplit(tt,"\\+")
        work <- work[work$label %in% unlist(SS_list),]

        #j=1
        temp <- work[,c(1,2,3+(j-1)*4,4+(j-1)*4,5+(j-1)*4,6+(j-1)*4)]
        model_title <- substr(names(temp)[3],1,nchar(names(temp)[3])-5)
        names(temp) <- c("variable","order","model_coef","coef_lowCI","coef_highCI","coef_pvalue")
        temp$clock <- clock_list[i]

        dataforplots <- rbind(dataforplots,temp)
        #names(dataforplots)
      }

      if(namelist== "main"){plot_clock_list <- c("HannumAge","HorvathAge","WeidnerAge","LinAge","VidalBraloAge","SkinBloodAge","ZhangAge",
                                                  " ","YangCell","PhenoAge","GrimAgeMort","HorvathTelo","GrimAge2Mort",
                                                  "   ","DunedinPoAm")
      dataforplots <- rbind(dataforplots,c(as.character(dataforplots[1,1]),NA,NA,NA,NA,NA," "))
      dataforplots <- rbind(dataforplots,c(as.character(dataforplots[1,1]),NA,NA,NA,NA,NA,"   "))
      } else {plot_clock_list <- clock_list}

      dataforplots <- as.data.frame(dataforplots)

      dataforplots$clock <- factor(dataforplots$clock,levels=plot_clock_list)
      dataforplots$variable <- factor(dataforplots$variable,levels=unlist(SS_list))
      dataforplots$sign <- ""
      dataforplots$coef_pvalue <- as.numeric(dataforplots$coef_pvalue)
      dataforplots$sign[dataforplots$coef_pvalue<(0.05/(6*6))&!is.na(dataforplots$coef_pvalue)] <- "**"
      dataforplots$sign2[dataforplots$coef_pvalue<(0.05)&!is.na(dataforplots$coef_pvalue)] <- "*"
      
      for (temp in c("model_coef","coef_lowCI","coef_highCI")){
        dataforplots[,temp] <- as.numeric(dataforplots[,temp])
        dataforplots[dataforplots$clock %in% c("HorvathTelo"),temp] <- dataforplots[dataforplots$clock %in% c("HorvathTelo"),temp] *10
        dataforplots[dataforplots$clock %in% c("YangCell","DunedinPoAm"),temp] <- dataforplots[dataforplots$clock %in% c("YangCell","DunedinPoAm"),temp] *100
      }

      if(namelist== "main"){
        if (tt %in% c("anyoneSS")){limit_y <- c(-12,15)} else
          if (tt %in% c("noneedmoreSS")){limit_y <- c(-8,10)} else
          if (tt %in% c("HavefriendSS")){limit_y <- c(-20,15)} else
            if (tt %in% c("Have5friendSS")){limit_y <- c(-8,8)} else
              if (tt %in% c("anyonefinSS")){limit_y <- c(-8,12)} else
                if (tt %in% c("married")){limit_y <- c(-7,6)} else {
                  limit_y <- c(min(dataforplots$coef_lowCI,na.rm = T),max(dataforplots$coef_highCI,na.rm=T)) }
      } else {      if (tt %in% c("anyoneSS")){limit_y <- c(-2,2.5)} else
        if (tt %in% c("HavefriendSS")){limit_y <- c(-2.5,3)} else
          if (tt %in% c("Have5friendSS")){limit_y <- c(-1.5,1.5)} else
            if (tt %in% c("anyonefinSS")){limit_y <- c(-1,2)} else
              if (tt %in% c("married")){limit_y <- c(-0.8,0.8)} else {
                limit_y <- c(min(dataforplots$coef_lowCI,na.rm = T),max(dataforplots$coef_highCI,na.rm=T)) }
      }
      
      #limit_y <- c(-20,18)

      print(c(min(dataforplots$coef_lowCI,na.rm = T),max(dataforplots$coef_highCI,na.rm=T)))

      plot <- ggplot(data=dataforplots,aes(x=clock,y=model_coef))+
       # ylim(limit_y)+
        coord_cartesian(ylim = limit_y)+
        theme_classic()+ xlab(paste0(S_pp," participants"))+
        geom_errorbar(aes(x=clock,y=model_coef,ymin=coef_lowCI,ymax=coef_highCI,color=variable),
                      position = position_dodge(0.3), width= 0.4)+
        geom_point(aes(x=clock,y=model_coef,color=variable),position = position_dodge(0.3),size=2.5)+
        scale_color_manual(values = c("black","orange","blue","red"))+ geom_hline(yintercept = 0)+
        geom_vline(xintercept = 8,color="purple",linetype="dashed")+geom_vline(xintercept = 14,color="purple",linetype="dashed")+
        labs(y="Model coefficients") + theme(axis.text.x = element_text(angle=30,hjust = 1),legend.position = "none")+
        geom_text(aes(y=coef_highCI,label=sign), vjust=-0.01, color="red",size = 10)+
        geom_text(aes(y=coef_highCI,label=sign2), vjust=-0.01, color="red",size = 10)
      assign(paste0("tempplot_",S_pp),plot)

    }

    grid.arrange(tempplot_blackwomen,tempplot_blackmen,tempplot_whitewomen,tempplot_whitemen,
                 tempplot_mexicanwomen,tempplot_mexicanmen,nrow=3,
                 top=paste0("Plots for ",tt," (",model_title,")") )
  }
  dev.off()
}

####################################################################################################
##  block 6.6: sample descriptive for categorical variables without missing
####################################################################################################
variable_list_sample_categorical <- c("dead","","anyoneSS","NooneSS","","noneedmoreSS","needmoreSS_nooneSS","",
                                      "HavefriendSS","nofriendsSS","","Have5friendSS","Havelt5friendSS","",
                                      "anyonefinSS","NoonefinSS","","married","unmarried","",
                                      "nativity","forborn","",
                                      "lths","hs","somecoll","coll","",
                                      "pir_below1","pir_1_2","pir_2_5","pir_above5","",
                                      "hiwhite","lowwhite","hiblue","lowblue","nowork","",
                                      "nonSmoker_interref","formerSmoker_packyrslt30","formerSmoker_packyrs30",
                                      "currentSmoker_packyrslt30","currentSmoker_packyrs30","",
                                      "abstainer","drinker","",
                                      "sedentary","active","",
                                      "hei_quantile1","hei_quantile2","hei_quantile3","hei_quantile4","hei_quantile5")

## descriptive tables
descriptive_categorical_bygroup <- as.data.frame(variable_list_sample_categorical)

rawanalysis_multiimpu$NooneSS <- 0
rawanalysis_multiimpu$NooneSS[rawanalysis_multiimpu$anyoneSS==0] <- 1
rawanalysis_multiimpu$needmoreSS_nooneSS[!is.na(rawanalysis_multiimpu$anyoneSS)] <- 0
rawanalysis_multiimpu$needmoreSS_nooneSS[!is.na(rawanalysis_multiimpu$anyoneSS)&rawanalysis_multiimpu$noneedmoreSS==0] <- 1
rawanalysis_multiimpu$nofriendsSS <- 0
rawanalysis_multiimpu$nofriendsSS[rawanalysis_multiimpu$HavefriendSS==0] <- 1
rawanalysis_multiimpu$Havelt5friendSS <- 0
rawanalysis_multiimpu$Havelt5friendSS[rawanalysis_multiimpu$Have5friendSS==0] <- 1
rawanalysis_multiimpu$NoonefinSS <- 0
rawanalysis_multiimpu$NoonefinSS[rawanalysis_multiimpu$anyonefinSS==0] <- 1

for (tt in c("blackwomen","blackmen","whitewomen","whitemen","mexicanwomen","mexicanmen")){
  i=1;descriptive_categorical <- NULL

  for (var in variable_list_sample_categorical){
    #var="noneedmoreSS"
    explore_dataset <- get("rawanalysis_multiimpu") ## get imputed dataset
    # explore_dataset <- get("rawanalysis") ## get completed cases dataset

    explore_dataset <- explore_dataset[explore_dataset[,tt]==1,]

    if (var==""){
      temp_2 <- as.data.frame(matrix(nrow=1,ncol = 4))
    } else{
      print(var)

      temp_1 <- as.data.frame(table(explore_dataset[,var],useNA = "a"))
      temp_1$Var1 <- as.character(temp_1$Var1)
      temp_1$Var1[temp_1$Var1==1] <- "Yes"
      temp_1$Var1[temp_1$Var1==0] <- "No"
      if (!"Yes" %in% temp_1$Var1){temp_1 <- rbind(temp_1,c("Yes",0))}

      temp_1$Per <- as.numeric(temp_1$Freq)/nrow(explore_dataset[!is.na(explore_dataset[,var]),])

      temp_1$Var1 <- factor(temp_1$Var1,levels = c("Yes"))

      temp_1 <- temp_1[order(temp_1$Var1),]
      temp_1 <- temp_1[!is.na(temp_1$Var1),]

      temp_2 <- cbind(c(var,rep("",nrow(temp_1)-1)),temp_1)
    }

    if (!is.null(names(descriptive_categorical))){
      names(temp_2) <- names(descriptive_categorical)}

    descriptive_categorical <- rbind(descriptive_categorical,temp_2)

    i=i+1
  }
  descriptive_categorical$`c(var, rep("", nrow(temp_1) - 1))` <- NULL
  descriptive_categorical$Var1 <- NULL
  names(descriptive_categorical) <- paste0( tt,"_", c("count","percentage"))


  descriptive_categorical_bygroup <- cbind(descriptive_categorical_bygroup,descriptive_categorical)
}

row.names(descriptive_categorical_bygroup) <- NULL
wb <- createWorkbook()
addWorksheet(wb, "Table 1_nonmissing_raw")
writeData(wb,x=descriptive_categorical_bygroup,sheet = "Table 1_nonmissing_raw",
          rowNames = F)

####################################################################################################
##  block 6.7: sample descriptive for continuous variables
####################################################################################################
variable_list_sample_continuous <- c("age","agesq","lbxlypct","lbxmopct","lbxnepct","lbxeopct","lbxbapct")
## descriptive tables
descriptive_continuous_bygroup <- as.data.frame(variable_list_sample_continuous)

for (tt in c("blackwomen","blackmen","whitewomen","whitemen",
             "mexicanwomen","mexicanmen")){
  i=1; descriptive_continuous <- NULL
  for (var in variable_list_sample_continuous){
    #var <- "age"
    
    explore_dataset <- get("rawanalysis_multiimpu") ## get imputed dataset
    # explore_dataset <- get("rawanalysis") ## get completed cases dataset

    explore_dataset <- explore_dataset[!is.na(explore_dataset$anyoneSS),]
    explore_dataset <- explore_dataset[explore_dataset[,tt]==1,]

    print(var)

    if (length(table(is.na(explore_dataset[,var])))>1) {
      temp_0 <- table(!is.na(explore_dataset[,var]))[1]/nrow(explore_dataset)
    } else {temp_0 <- 0}

    temp_1 <- mean(explore_dataset[,var],na.rm=TRUE)
    temp_2 <- sd(explore_dataset[,var],na.rm=TRUE)
    # temp_3 <- quantile(explore_dataset[,var],na.rm=TRUE,
    #                    c(.25,.5,.75))
    temp_4 <- cbind(var,table(is.na(explore_dataset[,var]))[1],#as.data.frame(temp_0),
                    as.data.frame(temp_1),as.data.frame(temp_2)#,
                    #t(as.data.frame(temp_3)),
    )

    descriptive_continuous <- rbind(descriptive_continuous,temp_4)
    i=i+1
  }
  descriptive_continuous$var <- NULL
  names(descriptive_continuous) <- paste0( tt,"_", c(#"missing_rate",
    "count","mean","se"#,
    #"lower_quartile","median","upper_quartile",
  ))


  descriptive_continuous_bygroup <- cbind(descriptive_continuous_bygroup,descriptive_continuous)
}

addWorksheet(wb, "Table 1_continuous")
writeData(wb,x=descriptive_continuous_bygroup,sheet = "Table 1_continuous",
          rowNames = F)
saveWorkbook(wb, file=paste0("Results/Table 1",
                             format(Sys.Date(),"_%m_%d_%Y"), ".xlsx"), overwrite = T)
####################################################################################################


####################################################################################################
##  block 6.8: printing out tables for SS
####################################################################################################
SS_tablelist <- c("anyoneSS","noneedmoreSS","HavefriendSS","Have5friendSS","anyonefinSS","married")
SS_tablelistPP <- c("blackwomen","blackmen","whitewomen","whitemen",
                      "mexicanwomen","mexicanmen")

wb <- createWorkbook()
p_adjust_namelist <- c()
p_adjust_list <- c()
for (tt in SS_tablelist){
  #tt <- "anyoneSS"
  dataforplots_final <- NULL

  for (S_pp in SS_tablelistPP){
    #S_pp <- "blackwomen"

    for (j in c(4)){
      #j=4;5
      dataforplots <-NULL

      for (i in 1:length(clock_list)){
        #i=1;i=13
        print(i)
        print(clock_list[i])
        work <-  readxl::read_excel(paste0("Results/SS (",tt,") associated with clocks among ",S_pp," participants ",namelist,
                                           format(Sys.Date(),"_%m_%d_%Y"), ".xlsx"),sheet = i)
        work <- as.data.frame(work[3,])

        p_adjust_list <- c(p_adjust_list,as.numeric(work[,6+(j-1)*4]))
        p_adjust_namelist <- c(p_adjust_namelist,paste0(tt,"_",S_pp,"_",clock_list[i]))
        temp <- work[,c(1,2,3+(j-1)*4,4+(j-1)*4,5+(j-1)*4)]
        names(temp) <- c("SSexposure","order","model_coef","coef_lowCI","coef_highCI")
        temp$clock <- clock_list[i]

        dataforplots <- rbind(dataforplots,temp)
      }
      dataforplots$set <- S_pp
      dataforplots$gap <- NA

      if (is.null(dataforplots_final)){
        dataforplots_final <- dataforplots
      } else {dataforplots_final<- cbind(dataforplots_final,dataforplots)}
    }
  }
  addWorksheet(wb, paste0(tt))
  writeData(wb,x=dataforplots_final,sheet = paste0(tt),
            rowNames = F)
}

ll <- p.adjust(p_adjust_list, method = "fdr", n = length(p_adjust_list))
p_adjust_list[ll<=0.05]
(1:1404)[ll<=0.05]
p_adjust_namelist[ll<=0.05]
hist(p_adjust_list)
#ll[p_adjust_namelist=="HavefriendSS_blackmen_WeidnerAge"]
#ll[p_adjust_namelist=="HavefriendSS_blackmen_SkinBloodAge"]
#range(ll);table(ll<=0.05)




  dataforplots_final <- NULL
  for (S_pp in SS_tablelistPP){
    #S_pp <- "blackwomen"

    dataforplots <-NULL
    for (tt in SS_tablelist){
      #tt <- "anyoneSS"
        work <-  readxl::read_excel(paste0("Results/survival model for SS ",tt," among ",S_pp," participants ",namelist,
                                           format(Sys.Date(),"_%m_%d_%Y"), ".xlsx"),sheet = 1)
        work <- as.data.frame(work[45,])
        work <- work[,c(4:6)]
        names(work) <- c("HR","HR_lowCI","HR_highCI")
        work$set <- tt
        work$pp <- S_pp

        dataforplots <- rbind(dataforplots,work)
  }
      dataforplots$gap <- NA

      if (is.null(dataforplots_final)){
        dataforplots_final <- dataforplots
      } else {dataforplots_final<- cbind(dataforplots_final,dataforplots)}
}
  addWorksheet(wb, "Survival")
  writeData(wb,x=dataforplots_final,sheet = "Survival",
            rowNames = F)
saveWorkbook(wb, file=paste0("Results/model",j," print tables.xlsx"), overwrite = T)
####################################################################################################


####################################################################################################
tt_toc <- proc.time()
tt_tic_toc_time <- tt_toc-tt_tic
tt_tic_toc_time/3600
####################################################################################################




