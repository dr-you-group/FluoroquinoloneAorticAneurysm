# Copyright 2023 Observational Health Data Sciences and Informatics
#
# This file is part of FluoroquinoloneAorticAneurysm
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# Load necessary library and Set default environment

library(dplyr)
library(Andromeda)
library(ggplot2)
library(ggsci)
# install.packages("jsonlite")
# install.packages("ggplot2")

# To install EvidenceSynthesis
# system("sudo apt install cmake")
# devtools::install_github("ohdsi/EvidenceSynthesis")

Sys.setenv(DATABASECONNECTOR_JAR_FOLDER="/home/rstudio/jdbcDrivers")
Sys.setenv(FluoroquinoloneAorticAneurysm_output_folder="/home/rstudio/output/FluoroquinoloneAorticAneurysm_output_folder")
Sys.setenv(FluoroquinoloneAorticAneurysm_temp_folder="/home/rstudio/temp")
resultFolder <- Sys.getenv("FluoroquinoloneAorticAneurysm_output_folder")
tempFolder <- Sys.getenv("FluoroquinoloneAorticAneurysm_temp_folder")

# required packages
# "SqlRender"
#"ParallelLogger"
source("./extras/HelpersForReporting.R")

# Study setting
tarPrime <- 60 # The Time-at-risk of primary analysis is defined as 60.
equipoiseBounds <- c(0.3,0.7)
targetColorR <- 255/255
targetColorG <- 99/255
targetColorB <- 71/255

comparatorColorR1st <- 30/255
comparatorColorG1st <- 144/255
comparatorColorB1st <- 255/255

comparatorColorR2nd <- 153/255
comparatorColorG2nd <- 255/255
comparatorColorB2nd <- 51/255

# targetColor <- rgb(255/255,99/255,71/255, alpha = 0.8)
# comparatorColor1st <- rgb(30/255,144/255,255/255, alpha = 0.8)
# comparatorColor2nd <- rgb(144/255,30/255,255/255, alpha = 0.8)
#
# targetColorMiddle <- rgb(255/255,99/255,71/255, alpha = 0.5)
# comparatorMiddle1st <- rgb(30/255,144/255,255/255, alpha = 0.5)
# comparatorMiddle2nd <- rgb(144/255,30/255,255/255, alpha = 0.5)
#
# targetColorFill <- rgb(255/255,99/255,71/255, alpha = 0.3)
# comparatorColor1st <- rgb(30/255,144/255,255/255, alpha = 0.3)
# comparatorColor2nd <- rgb(144/255,30/255,255/255, alpha = 0.3)

# image(1:6, 1, as.matrix(1:6), col=c(targetColor, targetColorFill, comparatorColor1, comparatorColorFill1, comparatorColor2, comparatorColorFill2), axes=FALSE , xlab="", ylab="")

# OHDSI shinydb read-only credentials
connectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "postgresql",
  server = keyring::key_get("sosFqAaServer"),
  user = keyring::key_get("sosFqAaUser"),
  password = keyring::key_get("sosFqAaPassword"))
connection <- DatabaseConnector::connect(connectionDetails)

targetDatabase <- "shinydb"
resultsSchema <- "quinoloneaa"


# List all tables from resultsSchema
sql <- "SELECT * FROM information_schema.tables
        WHERE table_schema = '@results_schema'"
sql <- SqlRender::render(sql,
                         results_schema = resultsSchema)
sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
resultTables<- DatabaseConnector::querySql(connection, sql)
colnames(resultTables) <- colnames(resultTables) %>%
  SqlRender::snakeCaseToCamelCase() #snake to camel
resultTableNames <- resultTables %>% pull(tableName)

#Filter only cm tables
cmTableNames <- resultTables %>%
  filter(grepl("^cm",tableName)) %>%
  pull(tableName)
cmTableNames

#### Retrieve all cm results to local andromeda object ####
# cmResultSuite <- Andromeda::andromeda()
# for (targetTable in cmTableNames){
#   print(sprintf("%s table is now being loaded", targetTable))
#   sql <- "SELECT * FROM @results_schema.@target_table;"
#   sql <- SqlRender::render(sql,
#                            results_schema = resultsSchema,
#                            target_table = targetTable
#                            )
#   sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
#   DatabaseConnector::querySqlToAndromeda(connection = connection,
#                                          sql,
#                                          andromeda = cmResultSuite,
#                                          andromedaTableName = targetTable,
#                                          snakeCaseToCamelCase = FALSE,
#                                          appendToTable = FALSE)
# }
# Andromeda::saveAndromeda(cmResultSuite, file.path(tempFolder, "cm_result_suite"))
cmResultSuite <- Andromeda::loadAndromeda(file.path(tempFolder, "cm_result_suite"))

############################
####TCO, Analysis, DB list####

# Cohort Definition
cohortDefinition <- sosPullTable(connection = connection,
                                 resultsSchema = resultsSchema,
                                 targetTable = "cg_cohort_definition", #"cd_cohort"
                                 limit = 0)

db <- sosPullTable(connection = connection,
                   resultsSchema = resultsSchema,
                   targetTable = "database_meta_data", #"cd_cohort"
                   limit = 0)

cohortTidy <- cohortDefinition %>%
  select(cohortDefinitionId, cohortName)

cohortTidy$cohortName <- cohortTidy %>%
  pull(cohortName) %>%
  stringr::str_remove("^\\[.*\\] ") %>% #remove '[SOS Phenotype]' prefix
  stringr::str_remove(" -  in cohorts.*females ") %>%
  stringr::str_remove(" -  first ever occurence with at least 365 days prior observation and 1 days follow up observation, males, females") %>%
  stringr::str_remove(", occurs after 2010-01-01 and before 2019-12-31")

cohortTidy$cohortNameAbbreviation <- c("AA", "AD", "AAorAD", "FQ_General", "CPH", "UTI_General", "TMP_General", "FQ", "CPH", "TMP", "UTI")

# TCO
returnCamelDf (targetTable = "cm_target_comparator_outcome",
               andromedaObject = cmResultSuite)
tco <- targetComparatorOutcome %>%
  filter (outcomeOfInterest == 1) %>%
  select (targetId, comparatorId, outcomeId)
tco$isPrimaryTco <- ifelse(tco$outcomeId == 1782489, 1, 0)

tcoTidy <- tco %>% inner_join(cohortTidy,
                              by = c("targetId" = "cohortDefinitionId")) %>%
  rename(targetName = cohortName, targetAbbreviation = cohortNameAbbreviation) %>%
  inner_join(cohortTidy,
             by = c("comparatorId" = "cohortDefinitionId")) %>%
  rename(comparatorName = cohortName, comparatorAbbreviation = cohortNameAbbreviation) %>%
  inner_join(cohortTidy,
             by = c("outcomeId" = "cohortDefinitionId")) %>%
  rename(outcomeName = cohortName, outcomeAbbreviation = cohortNameAbbreviation)

#Assin colors to target and comparator
tcoTidy$targetColorR <- targetColorR
tcoTidy$targetColorG <- targetColorG
tcoTidy$targetColorB <- targetColorB
tcoTidy$comparatorColorR <- ifelse(tcoTidy$comparatorAbbreviation=="TMP", comparatorColorR1st, tcoTidy$tcOrder <- ifelse(tcoTidy$comparatorAbbreviation=="CPH", comparatorColorR2nd, NA))
tcoTidy$comparatorColorG <- ifelse(tcoTidy$comparatorAbbreviation=="TMP", comparatorColorG1st, tcoTidy$tcOrder <- ifelse(tcoTidy$comparatorAbbreviation=="CPH", comparatorColorG2nd, NA))
tcoTidy$comparatorColorB <- ifelse(tcoTidy$comparatorAbbreviation=="TMP", comparatorColorB1st, tcoTidy$tcOrder <- ifelse(tcoTidy$comparatorAbbreviation=="CPH", comparatorColorB2nd, NA))

# DB list
dbs <- sosPullTable(connection = connection,
                    resultsSchema = resultsSchema,
                    targetTable = "database_meta_data",
                    limit = 0)
dbTidy <- dbs %>% select(cdmSourceName, cdmSourceAbbreviation, databaseId)

# Find Primary analysis
returnCamelDf (targetTable = "cm_analysis",
               andromedaObject = cmResultSuite)

analysisTidy <- analysis
analysisTidy$isPrimaryAnalysis <- 0

for(i in seq(nrow(analysisTidy))){
  a <- jsonlite::fromJSON(analysisTidy$definition[i])
  if(a$createStudyPopArgs$riskWindowEnd==tarPrime) analysisTidy$isPrimaryAnalysis[i] <- 1
}
analysisTidy <- analysisTidy %>%
  select(analysisId, description, isPrimaryAnalysis)
analysisTidy <- analysisTidy %>% rename(descriptionAnalysis = description)

analysisIdPrime = analysisTidy$analysisId[analysisTidy$isPrimaryAnalysis==1]#analysisIdPrime = c(2)

# Results of Diagnostics
returnCamelDf (targetTable = "cm_diagnostics_summary",
               andromedaObject = cmResultSuite)

resultList <- diagnosticsSummary %>%
  left_join(dbTidy,
            by = "databaseId") %>%
  left_join(analysisTidy,
            by = "analysisId") %>%
  left_join(tcoTidy,
            by = c("targetId", "comparatorId", "outcomeId")
  )

####Distribution of Preference score####
returnCamelDf (targetTable = "cm_preference_score_dist",
               andromedaObject = cmResultSuite)

tcdList <- resultList %>% select(analysisId, targetId, comparatorId, databaseId,
                                 equipoise, equipoiseDiagnostic, # unblind,
                                 cdmSourceName, cdmSourceAbbreviation,
                                 descriptionAnalysis, isPrimaryAnalysis,
                                 targetName, targetAbbreviation,
                                 comparatorName, comparatorAbbreviation,
                                 targetColorR, targetColorG, targetColorB,
                                 comparatorColorR, comparatorColorG, comparatorColorB) %>%
  unique()

# preferenceScoreDistWithDescription <- preferenceScoreDist %>%
#   left_join(tcdList,
#             by = c("analysisId", "targetId", "comparatorId", "databaseId"))
# nrow(preferenceScoreDist) == nrow(preferenceScoreDistWithDescription)

for(i in seq(nrow(tcdList))){
  if(!tcdList$isPrimaryAnalysis[i]) next #Only plot PS distribution for primary analysis

  databaseId <- tcdList[i,]$databaseId
  targetId <- tcdList[i,]$targetId
  comparatorId <- tcdList[i,]$comparatorId
  analysisId <- tcdList[i,]$analysisId
  targetColorR <- tcdList[i,]$targetColorR
  targetColorG <- tcdList[i,]$targetColorG
  targetColorB <- tcdList[i,]$targetColorB
  comparatorColorR <- tcdList[i,]$comparatorColorR
  comparatorColorG <- tcdList[i,]$comparatorColorG
  comparatorColorB <- tcdList[i,]$comparatorColorB

  ps = preferenceScoreDist %>% filter(databaseId == !!databaseId,
                                      targetId == !!targetId,
                                      comparatorId == !!comparatorId,
                                      analysisId == !!analysisId)
  targetName <- tcdList[i,]$targetAbbreviation
  comparatorName <- tcdList[i,]$comparatorAbbreviation
  dbName <- tcdList[i,]$cdmSourceAbbreviation

  if(nrow(ps)==0) next

  if(!file.exists(file.path(resultFolder,"ps"))) dir.create(file.path(resultFolder,"ps"))
  plotPs(ps,
         targetName= targetName,
         comparatorName = comparatorName,
         showEquiposeLabel = TRUE,
         equipoiseBounds = equipoiseBounds,
         fileName = file.path(resultFolder,"ps",
                              sprintf("ps_t_%s_c_%s_a_%s_%s.tiff",targetName, comparatorName, analysisId, dbName)),
         targetColorR = targetColorR,
         targetColorG = targetColorG,
         targetColorB = targetColorB,
         comparatorColorR = comparatorColorR,
         comparatorColorG = comparatorColorG,
         comparatorColorB = comparatorColorB
  )
}

####Survival analysis####
returnCamelDf (targetTable = "cm_kaplan_meier_dist",
               andromedaObject = cmResultSuite)

#Filter resluts for survival plots
resultListForSurvival <- resultList %>%
  filter(isPrimaryAnalysis == 1) %>%
  filter(isPrimaryTco == 1) %>%
  filter(unblind == 1) #only results passing the diagnostics

yLimUpperBound =
  1- kaplanMeierDist %>%
  filter(analysisId %in% unique(resultListForSurvival$analysisId)) %>%
  filter(outcomeId %in% unique(resultListForSurvival$outcomeId)) %>%
  summarise(yLimUpperBound = min(targetSurvivalLb, comparatorSurvivalLb))
yLimUpperBound = round(yLimUpperBound*1.3, digits = 3)

for (i in seq(nrow(resultListForSurvival))){

  analysisId = resultListForSurvival$analysisId[i]
  targetId = resultListForSurvival$targetId[i]
  comparatorId = resultListForSurvival$comparatorId[i]
  targetName = resultListForSurvival$targetAbbreviation[i]
  comparatorName = resultListForSurvival$comparatorAbbreviation[i]
  outcomeId = resultListForSurvival$outcomeId[i]
  outcomeName = resultListForSurvival$outcomeAbbreviation[i]
  databaseId <- resultListForSurvival$databaseId[i]
  databaseName <- resultListForSurvival$cdmSourceAbbreviation[i]

  targetColorR <- resultListForSurvival$targetColorR[i]
  targetColorG <- resultListForSurvival$targetColorG[i]
  targetColorB <- resultListForSurvival$targetColorB[i]
  comparatorColorR <- resultListForSurvival$comparatorColorR[i]
  comparatorColorG <- resultListForSurvival$comparatorColorG[i]
  comparatorColorB <- resultListForSurvival$comparatorColorB[i]

  kaplanMeier <- kaplanMeierDist[kaplanMeierDist$databaseId==databaseId&
                                   kaplanMeierDist$analysisId==analysisId&
                                   kaplanMeierDist$outcomeId==outcomeId&
                                   kaplanMeierDist$targetId == targetId&
                                   kaplanMeierDist$comparatorId == comparatorId
                                 ,]
  if(nrow(kaplanMeier)==0) next
  kpResult <- result[result$databaseId==databaseId&
                       result$analysisId==analysisId&
                       result$outcomeId==outcomeId&
                       result$targetId == targetId&
                       result$comparatorId == comparatorId
                     ,]

  kaplanMeier$targetSurvival <- 1-kaplanMeier$targetSurvival
  kaplanMeier$targetSurvivalLb <-1-kaplanMeier$targetSurvivalLb
  kaplanMeier$targetSurvivalUb <-1-kaplanMeier$targetSurvivalUb
  kaplanMeier$comparatorSurvival <-1-kaplanMeier$comparatorSurvival
  kaplanMeier$comparatorSurvivalLb <-1-kaplanMeier$comparatorSurvivalLb
  kaplanMeier$comparatorSurvivalUb <-1-kaplanMeier$comparatorSurvivalUb
  if(is.na(unique(kpResult$p))) next
  if(length(unique(kpResult$p))>1) next

  if(unique(kpResult$p)<0.001){
    pValue =  sprintf("italic(P) < 0.001")
    if (journalTheme=="jama") pValue = sprintf("italic(P)<.001")
    pValue = sprintf("italic(P)<%s",".001")
  }else {
    #if (journalTheme=="jama") pNum <- sub("^(-?)0.", "\\1.", sprintf("%.3f",unique(kpResult$p)))
    pValue =  sprintf("italic(P)==%#.3f",
                      unique(kpResult$p))

  }


  p <-plotKaplanMeier(kaplanMeier,
                      targetName,
                      comparatorName,
                      ylims = c(0,0.003),#c(0,round(max(kaplanMeier$comparatorSurvivalUb,kaplanMeier$targetSurvivalUb),2)+0.02),
                      xBreaks = NULL,#seq(from = 0,to = 2000, by = 250),#c(0,100,200,300),
                      targetColorR = targetColorR,
                      targetColorG = targetColorG,
                      targetColorB = targetColorB,
                      comparatorColorR = comparatorColorR,
                      comparatorColorG = comparatorColorG,
                      comparatorColorB = comparatorColorB,
                      pValue = pValue,
                      title = sprintf("%s vs %s for %s in %s ", targetName, comparatorName, outcomeName, databaseName))

  if(!file.exists(file.path(resultFolder,"kmplot"))) dir.create(file.path(resultFolder,"kmplot"))
  ggplot2::ggsave(file.path(resultFolder,"kmplot",sprintf("km_plot_%s_t%d_c%d_o%d_a%d.eps",
                                                          databaseId,
                                                          targetId,
                                                          comparatorId,
                                                          outcomeId,
                                                          analysisId)), p, device = "eps", width = 16, height = 12, units = "cm", dpi = 400)
  ggplot2::ggsave(file.path(resultFolder,"kmplot",sprintf("km_plot_%s_t%d_c%d_o%d_a%d.pdf",
                                                          databaseId,
                                                          targetId,
                                                          comparatorId,
                                                          outcomeId,
                                                          analysisId)), p, device = "pdf", width = 16, height = 12, units = "cm", dpi = 400)
  ggplot2::ggsave(file.path(resultFolder,"kmplot",sprintf("km_plot_%s_t%d_c%d_o%d_a%d.tiff",
                                                          databaseId,
                                                          targetId,
                                                          comparatorId,
                                                          outcomeId,
                                                          analysisId)), p, device = "tiff", width = 16, height = 12, units = "cm", dpi = 400)
  #}
}

####Meta-analysis####
returnCamelDf (targetTable = "cm_likelihood_profile",
               andromedaObject = cmResultSuite)

# returnCamelDf (targetTable = "cm_result",
#                andromedaObject = cmResultSuite)

listPrimaryAll <- resultList %>%
  filter(isPrimaryAnalysis == 1) %>%
  filter(isPrimaryTco == 1)

analysisPrimaryAll <- listPrimaryAll %>%
  select(targetId, comparatorId, outcomeId, analysisId,
         targetAbbreviation, comparatorAbbreviation, outcomeAbbreviation) %>%
  unique()

for (i in seq(nrow(analysisPrimaryAll))){

  analysisId = analysisPrimaryAll$analysisId[i]
  targetId = analysisPrimaryAll$targetId[i]
  comparatorId = analysisPrimaryAll$comparatorId[i]
  outcomeId = analysisPrimaryAll$outcomeId[i]
  targetName = analysisPrimaryAll$targetAbbreviation[i]
  comparatorName = analysisPrimaryAll$comparatorAbbreviation[i]
  outcomeName = analysisPrimaryAll$outcomeAbbreviation[i]

  sampleLP <- likelihoodProfile %>% filter(targetId == !!targetId,
                                           comparatorId == !!comparatorId,
                                           outcomeId == !!outcomeId,
                                           analysisId == !!analysisId
  )
  approx <- sampleLP %>% select("logRr", "logLikelihood", "databaseId")
  colnames(approx)[1:2] <- c("point", "value")
  approximations <- splitTable(approx, 'databaseId')
  estimate <- EvidenceSynthesis::computeBayesianMetaAnalysis(approximations)

  if(!file.exists(file.path(resultFolder,"meta"))) dir.create(file.path(resultFolder,"meta"))

  EvidenceSynthesis::plotMetaAnalysisForest(
    data = approximations,
    labels = dbs$cdmSourceAbbreviation[match(names(approximations), dbs$databaseId)],
    estimate = estimate,
    xLabel = "HR",
    summaryLabel = sprintf("Summary\n(%s vs %s)", targetName, comparatorName),
    showLikelihood = TRUE,
    fileName  = file.path(resultFolder, "meta", sprintf("t%sc%so%sa%s_all.png", targetName, comparatorName, outcomeId, analysisId))
  )
}

#primary and those passing the diagnostics
for (i in seq(analysisPrimaryAll)){

  analysisId = analysisPrimaryAll$analysisId[i]
  targetId = analysisPrimaryAll$targetId[i]
  comparatorId = analysisPrimaryAll$comparatorId[i]
  outcomeId = analysisPrimaryAll$outcomeId[i]
  targetName = analysisPrimaryAll$targetAbbreviation[i]
  comparatorName = analysisPrimaryAll$comparatorAbbreviation[i]
  outcomeName = analysisPrimaryAll$outcomeAbbreviation[i]

  databaseIds <- resultList %>%
    filter(targetId == !!targetId,
           comparatorId == !!comparatorId,
           outcomeId == !!outcomeId,
           analysisId == !!analysisId) %>%
    filter(unblind == 1) %>%
    pull(databaseId)

  sampleLP <- likelihoodProfile %>%
    filter(targetId == !!targetId,
           comparatorId == !!comparatorId,
           outcomeId == !!outcomeId,
           analysisId == !!analysisId
    ) %>%
    filter (databaseId %in% !!databaseIds) #only those passing diagnostics

  approx <- sampleLP %>% select("logRr", "logLikelihood", "databaseId")
  colnames(approx)[1:2] <- c("point", "value")
  approximations <- splitTable(approx, 'databaseId')
  estimate <- EvidenceSynthesis::computeBayesianMetaAnalysis(approximations)

  if(!file.exists(file.path(resultFolder,"meta"))) dir.create(file.path(resultFolder,"meta"))

  EvidenceSynthesis::plotMetaAnalysisForest(
    data = approximations,
    labels = dbs$cdmSourceAbbreviation[match(names(approximations), dbs$databaseId)],
    estimate = estimate,
    xLabel = "Hazard Ratio",
    summaryLabel = sprintf("Summary\n(%s vs %s)", targetName, comparatorName),
    showLikelihood = TRUE,
    fileName  = file.path(resultFolder, "meta", sprintf("t%sc%so%sa%s_passing_diag.png", targetName, comparatorName, outcomeId, analysisId))
  )
}
