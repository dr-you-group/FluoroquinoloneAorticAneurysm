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
# install.packages("jsonlite")

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
dbTidy

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
