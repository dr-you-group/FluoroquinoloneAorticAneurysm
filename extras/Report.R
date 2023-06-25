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
####Find Primary analysis####
cmAnalysis <- cmResultSuite$cm_analysis %>% data.frame()
colnames(cmAnalysis) <- SqlRender::snakeCaseToCamelCase(colnames(cmAnalysis))
for(i in seq(nrow(cmAnalysis))){
  a <- jsonlite::fromJSON(cmAnalysis$definition[i])
  if(a$createStudyPopArgs$riskWindowEnd==tarPrime) analysisIdPrime <- cmAnalysis$analysisId[i]
}
analysisIdPrime #analysisIdPrime = 2

####Distribution of Preference score####
returnCamelDf (targetTable = "cm_preference_score_dist",
               andromedaObject = cmResultSuite)

tcdList <- unique(preferenceScoreDist %>%
         filter(analysisId ==!!analysisIdPrime) %>%
         select(analysisId, targetId, comparatorId, databaseId))
for(i in seq(nrow(tcdList))){
  databaseId <- tcdList[i,]$databaseId
  targetId <- tcdList[i,]$targetId
  comparatorId <- tcdList[i,]$comparatorId
  analysisId <- tcdList[i,]$analysisId
  ps = preferenceScoreDist %>% filter(databaseId == !!databaseId,
                                      targetId == !!targetId,
                                      comparatorId == !!comparatorId,
                                      analysisId == !!analysisId)
  if(nrow(ps)==0) next

  if(!file.exists(file.path(resultFolder,"ps"))) dir.create(file.path(resultFolder,"ps"))
  plotPs(ps,
         targetName= targetId, #targetName should be added
         comparatorName = comparatorId, #comparatorName should be added
         showEquiposeLabel = TRUE,
         equipoiseBounds = equipoiseBounds,
         fileName = file.path(resultFolder,"ps",
                              sprintf("ps_t_%s_c_%s_a_%s_%s.tiff",targetId, comparatorId, analysisId,databaseId)))
}

####List of results with diagnostics####
#analysis_id, #target_id #comparator_id, #outcome_id #database_id, #unblind
cmResultSuite$cm_diagnostics_summary %>%
  filter(unblind==1)
