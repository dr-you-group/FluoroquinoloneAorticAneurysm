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
Sys.setenv(DATABASECONNECTOR_JAR_FOLDER="/home/rstudio/jdbcDrivers")

# required packages
# "SqlRender"
 #"ParallelLogger"

# OHDSI shinydb read-only credentials
connectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "postgresql",
  server = keyring::key_get("sosFqAaServer"),                
  user = keyring::key_get("sosFqAaUser"),
  password = keyring::key_get("sosFqAaPassword"))
connection <- DatabaseConnector::connect(connectionDetails)

targetDatabase <- "shinydb"
resultsSchema <- "quinoloneaa"

#Filter only cm tables
cmTableNames <- resultTables %>% 
  filter(grepl("^cm",tableName)) %>%
  pull(tableName)
cmTableNames

####Helper Functions####
#ParamsCheck
checkIsClass<- function(parameter,classes) {
  name = deparse(substitute(parameter))
  if (!inherits(x = parameter, what = classes)) {
    ParallelLogger::logError(paste0(name, ' should be of class:', classes))      
    stop(paste0(name, ' is wrong class'))
  }
  return(TRUE)
}

#Pull result from Postresql 
sosPullTable <- function(resultsSchema, 
                         targetTable,
                         limit = 0 #0 means no limit;
                         ){
  checkIsClass(limit,c('integer', 'numeric'))
  sql <- "SELECT {@limit_true} ? {TOP @limit} *
          FROM @results_schema.@target_table;"
  sql <- SqlRender::render(sql, 
                           results_schema = resultsSchema,
                           target_table = targetTable,
                           limit_true = ifelse(limit,TRUE,FALSE),
                           limit = limit)
  sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
  result <- DatabaseConnector::querySql(connection, sql)
  colnames(result) <- colnames(result) %>%
    SqlRender::snakeCaseToCamelCase() #snake to camel
  return(result)
}
############################


resultDemo<- sosPullTable(resultsSchema = resultsSchema,
             targetTable = "cm_result",
             limit = 100)
resultDemo
