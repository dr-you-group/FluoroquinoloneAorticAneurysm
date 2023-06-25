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
sosPullTable <- function(connection,
                         resultsSchema,
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

#Get balance
getBalance <- function(connection,
                       databaseId,
                       # dataFolder,
                       targetId,
                       comparatorId,
                       analysisId,
                       outcomeId = NULL){

  sql <- "SELECT *
          FROM @results_schema.@target_table
          WHERE database_id = @database_id
          AND target_id = @target_id
          AND comparator_id = @comparator_id
          AND analysis_id = @analysis_id"
  sql <- SqlRender::render(sql,
                           results_schema = resultsSchema,
                           target_table = targetTable,
                           database_id = databaseId,
                           target_id = targetId,
                           comparator_id = comparatorId,
                           analysis_id = analysisId)
  sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
  balance <- DatabaseConnector::querySql(connection, sql)


  sql <- "SELECT {@limit_true} ? {TOP @limit} *
          FROM @results_schema.@target_table;"
  sql <- SqlRender::render(sql,
                           results_schema = resultsSchema,
                           target_table = targetTable,
                           limit_true = ifelse(limit,TRUE,FALSE),
                           limit = limit)
  sql <- SqlRender::translate(sql, targetDialect = connection@dbms)
  balance <- DatabaseConnector::querySql(connection, sql)

  colnames(balance)<-SqlRender::snakeCaseToCamelCase(colnames(balance))
  colnames(covariate)<-SqlRender::snakeCaseToCamelCase(colnames(covariate))

  if(is.null(outcomeId)){
    balance <- balance[balance$analysisId == analysisId, ]
  }else{
    balance <- balance[balance$analysisId == analysisId & balance$outcomeId == outcomeId, ]
  }

  covariate <- covariate[covariate$analysisId == analysisId,]
  balance <- merge(balance, covariate[,c("covariateId", "covariateAnalysisId", "covariateName")])
  balance <- balance[ c("covariateId",
                        "covariateName",
                        "covariateAnalysisId",
                        "targetMeanBefore",
                        "comparatorMeanBefore",
                        "stdDiffBefore",
                        "targetMeanAfter",
                        "comparatorMeanAfter",
                        "stdDiffAfter")]
  colnames(balance) <- c("covariateId",
                         "covariateName",
                         "analysisId",
                         "beforeMatchingMeanTreated",
                         "beforeMatchingMeanComparator",
                         "beforeMatchingStdDiff",
                         "afterMatchingMeanTreated",
                         "afterMatchingMeanComparator",
                         "afterMatchingStdDiff")
  balance$absBeforeMatchingStdDiff <- abs(balance$beforeMatchingStdDiff)
  balance$absAfterMatchingStdDiff <- abs(balance$afterMatchingStdDiff)
  return(balance)
}

plotCovariateBalanceScatterPlot <- function(balance,
                                            beforeLabel = "Before stratification",
                                            afterLabel = "After stratification",
                                            limits = NULL) {
  if(is.null(limits)){limits <- c(min(c(balance$absBeforeMatchingStdDiff, balance$absAfterMatchingStdDiff),
                                      na.rm = TRUE),
                                  max(c(balance$absBeforeMatchingStdDiff, balance$absAfterMatchingStdDiff),
                                      na.rm = TRUE))}
  theme <- ggplot2::element_text(colour = "#000000", size = 12)
  plot <- ggplot2::ggplot(balance, ggplot2::aes(x = absBeforeMatchingStdDiff, y = absAfterMatchingStdDiff)) +
    ggplot2::geom_point(color = rgb(0, 0, 0.8, alpha = 0.3), shape = 16, size = 2) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::scale_x_continuous(beforeLabel, limits = limits) +
    ggplot2::scale_y_continuous(afterLabel, limits = limits) +
    ggplot2::theme(text = theme)

  return(plot)
}

plotPs <- function(ps, targetName, comparatorName,
                   showEquiposeLabel = TRUE, equipoiseBounds = c(0.3,0.7),
                   fileName = NULL) {
  psOrigin <- ps
  ps <- rbind(data.frame(x = ps$preferenceScore, y = ps$targetDensity, group = targetName),
              data.frame(x = ps$preferenceScore, y = ps$comparatorDensity, group = comparatorName))
  ps$group <- factor(ps$group, levels = c(as.character(targetName), as.character(comparatorName)))
  theme <- ggplot2::element_text(colour = "#000000", size = 12, margin = ggplot2::margin(0, 0.5, 0, 0.1, "cm"))
  plot <- ggplot2::ggplot(ps,
                          ggplot2::aes(x = x, y = y, color = group, group = group, fill = group)) +
    ggplot2::geom_density(stat = "identity") +
    ggplot2::scale_fill_manual(values = c(rgb(0.8, 0, 0, alpha = 0.5),
                                          rgb(0, 0, 0.8, alpha = 0.5))) +
    ggplot2::scale_color_manual(values = c(rgb(0.8, 0, 0, alpha = 0.5),
                                           rgb(0, 0, 0.8, alpha = 0.5))) +
    ggplot2::scale_x_continuous("Preference score", limits = c(0, 1)) +
    ggplot2::scale_y_continuous("Density") +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_blank(),
                   panel.grid.minor = ggplot2::element_blank(),
                   legend.position = "top",
                   legend.text = theme,
                   axis.text = theme,
                   axis.title = theme)
  if (showEquiposeLabel) {
    labelsLeft <- c()
    labelsRight <- c()
    if (showEquiposeLabel) {
      equiIndex <- psOrigin$preferenceScore>=equipoiseBounds[1] & psOrigin$preferenceScore<=equipoiseBounds[2]
      equipoise <- mean (sum(psOrigin$targetDensity[equiIndex]), sum(psOrigin$comparatorDensity[equiIndex]))/100
      labelsRight <- c(labelsRight, sprintf("%2.1f%% is in equipoise",
                                            equipoise * 100))
    }
    if (length(labelsLeft) > 0) {
      dummy <- data.frame(text = paste(labelsLeft, collapse = "\n"))
      plot <- plot + ggplot2::geom_label(x = 0, #y = max(d$y) * 1.24,
                                         hjust = "left", vjust = "top", alpha = 0.8,
                                         ggplot2::aes(label = text), data = dummy, size = 3.5)
    }
    if (length(labelsRight) > 0) {
      dummy <- data.frame(text = paste(labelsRight, collapse = "\n"))
      plot <- plot + ggplot2::annotate("label", x = 1, y = max(ps$y) * 1,
                                       hjust = "right", vjust = "top",
                                       alpha = 0.8,
                                       label = labelsRight,
                                       #ggplot2::aes(label = labelsRight),
                                       #ggplot2::aes(label = text), data = dummy,
                                       size = 3.5)
      # plot <- plot + ggplot2::geom_label(x = 1, y = max(ps$y) * 1.24,
      #                                    hjust = "right", vjust = "top",
      #                                    alpha = 0.8,
      #                                    ggplot2::aes(label = labelsRight),
      #                                    ggplot2::aes(label = text), data = dummy,
      #                                    size = 3.5)
    }
  }
  if (!is.null(fileName))
    ggplot2::ggsave(fileName, plot, width = 5, height = 3.5,
                    dpi = 400)
  return(plot)
}

returnCamelDf <- function(targetTable, andromedaObject){
  sql <- "SELECT * FROM @target_table"
  sql <- SqlRender::render(sql,
                           target_table = targetTable)
  newData <- RSQLite::dbGetQuery(andromedaObject, sql)

  colnames(newData) <- SqlRender::snakeCaseToCamelCase(colnames(newData))
  camelCaseName <- targetTable %>%
    stringr::str_remove("^cm_") %>%
    SqlRender::snakeCaseToCamelCase()

  assign(camelCaseName, newData, envir = .GlobalEnv)

  invisible(NULL)
}

