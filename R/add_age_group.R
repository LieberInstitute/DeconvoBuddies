

#' Add Age Groups and adjusted ages to pd
#'
#' @param pd colData for rse object
#' @param age_col name of column containing numeric age information
#'
#' @return colData appended with age_group and adjusted ages
#' @export
#'
#' @examples
#'
#'library(SummarizedExperiment)
#'pd_test <- colData(rse_bulk_test)
#'pd_test2 <- add_age_groups(pd_test)
#'


add_age_groups <- function(pd, age_col = "Age"){

  # Add age group
  age_mins <- map(age_categories,"min")
  pd_ages <- pd[[age_col]]
  pd$ageGroup <- cut(pd_ages, c(age_mins,Inf), labels = names(age_mins))

  # Add birth and prenatal cols
  pd$Prenatal <- as.integer(pd_ages < 0)
  pd$Birth <- pd_ages
  pd$Birth[pd_ages < 0] <- 0

  # Add Adjusted ages
  adj_age <- purrr::map_dfc(age_categories[c("Infant","Child","Teen","Adult")],
                           ~.get_adjust_age(ages = pd_ages, age_limits = .x))
  pd <- cbind(pd, adj_age)
  return(pd)
}

.get_adjust_age <- function(ages, age_limits){
  ages <- ages - age_limits$max
  ages[ages < 0] <- 0
  return(ages)
}
