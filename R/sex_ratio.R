#' Calculate Sex Ratio
#'
#' @param data 
#' @param by
#' @param tips 
#' @param clusters 
#' @param strata 
#' @param id 
#' @param dob 
#' @param intv 
#' @param weight 
#'
#' @return
#' @export
#'
#' @examples
calc_srb <- function(data,
                     by = NULL,
                     tips = c(0, 3),
                     sex = "b4",
                     agegr = NULL,
                     period = NULL,
                     cohort = NULL,
                     clusters=~v021,
                     strata=~v024+v025,
                     id="bidx",
                     dob="b3",
                     intv = "v008",
                     weight= "v005",
                     origin = 1900,
                     scale = 12,
                     varmethod = "none"){
  
  data$id <- data[[id]]
  data$sex <- data[[sex]]
  data$dob <- data[[dob]]
  data$intv <- data[[intv]]
  data$weights <- data[[weight]] / mean(data[[weight]])
  
  if(is.null(by))
    by <- ~1
  
  # create a model.frame with needed variables
  vars <- unique(unlist(lapply(c(by, strata, clusters), all.vars)))
  f <- formula(paste("~", paste(vars, collapse = "+"), "+", sex))
  mf <- model.frame(formula = f, data = data, na.action = na.pass,
                    id = id, weights = weights, dob = dob, intv = intv)
  
  mf$tstop <- mf$`(intv)`
  mf$tstart <- mf$`(dob)`
  mf$birth <- 1
  
  #TODO: Detect if data is BR in wide format and reshape
  
  #calculate weighted number of births by sex and etc.
  aggr <- demog_pyears(f, mf, period=period, agegr=agegr, cohort=cohort, tips=tips,
                       event="birth", weights="(weights)", origin=origin, scale=scale)$data 

  # reshape to have male births and female births in 2 columns
  aggr$pyears <- NULL
  aggr$n <- NULL
  aggr <- reshape2::dcast(aggr, paste0("... ~", sex), value.var = "event")
  
  #calculate SRB
  if (varmethod == "none") {
    aggr$sex_ratio <- aggr$male / (aggr$male + aggr$female)
  }
  
  return(aggr)
}
