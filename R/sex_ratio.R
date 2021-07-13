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
                     tips = seq(0, 50, 5),
                     sex = "b4",
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
                     varmethod = "none",
                     counts = FALSE,
                     clustcounts = FALSE){

  data$id <- data[[id]]
  data$sex <- data[[sex]]
  data$dob <- data[[dob]]
  data$intv <- data[[intv]]
  data$weights <- data[[weight]] / mean(data[[weight]])

  if(is.null(by))
    by <- ~1

  # apply cohort, tips, and scale
  if(!is.null(cohort)){

    data$cohort <- cut(data[[dob]], (cohort-origin)*scale, .epis_labels(cohort),
                       include.lowest=TRUE, right=FALSE)
    
    by <- update(by, ~. + cohort)
  }

  if(!is.null(tips)){
    data$tips <- cut(x = data[[dob]] - data[[intv]],
                     breaks = -rev(tips)*scale, 
                     labels = rev(.epis_labels(tips)),
                     include.lowest = TRUE,
                     right = FALSE)
    by <- update(by, ~. + tips)
  }

  if(!is.null(period)){
    data$period <- tcut(data[[dob]], (period-origin)*scale, .epis_labels(period))
    by <- update(by, ~. + period)
  }

  # create a model.frame with needed variables
  vars <- unique(unlist(lapply(c(by, strata, clusters), all.vars)))
  f <- formula(paste("~", paste(vars, collapse = "+"), "+", sex))
  mf <- model.frame(formula = f, data = data, na.action = na.pass,
                    id = id, weights = weights, dob = dob, intv = intv)


  #TODO: Detect if data is BR in wide format and reshape

  # tabulate births by sex
  vars <- c(vars, sex)
  aggr <- dplyr::group_by(mf, dplyr::across(vars)) %>%
    dplyr::summarise(
      births = sum(`(weights)`)
    ) %>%
    data.frame()

  # reshape to have male births and female births in 2 columns
  aggr <- reshape2::dcast(aggr, paste0("... ~", sex), value.var = "births")
  aggr$male[is.na(aggr$male)] <- 0
  aggr$female[is.na(aggr$female)] <- 0
  aggr$total <- aggr$male + aggr$female

  ## construct interaction of all factor levels that appear
  byvar <- intersect(c(all.vars(by), "period", "cohort", "tips"),
                     names(aggr))
  aggr$byf <- interaction(aggr[byvar], drop=TRUE)

  ## prediction for all factor levels that appear
  pred <- data.frame(aggr[c(byvar, "byf")])[!duplicated(aggr$byf),]
  pred <- pred[order(pred$byf), ]

  if(counts || varmethod == "none"){
    mc <- model.matrix(~-1+byf, aggr)
    clong <- aggr[c("male", "female")]
    pred[c("male", "female")] <- t(mc) %*% as.matrix(clong)
  }

  #calculate SRB
  if (varmethod == "none") {
    pred$sex_ratio <- pred$male / pred$total
    pred$byf <- NULL

    if(!counts)
      pred[c("male", "female")] <- NULL

  } else if(varmethod == "lin") {

    des <- survey::svydesign(ids=clusters, strata=strata, data=aggr, weights=~1)
    class(des) <- c("svypyears", class(des))

    ## fit model
    f <- if(length(levels(aggr$byf)) == 1)
      male ~ offset(log(total))
    else
      male ~ -1 + byf + offset(log(total))

    mod <- survey::svyglm(f, des, family=quasipoisson)

    ## prediction for all factor levels that appear
    pred$total <- 1

    srb <- predict(mod, pred, type="response", vcov=TRUE)
    v <- vcov(srb)
    dimnames(v) <- list(pred$byf, pred$byf)

    pred$srb <- as.numeric(srb)
    pred$se_srb <- sqrt(diag(v))
    pred[c("byf", "total")] <- NULL
    attr(pred, "var") <- v
  } else if(varmethod %in% c("jkn", "jk1")) {

    ## Convert to array with events and PYs for each cluster
    ## reshape2::acast is MUCH faster than stats::reshape
    male_clust <- reshape2::acast(aggr, update(clusters, byf ~ .), value.var="male")
    total_clust <- reshape2::acast(aggr, update(clusters, byf ~ .), value.var="total")

    if(varmethod == "jkn"){
      aggr$strataid <- as.integer(interaction(aggr[all.vars(strata)], drop=TRUE))
      strataid <- drop(reshape2::acast(unique(aggr[c(all.vars(clusters), "strataid")]),
                                       update(clusters,  1 ~ .), value.var="strataid"))
    } else
      strataid <- NULL

    estdf <- jackknife(male_clust, total_clust, strataid)

    pred$srb <- estdf$est
    pred$se_srb <- estdf$se
    attr(pred, "var") <- vcov(estdf)
    pred$byf <- NULL

    #not sure what this does?
    if(clustcounts){
      attr(pred, "male_clust") <- male_clust
      attr(pred, "total_clust") <- total_clust
      attr(pred, "strataid") <- strataid
    }
  } else
    stop(paste0("varmethod = \"", varmethod, "\" is not recognized."))

  rownames(pred) <- NULL

  #covert to proportion male

  return(pred)
}
