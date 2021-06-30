#' @title Calculate mean children ever born (ceb)
#'
#' @description Calculate mean children ever born or other related statistics
#'   (like mean number of living children) from a DHS individual recode file.
#'
#' @param data \[`data.frame()`\]\cr
#'   Individual recode dataset.
#' @param by \[`formula()`\]\cr
#'   A formula specifying factor variables by which to stratify analysis.
#'   Default is 'NULL' to specify no stratification. Passed to [`survey::svyby()`].
#' @param agegr \[`numeric()`\]\cr
#'   The start of each age group (in years) to group each individual's age by
#'   and to stratify analysis by. Can be 'NULL' to specify no age-groups should
#'   be created and stratified by. Default is 15, 20, 25, 30, 35, 40, 45, 50.
#' @param clusters \[`formula()`\]\cr
#'   A formula specifying cluster ids to be passed as `ids` to survey::svydesign.
#'   Must correspond to variables in `data`. Default is '~v021'.
#' @param strata \[`formula()`\]\cr
#'   A formula specifying strata to be passed as `strata` to survey::svydesign.
#'   Must correspond to variables in `data`. Default is '~v024+v025'.
#' @param weight \[`character(1)`\]\cr
#'   Variable name for sample weight. Default is 'v005'.
#' @param ceb \[`character(1)`\]\cr
#'   Variable name for children ever born related variable to calculate the mean
#'   & standard error of. Default is 'v201'.
#' @param dob \[`character(1)`\]\cr
#'   Variable name for date of birth of each interviewee. Default is 'v008'.
#'   Required when `agegr` is not 'NULL'.
#' @param intv \[`character(1)`\]\cr
#'   Variable name for interview date of each interviewee. Default is 'v008'.
#'   Required when `agegr` is not 'NULL'.
#' @param varmethod \[`character(1)`\]\cr
#'   Method for variance calculation. Default is 'lin' for Taylor linearisation.
#'   Other options include 'JKn' for stratified designs & 'JK1' for stratified
#'   designs. See `survey::as.svrepdesign` for more information.
#' @param scale \[`numeric(1)`\]\cr
#'   Scale for dates inputs (`dob`, `intv`) to calendar years. 12 for CMC inputs.
#'
#' @seealso [survey::svydesign()], [survey::as.svrepdesign], [survey::svyby()],
#' & [survey::svymean()] which this function uses internally.
#'
#' @return \[`data.frame()`\] consisting of estimates and standard errors.
#' The full covariance matrix of the estimates can be retrieved by `vcov(val)`.
#'
#' @examples
#' data(zzir)
#'
#' # children ever born
#' calc_ceb(zzir)
#' # children ever born and living
#' calc_ceb(zzir, ceb = "v218")
#'
#' @importFrom dplyr %>%
#' @export
calc_ceb <- function(data,
                     by = NULL,
                     agegr = 3:10*5,
                     clusters = ~v021,
                     strata = ~v024+v025,
                     weight = "v005",
                     awfactt = NULL,
                     ceb = "v201",
                     dob = "v011",
                     intv = "v008",
                     varmethod = "lin",
                     scale = 12,
                     awfactt_scale = 100,
                     counts=FALSE,
                     clustcounts = FALSE) {

  checkmate::assert_data_frame(data)
  checkmate::assert_formula(by, null.ok = TRUE)
  checkmate::assert_formula(clusters)
  checkmate::assert_formula(strata)
  checkmate::assert_string(weight)
  checkmate::assert_string(awfactt, null.ok = TRUE)
  checkmate::assert_string(ceb)
  checkmate::assert_string(varmethod) # additional assertions in survey::as.svrepdesign
  checkmate::assert_numeric(scale, lower = 1, len = 1)
  checkmate::assert_numeric(awfactt_scale, len = 1)
  checkmate::assert_logical(counts, len = 1)
  checkmate::assert_logical(clustcounts, len = 1)

  # check for expected columns in `data`
  check_vars <- c(all.vars(clusters), all.vars(strata), weight, ceb)
  if (!is.null(by)) check_vars <- c(check_vars, all.vars(by))
  if (!is.null(agegr)) check_vars <- c(check_vars, dob, intv)
  checkmate::assert_names(names(data), must.include = check_vars)

  data$children <- data[[ceb]]
  data$weights <- data[[weight]] / mean(data[[weight]])

  if (!is.null(awfactt)) {
    checkmate::assert_names(names(data), must.include = awfactt)
    data$awfactt <- data[[awfactt]] / awfactt_scale
  } else {
    data$awfactt <- 1
  }

  if(is.null(by))
    by <- ~1

  # create age groupings
  if(!is.null(agegr)){
    data$dob <- data[[dob]]
    data$intv <- data[[intv]]

    data$agegr <- cut(
      x = data[[intv]] - data[[dob]],
      breaks = agegr*scale, labels = .epis_labels(agegr),
      include.lowest = TRUE, right = FALSE
    )
    by <- update(by, ~. + agegr)
  }

  # tabulate numerators and denominator
  vars <- unique(unlist(lapply(c(by, strata, clusters), all.vars)))
  aggr <- dplyr::group_by(data, across(vars)) %>%
    dplyr::summarise(
      children = sum(children * weights),
      women = sum(weights * awfactt)
    ) %>%
    data.frame()

  ## construct interaction of all factor levels that appear
  byvar <- intersect(c(all.vars(by), "agegr", "cohort"),
                     names(aggr))
  aggr$byf <- interaction(aggr[byvar], drop=TRUE)

  ## prediction for all factor levels that appear
  pred <- data.frame(aggr[c(byvar, "byf")])[!duplicated(aggr$byf),]
  pred <- pred[order(pred$byf), ]

  if(counts || varmethod == "none"){
    mc <- model.matrix(~-1+byf, aggr)
    clong <- aggr[c("children", "women")]
    pred[c("children", "women")] <- t(mc) %*% as.matrix(clong)
  }

  if (varmethod == "none") {
    pred$mean <- pred$children / pred$women
    pred$byf <- NULL
    if (!counts)
      pred[c("children", "women")] <- NULL
  } else if (varmethod == "lin") {

    des <- survey::svydesign(ids=clusters, strata=strata, data=aggr, weights=~1)

    ## fit model
    f <- if(length(levels(aggr$byf)) == 1)
      children ~ offset(log(women))
    else
      children ~ -1 + byf + offset(log(women))

    mod <- survey::svyglm(f, des, family=quasipoisson)

    ## prediction for all factor levels that appear
    pred$women <- 1

    pred_mean <- predict(mod, pred, type="response", vcov=TRUE)
    v <- vcov(pred_mean)
    dimnames(v) <- list(pred$byf, pred$byf)

    pred$mean <- as.numeric(pred_mean)
    pred$se <- sqrt(diag(v))
    pred[c("byf", "women")] <- NULL
    attr(pred, "var") <- v
  } else if (varmethod %in% c("jkn", "jk1")) {

    ## Convert to array with events and PYs for each cluster
    ## reshape2::acast is MUCH faster than stats::reshape
    children_clust <- reshape2::acast(aggr, update(clusters, byf ~ .), value.var="children")
    women_clust <- reshape2::acast(aggr, update(clusters, byf ~ .), value.var="women")

    if (varmethod == "jkn"){
      aggr$strataid <- as.integer(interaction(aggr[all.vars(strata)], drop=TRUE))
      strataid <- drop(reshape2::acast(unique(aggr[c(all.vars(clusters), "strataid")]),
                                       update(clusters,  1 ~ .), value.var="strataid"))
    } else
      strataid <- NULL

    estdf <- jackknife(children_clust, women_clust, strataid)

    pred$mean <- estdf$est
    pred$se <- estdf$se
    attr(pred, "var") <- vcov(estdf)
    pred$byf <- NULL
    if(clustcounts){
      attr(pred, "children_clust") <- children_clust
      attr(pred, "women_clust") <- women_clust
      attr(pred, "strataid") <- strataid
    }
  } else {
    stop(paste0("varmethod = \"", varmethod, "\" is not recognized."))
  }

  rownames(pred) <- NULL

  return(pred)
}
