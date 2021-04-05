#' @title Calculate mean children ever born (ceb) & living children (lc)
#'
#' @description Use default arguments to calculate mean children ever born
#' related statistics from a DHS individual recode file.
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
#' @param ceb \[`character()`\]\cr
#'   Variable names for children ever born related columns to calculate the mean
#'   & standard error of. Value corresponds to name in `data`, name corresponds
#'   to new variable name in output. Default is 'c(ceb = "v201", lc = "v218")'.
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
#' @return \[`data.frame()`\] consisting of estimates and standard errors as
#' returned by `survey::svyby`. The full covariance matrix of the estimates can
#' be retreived by `vcov(val)`.
#'
#' @examples
#' data(zzir)
#' calc_ceb(zzir)
#'
#' @export
calc_ceb <- function(data,
                     by = NULL,
                     agegr = 3:10*5,
                     clusters = ~v021,
                     strata = ~v024+v025,
                     weight = "v005",
                     ceb = c(ceb = "v201", lc = "v218"),
                     dob = "v011",
                     intv = "v008",
                     varmethod = "lin",
                     scale = 12) {

  checkmate::assert_data_frame(data)
  checkmate::assert_formula(by, null.ok = TRUE)
  checkmate::assert_formula(clusters)
  checkmate::assert_formula(strata)
  checkmate::assert_string(weight)
  checkmate::assert_character(ceb, min.len = 1)
  checkmate::assert_named(ceb)
  checkmate::assert_string(varmethod) # additional assertions in survey::as.svrepdesign
  checkmate::assert_numeric(scale, lower = 1, len = 1)

  # check for expected columns in `data`
  check_vars <- c(all.vars(clusters), all.vars(strata), weight, ceb)
  if (!is.null(by)) check_vars <- c(check_vars, all.vars(by))
  if (!is.null(agegr)) check_vars <- c(check_vars, dob, intv)
  checkmate::assert_names(names(data), must.include = check_vars)

  # standardize variable names
  for (new_name in names(ceb)) {
    old_name <- ceb[new_name]
    data[[new_name]] <- data[[old_name]]
  }

  data$weights <- data[[weight]] / mean(data[[weight]])

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

  des <- survey::svydesign(
    ids = clusters,
    strata = strata,
    weights = ~weights,
    data = data,
  )
  if (varmethod != "lin") des <- survey::as.svrepdesign(des, type = varmethod)

  pred <- survey::svyby(
    formula = formula(paste("~", paste(names(ceb), collapse = "+"))),
    by = by,
    design = des,
    FUN = survey::svymean,
    covmat = TRUE
  )

  rownames(pred) <- NULL
  if (varmethod != "lin") {
    # survey::svyby returns se columns as 'se1' & 'se2' with non 'lin' varmethod
    names(pred)[grepl("^se", names(pred))] <- paste0("se.", names(ceb))
  }
  return(pred)
}
