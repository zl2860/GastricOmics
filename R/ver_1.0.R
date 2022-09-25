#' Quickly boxing continuous matrix
#'
#' @param m an input matrix
#' @param levels number of categories to be cut, depending on the quantiles (using ::cut())
#'
#' @return the output matrix
#' @export
#'
#' @examples data_boxing(matrix, 2)
matrix_boxing <- function(m, levels) {
  res <- apply(m, 2, function(x) as.numeric(cut(x, levels)))
  if (!is.null(colnames(res))) {
    colnames(res) <- paste(colnames(res), "_cat", sep = "")
  }
  return(res)
}


#' Title Quickly obtain pairwise interaction matrix
#'
#' @param df an input data.frame object
#'
#' @return interaction matrix
#' @export
#'
#' @examples pairwise_interation(df)
matrix_pairwise_interation <- function(df) {
  # df: an input data.frame object
  res <- model.matrix(~ .^2, df)
  return(res)
}


#' Title Get the list of interators for proteins base on the BioGrid database
#'
#' @param list an input list containing proteins to be analyzed
#'
#' @return an output list containing the original proteins with their interactors
#' @export
#'
#' @examples match_protein_interaction(pro_list)
match_protein_interaction <- function(list) {
  library(simpIntLists)
  # get full list of proteins
  data("HumanBioGRIDInteractionOfficial")
  protein_list_complete <- NULL
  for (i in 1:length(HumanBioGRIDInteractionOfficial)) {
    protein_list_complete[i] <- HumanBioGRIDInteractionOfficial[[i]]$name
  }

  # input: a list of proteins
  # output: a list of the proteins with their interactors
  idx <- which(protein_list_complete %in% list)
  target_list <- HumanBioGRIDInteractionOfficial[idx]
  return(target_list)
}

# 快速计算fold change


# 快速格式化ggplot
