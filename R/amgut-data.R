#' Microbiome family-level data from American Gut dataset
#'
#' Counts from 16S sequencing on the American Gut dataset.  Clinical covariates
#' also include
#'
#' @docType data
#'
#' @usage data(amgut)
#'
#' @format List with two objects: \code{clin.dat} is a dataframe with 100 rows and 4 columns containing the covariate data; \code{counts} is a matrix with 100 rows and 10 columns which contains the family-level counts of each subject.  The last column of \code{otu.counts} is a reference category containing the sum of all families not included in this dataset.
#'
#' @keywords datasets american gut microbiome
#'
#' @references McDonald, Daniel, et al. "American gut: an open platform for citizen science microbiome research." MSystems 3.3 (2018): e00031-18.
#'
#' @source \href{https://qiita.ucsd.edu/study/description/10317}{Qiita database}
#'
#' @examples
#' data(amgut)
#'
"amgut"
