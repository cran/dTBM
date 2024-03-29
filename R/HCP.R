#' HCP data
#'
#' The HCP data is obtained by preprocessing the data from Human Connectome Project (HCP); see https://wiki.humanconnectome.org/display/PublicData/.
#' 
#' The array "tensor" is a 68 × 68 × 136 binary tensor consisting of structural connectivity patterns among 68 brain regions for 136 individuals. All the individual images were preprocessed following a standard pipeline (Zhang et al., 2018), and the brain was parcellated to 68 regions-of-interest following the Desikan atlas (Desikan et al., 2006). The tensor entries encode the presence or absence of fiber connections between those 68 brain regions for each of the 136 individuals.
#' 
#' The data frame "attr" is a 136 × 573 matrix consisting of 573 personal features for 136 individuals. The full list of covariates can be found at:
#' https://wiki.humanconnectome.org/display/PublicData/
#'
#'@docType data
#'@usage data(HCP)
#'
#'@format  A list. Includes a 68-68-136 binary array named "tensor" and a 136-573 data frame named "attr".
#'
#'@keywords datasets
#'
#'
"HCP"