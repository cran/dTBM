#' Peru Legislation data
#' 
#' The Peru Legislation data is obtained by preprocessing the original data in Lee et al., 2017. 
#'
#' The data frame "attr_data" is a 116 x 2 matrix consisting the name and party affiliation of 116 legislators
#' in the top five parties. The legislators IDs are recorded in the row names of the matrix.
#' 
#' The data frame "laws_data" is a 5844 x 7 matrix recording the co-sponsorship of 116 legislators of 802 bills during the first half of 2006-2007 year.
#' 
#' The array "network_data" is a 116 x 116 x 116 binary tensor recording the presence of order-3 co-sponsorship among legislators based on "laws_data". 
#' Specfically, the entry (i,j,k) is 1 if the legislators (i,j,k) have sponsored the same bill, and the entry (i,j,k) is 0  otherwise.
#'
#'@docType data
#'@usage data(peru)
#'
#'@format  A list. Includes a 116-2 data frame named "attr_data", a 5844-7 data frame named "laws_data", and a 116-116-116 binary array named "network_data".
#'
#'@keywords datasets
#'
#'
"peru"