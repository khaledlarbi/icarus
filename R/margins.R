# copyright (C) 2014-2016 A.Rebecq

# TODO : add xtable parameters to marginsToTeX parameters
marginsToTeX <- function(marginMatrix, names=NULL, pct=FALSE, popTotal=NULL,
                         scaleboxTeX=NULL, file=NULL,
                         label=NULL, caption=NULL) {
  
  if (!requireNamespace("xtable", quietly = TRUE)) {
    stop("Package xtable needed for export of margins in LateX to work. Please install it.",
         call. = FALSE)
  }
  
  if(!is.matrix(marginMatrix)) {
    stop("marginsToTeX input type has to be matrix.")
  }
  
  if(!is.null(names)) {
    if(length(names) != nrow(marginMatrix)) {
      stop("Name length must equal number of rows in marginMatrix")
    }
    
    marginMatrix[,1] <- names
  }
  
  if(pct) {
    numericPart <- marginMatrix[,3:(ncol(marginMatrix))]
    numericPart <- as.numeric(numericPart)
    numericPart <- 100*numericPart
    numericPart -> marginMatrix[,3:(ncol(marginMatrix))]
  }
  
  # Write zeros as NA
  marginMatrix[as.numeric(marginMatrix) == 0] <- NA

  marginDF <- as.data.frame(marginMatrix)
  
  # Heuristic rule for scalebox
  if(is.null(scaleboxTeX)) {
    if(ncol(marginDF) >= 10) {
      scaleboxTeX <- 1.4 - ncol(marginDF) / 20
    }
    
    if(ncol(marginDF) >= 28) {
      stop("Automatic scaleboxing not configured for more than 28 margins.")
    } 
  }
  
  captionTeX <- caption
  if(!is.null(popTotal)) {
    captionTeX <- paste(caption, " -- total population : ", round(popTotal,0),sep="")
  }
  
  print(xtable::xtable(marginDF, caption=captionTeX, label=label), include.rownames = FALSE, include.colnames = FALSE,
               floating = TRUE, scalebox=scaleboxTeX, file=file
        )

}

#' Create empty margin matrix
#' @description 
#' Use this to create an empty margin matrix (which facilitates
#' the use of magrittr syntax to enter margins)
#' 
#' @examples 
#' library(magrittr)
#' N <- 230 ## population total
#' ## Horvitz Thompson estimator of the mean: 2.174
#' weightedMean(data_employees$movies, data_employees$weight, N)
#' ## Enter calibration margins:
#' margins <- newMarginMatrix() %>%
#'   addMargin("category", c(0.35, 0.40, 0.25)) %>%
#'   addMargin("sex", c(0.6, 0.4)) %>%
#'   addMargin("department", c(0.45, 0.55)) %>%
#'   addMargin("salary", 470000)
#' ## Compute calibrated weights with raking ratio method
#' wCal <- calibration(data=data_employees, marginMatrix=margins, colWeights="weight"
#'                     , method="raking", pct = TRUE, description=FALSE
#'                     , popTotal = N)
#' ## Calibrated estimate: 2.471917
#' weightedMean(data_employees$movies, wCal, N)
#' @export
newMarginMatrix <- function() {
  return(matrix(, nrow = 0, ncol = 1))
}


#' Consistency Check for Margins and Linear Combinations
#'
#' The function checks the consistency between a matrix X and a vector of margins "margins".
#' It returns TRUE if they are consistent, otherwise, it returns FALSE. The consistency is
#' checked by detecting linear dependencies in X and verifying whether the margins follow
#' the expected relationships for those dependencies.
#'
#' @param X A numeric matrix (n x p) where n is the number of observations and p is the 
#'          number of variables. It represents the data on which consistency with margins
#'          is to be checked.
#' @param margins A numeric matrix of size (p x 1). It represents the margins 
#'                (e.g., observed sums or marginal totals) for the variables in X.
#' @param tol A numeric value specifying the tolerance for checking consistency. Defaults to 
#'            1e-7. If the sum of squared differences between linear combinations and the 
#'            margins is less than this value, the consistency is considered true.
#'
#' @return A logical value: TRUE if the variables in X are consistent with the margins, 
#'         otherwise FALSE. If any inconsistencies are found, additional details about the 
#'         problematic variables and linear combinations are printed.
#'
#' @details
#' The function performs the following steps:
#' 1. Checks if X is a matrix with two dimensions (n x p) and if margins is a matrix with 
#'    one dimension (n x 1).
#' 2. Detects if there are any linear combinations in X using the `caret::findLinearCombos` 
#'    function.
#' 3. If linear dependencies are found, it calculates the weights (or coefficients) needed 
#'    to approximate the margins based on these dependencies.
#' 4. For each identified linear combination, it calculates the sum of squared differences 
#'    between the predicted and actual margins, and checks if these differences are within 
#'    the specified tolerance.
#' 5. If any inconsistencies are found, it prints detailed information about the problematic 
#'    linear combinations and the associated margins.
#'
#' @examples
#' ind <- c(rep(1,5), rep(0,5))
#' X <- cbind(ind, 1-ind, 1)
#' margins <- matrix(c(10,20,30), ncol = 1)
#' margins2 <- matrix(c(10,20,31), ncol = 1)
#' consistency_aux_var(X, margins)
#' consistency_aux_var(X, margins2)
#' consistency_aux_var <- function(X, margins, tol = 1e-7)
consistency_aux_var <- function(X, margins, tol = 1e-7){
  # The function "consistency_aux_var" checks the consistency between a matrix X 
  # and a vector of margins "margins". It returns TRUE if they are consistent,
  # otherwise, it returns FALSE.
  
  # Initialize "is_consistent" as TRUE, indicating that, by default,
  # the variables are considered consistent.
  is_consistent <- TRUE
  
  # Initialize "not_issues" as TRUE, to check if any issues are detected.
  not_issues <- TRUE
  
  # Check if "X" is a matrix with two dimensions (n x p).
  if((!is.matrix(X)) | (length(dim(X)) != 2)){
    stop("X must be a matrix with two dimensions.")
  }
  
  # Check if "margins" is a matrix with one dimension (n x 1).
  if((!is.matrix(margins)) | (length(dim(X)) != 2)){
    stop("margins must be a matrix with one dimension.")
  }
  
  if(ncol(X) != nrow(margins)){
    stop("Number of columns from X must be equals to the number of columns from margins")
  }
  
  # Check the rank of matrix X. The idea here is to detect if X has
  # any linearly dependent variables. This step is crucial to ensure
  # that the variables are linearly independent in X.
  lin_combo <- caret::findLinearCombos(X)
  
  # If linear combinations are detected, meaning if the length of
  # "lin_combo$remove" is greater than 0, perform the following steps.
  if(length(lin_combo$remove) > 0){
    # "tX" contains the transpose of X, after removing columns that are
    # linearly dependent.
    tX <- t(X)[-lin_combo$remove, ,drop = FALSE]
    
    # "margins_rest" contains the margins after removing the rows corresponding
    # to the linearly dependent columns of X.
    margins_rest <- margins[-lin_combo$remove, ,drop = FALSE]
    
    # "pond" is a weight vector that solves the linear system 
    # to approximate the remaining margins.
    pond <- MASS::ginv(tX) %*% margins_rest
    
    # Calculate the differences between the margins and the predicted values
    # based on the linear combination, for each identified linear combination.
    diff <- lapply(X = lin_combo$linearCombos,
                   FUN = function(indice){
                     # Calculate the difference between the linear combination and the expected margin.
                     diff <- (t(X)[indice, , drop = FALSE] %*% pond) - margins[indice, , drop = FALSE]
                     return(sum(diff^2))  # Return the sum of squared differences.
                   })
    
    # If the sum of squared differences is below the tolerance (tol), 
    # it means the linear combination is consistent with the margins.
    not_issues <- abs(unlist(diff)) < tol
    
    # If all consistency checks pass, update "is_consistent".
    is_consistent <- is_consistent & all(not_issues)
  }
  
  # If any inconsistencies are found, identify and display the problematic variables.
  if(any(!not_issues)){
    # Get the indices of the variables with inconsistencies.
    ind_issues <- which(!not_issues)
    
    # Display a message indicating that there are consistency issues.
    cat("Some variables are linear combinations of others, but the margins do not adhere to this linear relationship. Please check the data.\n")
    
    # Loop over each index of problematic variables and print the corresponding linear combinations
    lapply(ind_issues, function(i){
      # Access the linear combinations for the current index i
      lin_combo_var_issues <- lin_combo$linearCombos[[i]]
      # Print the problematic variables in the linear combination
      cat("Linear combinaison : ", lin_combo_var_issues[1], " ~ ", paste(lin_combo_var_issues[-1], collapse = " - "), "\n")
      coef <- MASS::ginv(X[, lin_combo_var_issues[-1], drop = FALSE]) %*% X[,lin_combo_var_issues[1],drop = FALSE]
      cat("Margin based on ", lin_combo_var_issues[1], " : ", margins[lin_combo_var_issues[1]], "\n")
      cat("Margin based on ", paste(lin_combo_var_issues[-1], collapse = " - "), " : ",
          t(coef) %*% margins[lin_combo_var_issues[-1], , drop = FALSE], "\n")
      cat("---------------------- \n")
    })
  }
  
  
  # Return TRUE if the variables are consistent, otherwise FALSE.
  return(is_consistent)
}