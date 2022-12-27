#' Dummy Extension Scoring
#' 
#' Scores a matrix or dataframe on latent factors, using Dummy Extension Scoring (Keenan, 2023). The factors to be scored can either be provided in the form of a psych::fa output, or a matrix of loadings. Alternatively, if only a data matrix is provided, the factor analysis will be performed within the dummy_score function, with parameter passed by '...' The user can select 'core.variables' to be factored, and also 'scoring.variables' to be used in the scoring step. The 'scoring.variables' can include factored variables, variables not factored, or any combination thereof.
#'
#' @param data A matrix or dataframe containing the data to be analyzed
#' @param fa (optional) An object of class 'fa' from the 'fa' function in the 'psych' package, used to define the factors to be scored
#' @param loadings (optional) Alternatively, a matrix of factor loadings
#' @param method Which Dummy-Extension method to use. Options are 'Rubin','Bartlett','Harman', and 'Anderson.' Rubin is similar to the Thurstone method in traditional factor analysis, and Bartlett and Harman are precisely the same. 'Anderson' is short for the Anderson/Rubin method that produces uncorrelated factor scores.
#' @param parallel A logical indicating whether to use multiple cores for scoring. This will only offer speed improvement for very large datasets, and will probably be slower for most datasets.
#' @param core.variables (optional) Can be either a vector of variable names or index numbers. Defines the 'core' variables from which the factor extension is done. If 'fa' and 'loadings' and 'corr.matrix' are all null, it will use all the variables in 'data'
#' @param scoring.variables (optional) Either a vector of variable names of index numbers. Defines the variables used in the scoring step, which can be 'core' variables or other variables, or any combination thereof. Defaults to all the variables in 'data'
#' @param corr.matrix (optional) A correlation matrix of the 'core' variables, used in the factor extension step
#' @param ... (optional) Arguments passed to the psych::fa function if no factor loadings are provided. See psych::fa for more details.
#' @return A list with elements 'scores,' 'original.loadings,' and 'extension.loadings.' The factor scores for every observation are contained in 'scores.' The loadings of the core variables are in 'original.loadings' and the category loadings of the scoring variables are in 'extension.loadings'
#' @importFrom MASS ginv
#' @importFrom fastDummies dummy_cols
#' @importFrom future plan multisession
#' @importFrom future.apply future_apply
#' @importFrom psych fa
#' @export
#' @examples
#' \dontrun{big.five.scores <- dummy_score(data = bigfive, method = 'Bartlett')$scores}

dummy_score <- function(data, fa=NULL, loadings=fa$loadings, method=c('Rubin','Bartlett','Harman','Anderson'), parallel=F,
                       core.variables=NULL, scoring.variables=NULL, corr.matrix=NULL, score=T, P.ef=NULL, reorder=F, rows=T, ...){
  
  method <- match.arg(method) # If no method specified, use the Rubin method
  
  ## Helper functions
  
  na.zeros <- function(x) replace(x,is.na(x), value=0) # This function replaces NA values with zeros
  in.data <- function(x) x[x %in% colnames(data)] # Returns the elements that are in the column names of 'data'
  `%||%` <- function(a,b) if (is.null(a)) b else a # Returns a if it is not NULL, and b otherwise
  inv <- MASS::ginv # We will use generalized inverses throughout
  setNames <- stats::setNames # For clarity
  
  # A wrapper to psych::target.rot that returns a consistent output no matter the number of factors
  target.rot <- function(x,keys) if(ncol(keys)==1) psych::target.rot(x,keys) else psych::target.rot(x,keys)$loadings

  # Bins a continuous vector into 5 equal-sized categories
  cut.props <- function(x) quantile(x,1:4/5) |> c(-Inf,Inf) |> sort() |> cut(x=x,labels=F) 
  
  # This function does ML-FA with correlation matrix R and starting values P, used only if reorder=T
  ml.fa <- function(R,f,P=NULL){
    for (i in 1:100){
      U.inv <- diag(1 / (1-rowSums(P^2)))
      G <- U.inv - (U.inv %*% P) %*% inv(diag(1,f) + t(P) %*% U.inv %*% P) %*% (t(P) %*% U.inv)
      P <- R %*% G %*% P %*% inv(t(P) %*% G %*% R %*% G %*% P + diag(1,f) - t(P) %*% G %*% P)
    }
    return(P)
  }
  
  # An inverse matrix square root, which is used in the Anderson method
  invMatSqrt <- function(x, correct=F){
    if (dim(x)[1] == 1) return(sign(x) * 1 / sqrt(abs(x)))
    e <- svd(x)
    changed <- diag(1 / sqrt(e$d))
    if (correct) changed[e$d < 1] <- 0 else changed[e$d < .000001] <- 0
    return(e$u %*% changed %*% t(e$u))
  }
  
  extension <- function(R.ev,R.vv,P.vf){
    U.inv <- diag(1 / (1 - rowSums(P.vf^2)))
    R.ev %*% inv(R.vv) %*% P.vf
  }
  
  # This function takes in a response vector X and scores it according to the method
  score.response <- function(X,P,R,f,U.inv,method=c('Rubin','Bartlett','Harman','Anderson')){
    not.na <- !is.na(X)
    P <- P[not.na,]; R <- R[not.na, not.na]; U.inv <- U.inv[not.na, not.na]; X <- X[not.na]
    W <- switch(method,
            Harman = P %*% inv(t(P) %*% P),
            Bartlett = U.inv %*% P %*% inv(t(P) %*% U.inv %*% P),
            Rubin = U.inv %*% P %*% inv(diag(1,f) + t(P) %*% U.inv %*% P),
            Anderson = U.inv %*% P %*% invMatSqrt(t(P) %*% U.inv %*% R %*% U.inv %*% P))
    #W <- W / sqrt(diag(t(W) %*% R %*% W))
    validity <- (diag(t(P) %*% W) / sqrt(diag(t(W) %*% R %*% W))) |> min()
    return(c(X %*% W, validity))
  }

  ## Now the actual content of the main function
  
  # If no dimension names for data, use integers
  dimnames(data) <- list(rownames(data) %||% 1:nrow(data),colnames(data) %||% 1:ncol(data))
  
  # Create numeric and factor versions of the data matrix; discretize any continuous variables for the factor version
  
  X.numeric <- data |> apply(2,as.numeric)
  X.factor <- X.numeric |> apply(2,\(x) if (F %in% (x %% 1 == 0)) cut.props(x) else x) |> 
    as.data.frame() |> lapply(as.factor) |> as.data.frame() |> setNames(colnames(X.numeric))
  
  # Depending on what has been provided as arguments, choose core and extension variables
  
  v <- core.variables %||% rownames(loadings) %||% colnames(corr.matrix) %||% colnames(data) |> in.data()
  e <- scoring.variables %||% colnames(data) |> in.data()
  
  # Create a matrix of dummy variables and center it
  X.ne <- X.factor[, e] |> fastDummies::dummy_cols(ignore_na=T, remove_selected_columns=T) |> scale(scale=F)
  
  # Calculate core-variable correlation matrix if not given, and replace NAs with zeros
  R.vv <- (if (is.null(corr.matrix)) stats::cor(X.numeric[,v], use='pairwise.complete.obs') else corr.matrix[v, v]) |> na.zeros()

  # If loadings not given, run a factor analysis, with arguments passed down

  P.vf <- if (is.null(loadings)) psych::fa(R.vv, ...) |> (`$`)('loadings') else loadings[v,] # P.vf is the core loadings matrix
  f <- ncol(P.vf) # Number of factors

  # Calculate extension correlation matrix if not given, and replace NAs with zeros
  R.ev <- stats::cor(X.ne, X.numeric[,v], use='pairwise.complete.obs') |> suppressWarnings() |> na.zeros()

  # Using the Gorsuch factor extension method, find loadings and uniquenesses for dummy variables
  
  #P.ef <- if (is.null(P.ef)) R.ev %*% inv(R.vv) %*% P.vf else P.ef
  P.ef <- if (is.null(P.ef)) extension(R.ev,R.vv,P.vf) else P.ef
  R.ee <- tcrossprod(P.ef) |> (`diag<-`)(1)           # R.ee = P.ef x t(P.ef) + U.ee
  U.inv <- diag(1 / (1 - rowSums(P.ef^2)))            # Inverse of uniqueness
  
  # If exploratory category re-ordering is requested (see Keenan (2023) for more info)

  if (reorder){
    for (i in 1:f){
      X.numeric.new <- X.numeric |> apply(2,\(x) if (F %in% (x %% 1 == 0)) cut.props(x) else x) |> 
        setNames(colnames(X.numeric))
      X.numeric.new <- lapply(v,\(x){
        P.ef |> (`[`)(rownames(P.ef) |> sapply(\(.) sub("\\s*\\_\\d+", "", x=.))==x,i,drop=F) |> 
          order() |> factor(X.numeric.new[,x],levels=_) |> as.integer()
      }) |> as.data.frame() |> setNames(v) 
      P.target <- cor(X.numeric.new[rows,],X.numeric[rows,v],use='pairwise.complete.obs') |> 
                     na.zeros() |> extension(R.vv,P.vf)
      R.vv.new <- X.numeric.new[rows,] |> cor(use='pairwise.complete.obs') |> na.zeros()
      P.vf.new <- R.vv.new |> ml.fa(f,P.target) |> target.rot(P.target)
      P.ef[,i] <- cor(X.ne,X.numeric.new,use='pairwise.complete.obs') |> 
                     na.zeros() |> extension(R.vv.new,P.vf.new) |> subset(select=i)
    }
  }
  
  # Now do the row-by-row scoring of the dummy variable matrix, in parallel if requested
  
  if (!score) return(list(scores=NULL, original.loadings=loadings, extension.loadings=P.ef))
  
  if (parallel){
    future::plan(future::multisession)
    scores <- X.ne |> future.apply::future_apply(1, score.response, P.ef, R.ee, f, U.inv, method) |> 
      matrix(ncol=f+1, byrow=T, dimnames=list(rownames(data),c(colnames(P.vf),'Validity')))
    future::plan(future::sequential)
    } 
  else scores <- X.ne |> apply(1, score.response, P.ef, R.ee, f, U.inv, method) |> 
    matrix(ncol=f+1, byrow=T, dimnames=list(rownames(data),c(colnames(P.vf),'Validity')))

  # Return the estimated scores, original loadings, and extension loadings
  
  return(list(scores=scores, original.loadings=P.vf, extension.loadings=P.ef))
}
