#' Dummy Extension Scoring
#' 
#' Scores a matrix or dataframe on latent factors, using Dummy Extension Scoring (Keenan, 2023). The factors to be scored can either be provided in the form of a psych::fa output, or a matrix of loadings. Alternatively, if only a data matrix is provided, a factor analysis will be performed using parameters passed by '...' The user can select 'core.variables' to be factored, and also 'scoring.variables' to be used in the scoring step. The 'scoring.variables' can include factored variables, variables not factored, or any combination thereof.
#'
#' @param data A matrix or dataframe containing the data to be analyzed
#' @param fa (optional) An object of class 'fa' from the 'fa' function in the 'psych' package, used to define the factors to be scored
#' @param loadings (optional) Alternatively, a matrix of factor loadings
#' @param method Which Dummy-Extension method to use. Options are 'Rubin','Bartlett','Harman', and 'Anderson.' Rubin is similar to the Thurstone method in traditional factor analysis, and Bartlett and Harman are precisely the same. 'Anderson' is short for the Anderson/Rubin method that produces uncorrelated factor scores.
#' @param parallel A T/F indicating whether to use multiple cores for scoring. This will only offer speed improvement for very large datasets, and will probably be slower for most datasets.
#' @param ... Additional arguments, including
#' 
#' \itemize{
#'   \item{Arguments to 'fa' function}{
#'     If no factor loadings are provided, they will be estimated from the data using the 'fa' function in the 'psych' package. 
#'     Arguments to the 'fa' function, such as 'nfactors' or 'rotate,' can be passed directly into 'dummy_score.' See psych::fa for more.
#'     }
#'   \item{core.variables}{
#'     Can be either a vector of variable names or index numbers. 
#'     Defines the 'core' variables from which the factor extension is done. 
#'     If 'fa' and 'loadings' and 'corr.matrix' are all null, it will default to all the variables in 'data'
#'     }
#'   \item{scoring.variables}{
#'    Either a vector of variable names or index numbers. 
#'    Defines the variables used in the scoring step, which can be 'core' variables or other variables, or any combination thereof. 
#'    Defaults to all the variables in 'data'
#'   }
#'   \item{core.rows}{
#'    A vector indicating which rows to include in the 'core' factor analysis. 
#'    If the data are from multiple years, for example, it may be desirable to perform the core analysis on one year, and extend to the rest.
#'  }
#'   \item{score}{
#'  A T/F indicating whether to estimate scores for every individual observation. If you just desire the category loadings, set score to F.
#'  }
#'  \item{corr.matrix}{
#'   A correlation matrix of the 'core' variables, used in the factor extension step. Saves time if this is already known.
#'   }
#'  \item{P.ef}{
#'  Category loadings, if they are already known (perhaps from a previous run of dummy_score).
#'  }
#'  \item{exploratory.reorder}{
#'  A T/F indicating whether to use exploratory category re-ordering. 
#'  First, category loadings are obtained through dummy-extension.
#'  Next, the proper order of the categories for each variable can be ascertained.
#'  Then, we can re-order the categories of variables to make them monotonically related to the factors.
#'  Finally, We can perform another factor analysis with these re-ordered variables, obtaining stronger loadings.
#'  This procedure is experimental and adds to computation time.
#'  }
#' }
#' 
#' @return A list with elements 'scores,' 'original.loadings,' and 'extension.loadings.' The factor scores for every observation are contained in 'scores.' The loadings of the core variables are in 'original.loadings' and the category loadings of the scoring variables are in 'extension.loadings'
#' @importFrom MASS ginv
#' @importFrom fastDummies dummy_cols
#' @importFrom future plan multisession
#' @importFrom future.apply future_apply
#' @importFrom psych fa
#' @importFrom Matrix Matrix crossprod
#' @importFrom stats quantile
#' @export
#' @examples
#' \dontrun{big.five.scores <- dummy_score(data = bigfive, method = 'Bartlett')$scores}

dummy_score <- function(data, fa=NULL, loadings=fa$loadings, method=c('Rubin','Bartlett','Harman','Anderson'), parallel=F, ...){
  
  method <- match.arg(method) # If no method specified, use the Rubin method
  extra.args <- list(...)
  
  ## Helper functions
  
  na.zeros <- function(x) replace(x,is.na(x), values=0) # This function replaces NA values with zeros
  in.data <- function(x) x[x %in% colnames(data)] # Returns the elements that are in the column names of 'data'
  `%||%` <- function(a, b) if (is.null(a)) b else a # Returns 'a' if it is not NULL, and 'b' otherwise
  set.colnames <- function(x, values) `colnames<-`(x, values) # For clarity
  inv <- MASS::ginv # We will use generalized inverses throughout
  
  # A wrapper to psych::target.rot that returns a consistent output no matter the number of factors
  target.rot <- function(x,keys) if(ncol(keys)==1) psych::target.rot(x,keys) else psych::target.rot(x,keys)$loadings
  
  # A wrapper to psych::fa that passes along extra arguments smoothly
  psych.fa <- function(r) do.call(psych::fa,append(list(r=r),extra.args)) |> (`$`)('loadings') |> suppressWarnings()

  # Bins a continuous vector into 5 equal-sized categories
  cut.props <- function(x) quantile(x,1:4/5) |> c(-Inf,Inf) |> sort() |> cut(x=x,labels=F) 
  
  # This function does ML-FA with correlation matrix R and starting values P; only used if reorder=T
  ml.fa <- function(R,f,P=NULL){
    d <- 1
    while (d > .0001){
      U.inv <- diag(1 / (1 - rowSums(P^2)))
      G <- U.inv - U.inv %*% P %*% inv(diag(1,f) + t(P) %*% U.inv %*% P) %*% (t(P) %*% U.inv)
      P.star <- R %*% (G %*% P) %*% inv(t(P) %*% G %*% R %*% (G %*% P) + diag(1,f) - t(P) %*% G %*% P)
      d <- (P - P.star) |> abs() |> max()
      P <- P.star
    }
    return(P)
  }
  
  # A function to do correlation, utilizing sparse matrices
  sparse.cor <- function(x,y=NULL){
    if(is.null(y)) y <- x
    x <- scale(x)
    y <- scale(y)
    x.cases <- ifelse(is.na(x),0,1) |> Matrix::Matrix(sparse=T)
    y.cases <- ifelse(is.na(y),0,1) |> Matrix::Matrix(sparse=T)
    pairwise.n <- crossprod(x.cases,y.cases) + 1
    x <- na.zeros(x) |> Matrix::Matrix(sparse=T)
    y <- na.zeros(y) |> Matrix::Matrix(sparse=T)
    R <- (crossprod(x,y) / pairwise.n) |> as.matrix()
  }
  
  # An inverse matrix square root, which is used in the Anderson method
  invMatSqrt <- function(x){
    if (dim(x)[1] == 1) return(sign(x) * 1 / sqrt(abs(x)))
    S <- svd(x)
    S$d[S$d < 0.01] <- 0.01
    return(S$u %*% diag(1 / sqrt(S$d)) %*% t(S$u))
  }
  
  # Factor Extension using the Gorsuch (1997) method
  extension <- function(X.nv, X.ne, P.vf, R.vv=NULL){
    R.vv <- R.vv %||% sparse.cor(X.nv)
    R.ev <- sparse.cor(X.ne,X.nv)
    S <- svd(R.vv) # Faster than using ginv(R.vv)
    S$d[S$d < .01] <- 0.01
    P.ef <- R.ev %*% (S$u %*% (diag(1/S$d) %*% crossprod(S$v,P.vf)))
  }
  
  # This function takes in a response vector X and scores it according to the method
  score.response <- function(X,P,R,f,U.inv,method=c('Rubin','Bartlett','Harman','Anderson')){
    not.na <- !is.na(X)
    if(sum(not.na)==0) return(c(rep(NA,f),NA))
    P <- P[not.na,]; R <- R[not.na, not.na]; U.inv <- U.inv[not.na, not.na]; X <- X[not.na]
    W <- switch(method,
            Harman = P %*% inv(t(P) %*% P),
            Bartlett = U.inv %*% P %*% inv(t(P) %*% U.inv %*% P),
            Rubin =  U.inv %*% P %*% inv(diag(1,f) + t(P) %*% U.inv %*% P),
            Anderson = U.inv %*% P %*% invMatSqrt(t(P) %*% U.inv %*% R %*% U.inv %*% P))
    validity <- (diag(t(W) %*% P) / sqrt(diag(t(W) %*% R %*% W))) |> min()
    return(c(X %*% W, validity))
  }

  ## Handling additional arguments
  
  # Set defaults for arguments passed in with '...'
  core.variables <- scoring.variables <- corr.matrix <- score <- P.ef <- exploratory.reorder <- core.rows <- NULL
  arg.defaults <- list(core.variables=NULL, scoring.variables=NULL, corr.matrix=NULL, 
                       score=T, P.ef=NULL, exploratory.reorder=F, core.rows=T)
  
  # Assign these extra arguments to variables in the main function environment
  lapply(names(arg.defaults), \(x) assign(x,extra.args[[x]] %||% arg.defaults[[x]], envir=parent.frame(2)))
  
  # Remove them from extra.args, so it can be passed to psych.fa cleanly
  extra.args <- extra.args[!(names(extra.args) %in% names(arg.defaults))]
  
  ## Now the actual content of the main function
  
  # If there are no column names for data, use integers
  colnames(data) <- colnames(data) %||% 1:ncol(data)
  
  # Depending on what has been provided as arguments, choose core and extension variables
  
  v <- core.variables %||% rownames(loadings) %||% colnames(corr.matrix) %||% colnames(data) |> in.data()
  e <- scoring.variables %||% colnames(data) |> in.data()
  
  # Create numeric and factor versions of the data matrix; discretize any continuous variables for the factor version
  
  X.numeric <- data |> apply(2,as.numeric)
  X.factor <- X.numeric |> apply(2,\(x) if (F %in% (x %% 1 == 0)) cut.props(x) else x) |> 
    as.data.frame() |> lapply(as.factor) |> as.data.frame() |> set.colnames(colnames(X.numeric))
  
  # Create a matrix of dummy variables and center it
  X.ne <- X.factor[, e] |> fastDummies::dummy_cols(ignore_na=T, remove_selected_columns=T) |> scale(scale=F)
  
  # Calculate core-variable correlation matrix if not given, and replace NAs with zeros
  R.vv <- corr.matrix[v, v] %||% sparse.cor(X.numeric[core.rows, v]) |> na.zeros()
  
  # If loadings not given, run a factor analysis, with arguments passed down

  P.vf <- loadings[v,] %||% psych.fa(R.vv) # P.vf is the core loadings matrix
  f <- ncol(P.vf) # Number of factors
  
  # Using the Gorsuch factor extension method, find loadings and uniquenesses for dummy variables
  
  P.ef <- P.ef %||% extension(X.numeric[ , v], X.ne, P.vf, R.vv)
  
  # If exploratory category re-ordering is requested-- see Keenan (2023)
  if (exploratory.reorder){
    for (i in 1:f){
      X.numeric.new <- X.numeric |> apply(2,\(x) if (F %in% (x %% 1 == 0)) cut.props(x) else x) |> 
        set.colnames(colnames(X.numeric))
      X.numeric.new <- lapply(v,\(x){
        P.ef |> (`[`)(rownames(P.ef) |> sapply(\(.) sub("\\s*\\_\\d+", "", x=.))==x,i,drop=F) |> 
          order() |> factor(X.numeric.new[,x],levels=_) |> as.integer()
      }) |> as.data.frame() |> as.matrix() |> set.colnames(v)
      P.target <- extension(X.numeric[core.rows,v],X.numeric.new[core.rows,],P.vf,R.vv)
      R.vv.new <- X.numeric.new[core.rows,] |> sparse.cor()
      P.vf.new <- R.vv.new |> ml.fa(f,P.target) |> target.rot(P.target)
      P.ef[,i] <- extension(X.numeric.new,X.ne,P.vf.new,R.vv.new) |> subset(select=i)
    }
  }
  
  R.ee <- tcrossprod(P.ef) |> (`diag<-`)(1)   # R.ee = P.ef x t(P.ef) + U.ee
  U.inv <- diag(1 / (1 - rowSums(P.ef^2)))    # Inverse of uniqueness
  
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
