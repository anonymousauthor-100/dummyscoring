#' Dummy Extension Scoring
#' 
#' Scores a matrix or dataframe on latent factors, using Dummy Extension Scoring (Keenan, 2023). The factors to be scored can either be provided in the form of a psych::fa object, or a matrix of factor loadings. Alternatively, if only a data matrix is provided, a factor analysis will be performed using parameters passed by '...' The user can select 'core.variables' to be factored, and also 'scoring.variables' to be used in the scoring step. The 'scoring.variables' can include factored variables, variables not factored, or any combination thereof.
#'
#' @param data A matrix or dataframe containing the data to be analyzed
#' @param fa (optional) An object of class 'fa' from the \code{fa()} function in the 'psych' package, used to define the factors to be scored
#' @param loadings (optional) Alternatively, a matrix of factor loadings
#' @param method Which Dummy-Extension method to use. Options are 'Ledermann','Bartlett','Harman', and 'Anderson.' Ledermann is similar to the Thurstone method in traditional factor analysis, and Bartlett and Harman are precisely the same. 'Anderson' is short for the Anderson/Rubin method that produces uncorrelated factor scores.
#' @param parallel A T/F indicating whether to do scoring in parallel. This is FALSE by default, but can significantly speed up runtime.
#' @param ... Additional arguments, including
#' 
#' \describe{
#'   \item{Arguments to \code{fa()} function}{
#'     If no factor loadings are provided, they will be estimated from the data using the \code{fa()} function in the 'psych' package. 
#'     Arguments to the \code{fa()} function, such as 'nfactors' or 'rotate,' can be passed directly into \code{dummy_score()}. See \code{psych::fa()} for more.
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
#' @importFrom future availableWorkers
#' @importFrom psych fa
#' @importFrom Rcpp evalCpp
#' @importFrom stats quantile
#' @useDynLib dummyscoring, .registration = TRUE
#' @export
#' @examples
#' \dontrun{big.five.scores <- dummy_score(data = bigfive, method = 'Bartlett')$scores}

dummy_score <- function(data, fa=NULL, loadings=fa$loadings, 
                        method=c('Ledermann','Bartlett','Harman','Anderson','Keenan','Thurstone','Heermann'), 
                        parallel=F, ...){

  ## Helper functions
  
  na.zeros <- function(x) replace(x, is.na(x), values=0) # Replaces NA values with zeros
  in.data <- function(x) x[x %in% colnames(data)] # Returns those elements in the column names of 'data'
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a # Returns 'a' if it exists, and 'b' otherwise
  catt <- function(x) if(verbose) cat(x) # Conditional print statement
  inv <- MASS::ginv # We will use generalized inverses throughout
  'print.dummyscore' |> assign(function(x, ...) cat('Done!\n'), parent.frame()) # Controls return message
  
  # A wrapper to psych::fa that passes along extra arguments smoothly
  psych.fa <- function(r) do.call(psych::fa, append(list(r=r), extra.args)) |> (`$`)('loadings') |> suppressWarnings()

  # Bins a continuous vector into 5 equal-sized categories
  cut.props <- function(x) quantile(x, 1:4/5) |> c(-Inf,Inf) |> sort() |> findInterval(x=x, all.inside=T) 
  
  # This function does ML-FA with correlation matrix R and starting values P; only used if reorder=T
  # Uses the EM algorithm from Rubin and Thayer (1982)
  ml.fa <- function(R, f, P){
    for (i in 1:20){
      U.inv <- diag(1 / (1 - rowSums(P^2)))
      G <- U.inv - U.inv %*% P %*% inv(diag(1, f) + t(P) %*% U.inv %*% P) %*% t(P) %*% U.inv
      P <- R %*% (G %*% P) %*% inv(t(P) %*% G %*% R %*% (G %*% P) + diag(1, f) - t(P) %*% G %*% P)
    }
    return(P)
  }
  
  # Factor Extension using the Gorsuch (1997) method
  extension <- function(X.nv, X.ne, P.vf, R.vv=NULL){
    R.vv <- R.vv %||% sparse.cor(X.nv)
    R.ev <- sparse.cor(X.ne, X.nv)
    S <- svd(R.vv)
    d <- (1 / S$d) |> (`[<-`)(S$d < 0.01, value=0) |> diag()
    P.ef <- R.ev %*% (S$u %*% (d %*% crossprod(S$v, P.vf)))   # P.ef = R.ev %*% inv(R.vv) %*% P.vf
  }
  
  category.reorder <- function(X.nv, X.ne, P.vf, P.ef, R.vv, R.ee, f, core.rows, v){
    for (i in 1:f){
      X.new <- X.nv |> apply(2, \(x) if (F %in% (x %% 1 == 0)) cut.props(x) else x) |> 
        (`dimnames<-`)(dimnames(X.nv))
      X.new <- lapply(v, \(x){
        P.ef |> (`[`)(rownames(P.ef) |> vapply(\(.) sub("\\s*\\_\\d+", "", x=.), character(1)) == x, i, drop=F) |> 
          order() |> factor(X.new[ , x],levels=_) |> as.integer()
      }) |> as.data.frame() |> as.matrix() |> (`colnames<-`)(v)
      R.new <- X.new[core.rows, ] |> sparse.cor()
      P.new <- extension(X.nv[core.rows, v], X.new[core.rows, ], P.vf, R.vv) |> ml.fa(R.new, f, P=_)
      P.ef[ , i] <- extension(X.new, X.ne, P.new, R.new) |> subset(select=i)
    }
    P.ef
  }
  
  ## Handling arguments
  
  method <- match.arg(method) # If no method specified, use the Ledermann method
  extra.args <- list(...)
  
  # Set defaults for arguments passed in with '...'
  core.variables <- scoring.variables <- corr.matrix <- score <- P.ef <- exploratory.reorder <- core.rows <- verbose <- NULL
  arg.defaults <- list(core.variables=NULL, scoring.variables=NULL, corr.matrix=NULL, 
                       score=T, P.ef=NULL, exploratory.reorder=F, core.rows=T, verbose=T)
  
  # Assign these extra arguments to variables in the main function environment
  for (i in names(arg.defaults)) assign(i, extra.args[[i]] %||% arg.defaults[[i]])
  extra.args <- extra.args[!(names(extra.args) %in% names(arg.defaults))]
  
  # If there are no column names for data, use integers
  colnames(data) <- colnames(data) %||% 1:ncol(data)
  
  # Depending on what has been provided as arguments, choose core and extension variables
  v <- core.variables %||% rownames(loadings) %||% colnames(corr.matrix) %||% colnames(data) |> in.data()
  e <- scoring.variables %||% colnames(data) |> in.data()
  
  ## Now the actual content of the main function
  
  # Create numeric and category versions of the data matrix; discretize continuous variables for the latter
  X.numeric <- data |> apply(2, as.numeric)
  X.binned <- if(identical(X.numeric, round(X.numeric))) X.numeric else{
    catt('Binning continuous data... \r')
    apply(X.numeric, 2, \(x) if (identical(x, round(x))) x else cut.props(x)) |> (`dimnames<-`)(dimnames(X.numeric))
  }
  
  # Calculate core-variable correlation matrix if not given; replace NAs with zeros
  catt('Estimating correlations... \r')
  R.vv <- corr.matrix[v, v] %||% sparse.cor(X.numeric[core.rows, v]) |> na.zeros()
  
  # If loadings not given, run a factor analysis, with arguments passed down
  catt('Factor-analyzing the data... \r')
  P.vf <- loadings[v,] %||% psych.fa(R.vv) # P.vf is the core loadings matrix
  f <- ncol(P.vf) # Number of factors
  
  # Create a matrix of dummy variables and center it; calculate correlation matrix
  catt('Estimating category loadings... \r')
  X.list <- dummy.extension(X.binned[, e], X.numeric[, v], P.vf, R.vv)
  X.ne <- X.list$X.ne
  X.ne[X.ne == 0] <- NA
  R.ee <- X.list$R.ee
  X.ne <- X.ne[ , colnames(X.ne)[colnames(X.ne) %in% rownames(P.ef)] %||% T]

  # Using the Gorsuch factor extension method, find loadings for dummy variables
  P.ef <- P.ef %||% X.list$P.ef |> (`[`)(colnames(X.ne), T)
  
  # If exploratory category re-ordering is requested-- see Keenan (2023)
  if (exploratory.reorder){
    catt('Re-ordering the categories...         \r')
    P.ef <- category.reorder(X.numeric, X.ne, P.vf, P.ef, R.vv, R.ee, f, core.rows, v)
  } 
  U.inv <- diag(1 / (1 - rowSums(P.ef^2)))    # Inverse of uniqueness
  
  # If scoring not requested, just return the core and extension loadings
  if (!score) return(list(scores=NULL, original.loadings=loadings, extension.loadings=P.ef))
  
  # Do row-by-row scoring of the dummy variable matrix, in parallel if requested
  catt('Scoring observations...                    \n')
  ncores <- if (parallel) future::availableWorkers() |> length() else 1
  scores <- score.matrix(X.ne, P.ef, U.inv, R.ee, method, ncores, verbose) |> (`rownames<-`)(rownames(data))
  colnames(scores) <- c(colnames(P.vf), paste(rep('Validity', f), 1:f, sep='_'))

  gc() # Because I can
  
  # Return the estimated scores, original loadings, and extension loadings
  result <- list(scores=scores, original.loadings=P.vf, extension.loadings=P.ef) |> (`class<-`)("dummyscore")
  return(result)
}
