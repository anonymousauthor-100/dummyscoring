#include <RcppParallel.h>
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

class SimplerProgressBar: public ProgressBar{
  public: 
    SimplerProgressBar()  { reset(); }
    ~SimplerProgressBar() {}

  public:

    void display() {
      Rprintf("0%%   10   20   30   40   50   60   70   80   90   100%%\n");
      Rprintf("[----+----+----+----+----+----+----+----+----+----]\n");
      Rprintf("[");
      flush_console();
    }

    // will finalize display if needed
    void update(float progress) {
      _update_ticks_display(progress);
      if (_ticks_displayed >= _max_ticks)
        _finalize_display();
    }

    void end_display() {
      update(1);
      reset();
    }

    void reset() {
      _max_ticks = 49;
      _ticks_displayed = 0;
      _finalized = false;
    }


  protected: // ==== other instance methods =====

    // update the ticks display corresponding to progress
    void _update_ticks_display(float progress) {
      int nb_ticks = _compute_nb_ticks(progress);
      int delta = nb_ticks - _ticks_displayed;
      if (delta > 0) {
        _display_ticks(delta);
        _ticks_displayed = nb_ticks;
      }

    }

    void _finalize_display() {
      if (_finalized) return;

      Rprintf("]\n");
      flush_console();
      _finalized = true;
    }

    int _compute_nb_ticks(float progress) {
      return int(progress * _max_ticks);
    }

    void _display_ticks(int nb) {
      for (int i = 0; i < nb; ++i) {
        Rprintf("=");
        R_FlushConsole();
      }
    }

    // N.B: does nothing on windows
    void flush_console() {
#if !defined(WIN32) && !defined(__WIN32) && !defined(__WIN32__)
       R_FlushConsole();
#endif
    }

  private:
    int _max_ticks;   		// the total number of ticks to print
    int _ticks_displayed; 	// the nb of ticks already displayed
    bool _finalized;

};

// [[Rcpp::export(name = "scale.cpp")]]
NumericMatrix scale(const NumericMatrix &X){

  NumericMatrix X_scaled(X.nrow(), X.ncol());
  for(int i = 0; i < X.ncol(); ++i){
    NumericVector thiscol = X(_, i);
    LogicalVector isfinite = is_finite(thiscol);
    NumericVector finitevals = thiscol[isfinite];
    finitevals = (finitevals - mean(finitevals)) / sd(finitevals);
    thiscol[isfinite] = finitevals;
    X_scaled(_, i) = thiscol;
  }
  return X_scaled;
}

// [[Rcpp::export(name = "sparse.cor")]]
NumericMatrix sparsecor(NumericMatrix &X, NumericMatrix &Y){

  NumericMatrix R(X.ncol(), Y.ncol());
  arma::mat R_arma(R.begin(), R.nrow(), R.ncol(), false);

  NumericMatrix X_scaled = scale(X), Y_scaled = scale(Y);
  arma::mat X_arma(X_scaled.begin(), X.nrow(), X.ncol(), false),
    Y_arma(Y_scaled.begin(), Y.nrow(), Y.ncol(), false);
  X_arma.replace(NA_REAL, 0);
  Y_arma.replace(NA_REAL, 0);
  arma::sp_mat X_sp(X_arma), Y_sp(Y_arma);
  
  arma::sp_mat X_cases = X_sp, Y_cases = Y_sp;
  X_cases.transform([] (double val) {return 1;});
  Y_cases.transform([] (double val) {return 1;});
  arma::mat pairwise_n(X_cases.t() * Y_cases);
  pairwise_n.replace(0, 1);
  
  R_arma = arma::mat(X_sp.t() * Y_sp) / pairwise_n;
  
  return R;
}

// [[Rcpp::export(name = "score.matrix")]]
NumericMatrix score_matrix(arma::mat &X, arma::mat &P, arma::mat &Uinv, arma::mat &R, String method, int nthreads, bool showprogress) {
  
  int n = X.n_rows, c = P.n_cols;
  NumericMatrix score2(n, c*2);
  arma::mat score(score2.begin(), score2.nrow(), score2.ncol(), false);
  
  SimplerProgressBar pb;
  Progress p(n, showprogress, pb);
  
  #ifdef _OPENMP
  #pragma omp parallel for num_threads(nthreads) schedule(dynamic)
  #endif
    
  for(int i=0; i < n; ++i){
  
    p.increment();
    arma::mat Y = X.row(i);
    arma::uvec notna = find_finite(Y);
    arma::mat Z = Y.cols(notna), Uinv2 = Uinv(notna, notna), R2 = R(notna, notna), P2 = P.rows(notna), W = P2;
      
    if(method == "Thurstone")       W = pinv(R2) * P2;
    else if(method == "Ledermann")  W = Uinv2 * P2 * pinv(arma::eye(c,c) + P2.t() * Uinv2 * P2);
    else if(method == "Heermann")   W = pinv(R2) * P2 * pinv(P2.t() * pinv(R2) * P2);
    else if(method == "Bartlett")   W = Uinv2 * P2 * pinv(P2.t() * Uinv2 * P2);
    else if(method == "Harman")     W = P2 * pinv(P2.t() * P2);
    else if(method == "Anderson")   W = Uinv2 * P2 * real(sqrtmat(pinv(P2.t() * Uinv2 * P2)));
    else if(method == "Keenan")     W = P2 * pinv(arma::eye(c,c) + P2.t() * P2);
    
    arma::mat fvar = W.t() * R2 * W, fcov = P2.t() * W;
    arma::vec validity = fcov.diag() / sqrt(fvar.diag());
    score.row(i) = join_rows(Z * W, validity.t());
  }
  
  return score2;
}

// [[Rcpp::export(name = "dummy.extension")]]
List dummy_extension(NumericMatrix &X, NumericMatrix &X_v, NumericMatrix &P_vf, arma::mat &R_vv) {

  int v = X_v.ncol(), e = X.ncol(), n = X_v.nrow();

  arma::mat P_vf_arma(P_vf.begin(), v, P_vf.ncol(), false);
  NumericMatrix X_scaled = scale(X_v);
  arma::mat X_arma(X_scaled.begin(), n, v, false);
  X_arma.replace(NA_REAL, 0);
  arma::sp_mat X_scaled_sp(X_arma);
  arma::sp_mat X_scaled_counts = X_scaled_sp;
  X_scaled_counts.transform([] (double val) {return 1;});
  
  int bignumber = sum(is_finite(X))*10;
  IntegerVector rowindex(bignumber),
    colpointers(e * 10),
    rows = seq_len(X.nrow()) - 1;
  NumericVector scaled_values = rep(1.0, bignumber),
    centered_values = rep(1.0, bignumber),
    count_values = rep(1.0, bignumber);
  
  CharacterVector newcolnames(e * 10), oldcolnames = colnames(X);
  
  int counter = 0, ncols = 0;
  
  for(int i = 0; i < e; ++i){
  
    NumericVector thiscol = X(_,i),
      uniquevals = sort_unique(thiscol);
    uniquevals = uniquevals[!is_na(uniquevals)];
    
    for(int j = 0; j < uniquevals.size(); ++j){
    
      int uniqueval = uniquevals[j];
      LogicalVector isnonzero = is_finite(thiscol);
      isnonzero[is_na(isnonzero)] = false;
      
      NumericVector nonzerovals = thiscol[isnonzero];
      nonzerovals = ifelse(nonzerovals == uniqueval, 1, 0);
      NumericVector nonzerovals_centered = (nonzerovals - mean(nonzerovals)),
        nonzerovals_scaled = nonzerovals_centered / sd(nonzerovals);
      
      int nonzeros = nonzerovals.size();
      IntegerVector range = seq(counter, counter + nonzeros - 1);
      
      rowindex[range] = rows[isnonzero];
      scaled_values[range] = nonzerovals_scaled;
      centered_values[range] = nonzerovals_centered;
      colpointers[ncols + 1] = colpointers[ncols] + nonzeros;
      
      String oldcolname = oldcolnames[i];
      newcolnames[ncols] = oldcolname.push_back("_").push_back(std::to_string(j+1));
      
      counter += nonzeros;
      ncols += 1;
    }
  }
  
  newcolnames = newcolnames[seq_len(ncols) - 1];
  
  NumericMatrix X_ne(n, ncols),
    P_ef(ncols, P_vf.ncol()),
    R_ee(ncols, ncols),
    R_ev(ncols, v);
    
  colnames(X_ne) = newcolnames;
  P_ef.attr("dimnames") = List::create(newcolnames, colnames(P_vf));
    
  arma::mat X_ne_arma(X_ne.begin(), n, ncols, false),
    P_ef_arma(P_ef.begin(), ncols, P_ef.ncol(), false),
    R_ee_arma(R_ee.begin(), ncols, ncols, false),
    R_ev_arma(R_ev.begin(), ncols, v, false);
  
  arma::uvec colptr = as<arma::uvec>(colpointers[seq_len(ncols + 1) - 1]), 
    rowind = as<arma::uvec>(rowindex[seq_len(counter) - 1]);
  arma::colvec scaled_values2 = as<arma::colvec>(scaled_values[seq_len(counter) - 1]),
    centered_values2 = as<arma::colvec>(centered_values[seq_len(counter) - 1]),
    count_values2 = as<arma::colvec>(count_values[seq_len(counter) - 1]);
  
  arma::sp_mat X_ne_sp(rowind, colptr, scaled_values2, n, ncols),
    X_ne_cntr(rowind, colptr, centered_values2, n, ncols),
    counts(rowind, colptr, count_values2, n, ncols);
  
  X_ne_arma = arma::mat(X_ne_cntr) / 1;
  
  arma::mat counts_prod(counts.t() * counts);
  counts_prod.replace(0, 1);
  R_ee_arma = arma::mat(X_ne_sp.t() * X_ne_sp) / counts_prod;
  
  arma::mat R_ev_counts(counts.t() * X_scaled_counts);
  R_ev_counts.replace(0, 1);
  R_ev_arma = arma::mat(X_ne_sp.t() * X_scaled_sp) / R_ev_counts;

  P_ef_arma = R_ev_arma * pinv(R_vv) * P_vf_arma;

  return List::create(_["R.ee"] = R_ee, _["X.ne"] = X_ne, _["R.ev"] = R_ev, _["P.ef"] = P_ef);
}
