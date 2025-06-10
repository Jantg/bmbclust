#include <RcppArmadillo.h>

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]


arma::mat dMvn(arma::mat x,
              arma::rowvec mean,
              arma::mat sigma,
              bool log = false) {
  double log2pi = std::log(2.0 * M_PI);
  int xdim = x.n_cols;
  arma::mat noise = arma::randu(sigma.n_rows,sigma.n_cols);
  //arma::mat noise_symmat = trimatu(noise) + arma::trans(trimatu(noise))-arma::diagmat(noise.diag());
  arma::mat sigma_w_noise = sigma + arma::diagmat(noise.diag()*0.01);
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma_w_noise))));
  double rootisum = arma::sum(arma::log(rooti.diag()));
  double constant = -(static_cast<double>(xdim)/2.0) * log2pi+rootisum;
  arma::mat Z = rooti*arma::trans(x.each_row() - mean);
  arma::mat Kernel = arma::sum(arma::trans(Z%Z),1);
  arma::mat out = -.5*Kernel+constant;
  if (log == false) {
    out = arma::exp(out);
  }
  return out;
}

// [[Rcpp::export]]
arma::mat dMvn_multi(arma::mat x_tmp,
               arma::rowvec mean,
               arma::mat sigma,
               int K,
               bool log = false) {
  double log2pi = std::log(2.0 * M_PI);
  arma::mat multip = arma::ones(1,K);
  arma::mat x = arma::kron(multip,x_tmp);
  int xdim = x_tmp.n_cols;
  int N = x_tmp.n_rows;
  int clst_dim = sigma.n_cols/K;

  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  arma::vec diags = arma::log(arma::vectorise(rooti.diag()));
  arma::mat Z = rooti*arma::trans(x.each_row()-mean);
  arma::mat Kernel_all = arma::trans(Z%Z);
  arma::mat retmat = arma::zeros(N,K);
  
  for(int i=0;i<K;++i){
    arma::mat Kernel_k = arma::sum(Kernel_all.cols(i*clst_dim,((i+1)*clst_dim-1)),1);
    double rootisum_k = arma::sum(diags.subvec(i*clst_dim,((i+1)*clst_dim)-1));
    double constant_k = -(static_cast<double>(xdim)/2.0) * log2pi+rootisum_k;
    arma::mat out = -.5*Kernel_k+constant_k;
    retmat.col(i) = out;
  }
  if (log == false) {
    retmat = arma::exp(retmat);
  }
  return retmat;
}
