// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppEigen.h>
#include <memory>
#include <vector>
#include <cmath>
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
typedef Eigen::SimplicialLDLT<SparseMatrix<double> > SparseFactor;

static double one_minus_average(double z) {
  if (z >= .05) return 1.0 + std::expm1(-z)/z;
  double term = z/2.0, value = term;
  for (int k=2; k<=16; ++k) { term *= -z/(k+1.0); value += term; }
  return value;
}

static double exponential_variance(double z) {
  if (z >= .05) {
    double average = -std::expm1(-z)/z;
    return -std::expm1(-2*z)/(2*z) - average*average;
  }
  // Centered 8-point Gauss rule; expm1 retains the O(z) deviations.
  static const double nodes[8] = {
    .0198550717512319,.101666761293187,.237233795041836,.408282678752175,
    .591717321247825,.762766204958164,.898333238706813,.980144928248768};
  static const double weights[8] = {
    .0506142681451881,.111190517226687,.156853322938944,.181341891689181,
    .181341891689181,.156853322938944,.111190517226687,.0506142681451881};
  double shift = one_minus_average(z), value = 0;
  for (int q=0; q<8; ++q) {
    double deviation = std::expm1(-z*nodes[q]) + shift;
    value += weights[q]*deviation*deviation;
  }
  return value;
}

// Every time step applies the shifted FEM matrix solves. No eigenvalues or
// modal recurrence of the discrete operator are used in this routine.
// [[Rcpp::export]]
Rcpp::List rough_matrix_evolve_cpp(
    const Eigen::MappedSparseMatrix<double>& C_in, Rcpp::List matrices,
    Eigen::VectorXd residues, double remainder, Eigen::VectorXd load,
    Eigen::MatrixXd projection, Eigen::MatrixXd residual_gram,
    Eigen::VectorXd rates, Eigen::VectorXd mode_j,
    double tau, int n_steps, double tail_squared, bool return_states = false) {
  SparseMatrix<double> C = C_in;
  int n = C.rows(), m = residues.size(), modes = rates.size();
  if (n_steps < 1 || tau <= 0 || matrices.size()!=m || load.size()!=n ||
      projection.rows()!=n || projection.cols()!=modes || mode_j.size()!=modes ||
      residual_gram.rows()!=modes || residual_gram.cols()!=modes)
    Rcpp::stop("Invalid rough matrix solver dimensions.");
  SparseFactor mass_factor;
  mass_factor.compute(C);
  if (mass_factor.info()!=Eigen::Success) Rcpp::stop("Mass factorization failed.");
  std::vector<std::unique_ptr<SparseFactor> > factors;
  for (int i=0; i<m; ++i) {
    SparseMatrix<double> A = Rcpp::as<SparseMatrix<double> >(matrices[i]);
    factors.emplace_back(new SparseFactor());
    factors.back()->compute(A);
    if (factors.back()->info()!=Eigen::Success) Rcpp::stop("Shifted factorization failed.");
  }
  VectorXd u=VectorXd::Zero(n), rhs(n), next(n), mean(modes), delta(n);
  VectorXd amplitude(modes), average_shift(modes), variance(modes);
  for (int j=0; j<modes; ++j) {
    amplitude[j]=1.0/(mode_j[j]*rates[j]);
    average_shift[j]=one_minus_average(rates[j]*tau);
    variance[j]=exponential_variance(rates[j]*tau);
  }
  MatrixXd states(return_states ? n : 0, return_states ? n_steps+1 : 0);
  if (return_states) states.col(0)=u;
  double total=0, compensation=0;
  for (int k=1; k<=n_steps; ++k) {
    if ((k & 4095)==0) Rcpp::checkUserInterrupt();
    rhs.noalias()=C*u;
    rhs.noalias()+=tau*load;
    next=remainder*mass_factor.solve(rhs);
    for (int i=0; i<m; ++i) next.noalias()+=residues[i]*factors[i]->solve(rhs);
    u.swap(next);
    double variance_sum=0;
    for (int j=0; j<modes; ++j) {
      double z=rates[j]*tau*(k-1.0), decay=std::exp(-z);
      mean[j]=amplitude[j]*(-std::expm1(-z)+decay*average_shift[j]);
      variance_sum+=amplitude[j]*amplitude[j]*decay*decay*variance[j];
    }
    delta.noalias()=u-projection*mean;
    double value=tau*(delta.dot(C*delta)+mean.dot(residual_gram*mean)+variance_sum);
    if (!std::isfinite(value) || value < 0) Rcpp::stop("Invalid squared matrix error.");
    double adjusted=value-compensation, updated=total+adjusted;
    compensation=(updated-total)-adjusted;
    total=updated;
    if (return_states) states.col(k)=u;
  }
  return Rcpp::List::create(Rcpp::Named("error")=std::sqrt(total+tail_squared),
    Rcpp::Named("resolved_squared")=total, Rcpp::Named("tail_squared")=tail_squared,
    Rcpp::Named("final_state")=u, Rcpp::Named("states")=states);
}
