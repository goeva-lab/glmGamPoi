#ifndef FIT_BETA_H
#define FIT_BETA_H

#include <calc_helpers.h>
#include <deviance_eigen.h>
#include <fisher_scoring_steps_eigen.h>
#include <opt_max.h>

#include <cfloat>
#include <cmath>

#include <Eigen/Dense>
#include <utility>
template <class D> using EMB = Eigen::MatrixBase<D>;
using Eigen::ArrayXd;
using Eigen::VectorXd;

#include <LBFGSpp/include/LBFGS.h>

template <class D1, class D2>
inline double compute_pen_sum(const EMB<D1> &beta_min_ridge_t, const EMB<D2> &ridge_penalty_sq, const double n_samples) {
  return n_samples * (beta_min_ridge_t.transpose() * ridge_penalty_sq * beta_min_ridge_t)(0, 0);
}

/**
 * This method takes in a proposal for a step and checks if it actually
 * decreases the deviance of the model. If does not add it tries again
 * with half the step size, then a quarter and so on.
 *
 * If even after 100 steps the deviance (0.5^100 = 7.9e-31) has not decreased
 * it returns NaN.
 *
 * Note that the first two parameters are changed: beta_hat and mu_hat
 *
 * The function returns the new deviance.
 *
 */
template <class D1, class D2, class D3, class D4, class D5, class D6, class Fn1>
inline double decrease_deviance(
    // in-out params
    EMB<D1> &beta_hat, EMB<D2> &mu_hat,
    // const params
    const EMB<D3> &step, const EMB<D4> &model_matrix,
    // callback for computing gp deviance
    const Fn1 &cgp_sum,
    // more const params
    const EMB<D5> &exp_off, const EMB<D6> &counts,
    // misc params
    const double theta, const double dev_old, const double tolerance, const double max_rel_mu_change) {
  double speeding_factor = 1.0;
  int line_iter = 0;
  double dev = 0;
  beta_hat += step;
  const VectorXd mu_old = mu_hat;
  while (true) {
    mu_hat = calculate_mu<D2>(model_matrix, beta_hat, exp_off);
    dev = cgp_sum(counts, mu_hat, theta, beta_hat);
    const double conv_test = std::fabs(dev - dev_old) / (std::fabs(dev) + 0.1);
    const double mu_rel_change = mu_hat.cwiseQuotient(mu_old).maxCoeff();
    if ((dev < dev_old && mu_rel_change < max_rel_mu_change) || conv_test < tolerance) {
      break; // while loop
    } else if (line_iter >= 100) {
      // speeding factor is very small, something is going wrong here
      dev = std::numeric_limits<double>::quiet_NaN();
      break; // while loop
    } else {
      // Halfing the speed
      speeding_factor = speeding_factor / 2.0;
      beta_hat = beta_hat - (step * speeding_factor);
    }
    line_iter++;
  }
  return dev;
}

template <class D1, class D2, class D3, class Fn1> class BetaOptim {
private:
  const Fn1 &cgp_sum_;
  const EMB<D1> &model_matrix_;
  const EMB<D2> &counts_;
  const EMB<D3> &exp_off_;
  const double &theta_;
  inline double fn_(const VectorXd &beta_row) const {
    return cgp_sum_(counts_, calculate_mu<VectorXd>(model_matrix_, beta_row, exp_off_), theta_, beta_row);
  }

public:
  BetaOptim(const Fn1 &cgp_sum, const EMB<D1> &model_matrix, const EMB<D2> &counts, const EMB<D3> &exp_off, const double &theta)
      : cgp_sum_(cgp_sum), model_matrix_(model_matrix), counts_(counts), exp_off_(exp_off), theta_(theta) {};
  inline double operator()(const VectorXd &beta_row, VectorXd &grad) const {
    const auto out = fn_(beta_row);

    // we do not know the gradient here explicitly, approximating instead using a numerical method
    // we need to copy the beta values since the provided reference is const
    // the eps value (1e-3) is taken from the default value used by R for this algorithm
    VectorXd beta_c = beta_row;
    auto beta_begin = beta_c.begin();
    const auto beta_end = beta_c.end();
    auto grad_begin = grad.begin();

    for (; beta_begin != beta_end; ++beta_begin, ++grad_begin) {
      const double b_curr = *beta_begin;
      *beta_begin = b_curr + 1e-3;
      const double v1 = fn_(beta_c);
      *beta_begin = b_curr - 1e-3;
      const double v2 = fn_(beta_c);
      *grad_begin = (v1 - v2) / 2e-3;
      *beta_begin = b_curr;
    }

    return out;
  }
};

const double RELTOL = std::sqrt(DBL_EPSILON);
template <class D1, class D2, class D3, class Fn1>
inline void fitBeta_FS_optim_step(
    // in-out params
    VectorXd &beta_out, double &dev_out, int &iters_out,
    // const params
    const EMB<D1> &model_matrix, const EMB<D2> &counts, const EMB<D3> &exp_off,
    // callback for optim
    const Fn1 &cgp_sum,
    // other params
    const double theta, const int max_iter) {
  LBFGSpp::LBFGSParam params;
  params.max_iterations = max_iter;
  params.past = 1;
  params.delta = RELTOL;

  LBFGSpp::LBFGSSolver<double, LBFGSpp::LineSearchMoreThuente> solver(params);
  double fx;
  const BetaOptim fn(cgp_sum, model_matrix, counts, exp_off, theta);

  iters_out = solver.minimize(fn, beta_out, fx);
  if (iters_out == max_iter) {
    beta_out.fill(NAN);
    dev_out = NAN;
    return;
  }
  dev_out = fx;
}

template <class D1, class D2, class D3, class D4, class Fn1, class Fn2>
inline void fitBeta_FS_internal_step(
    // in-out params
    EMB<D1> &beta_out, double &dev_out, int &iters_out,
    // const params
    const EMB<D2> &model_matrix, const EMB<D3> &counts, const EMB<D4> &exp_off,
    // callbacks for steps
    const Fn1 &cgp_sum, const Fn2 &fisher_step,
    // other params
    const double theta, const double tolerance, const double max_rel_mu_change, const int max_iter, const bool try_recov_w_optim) {

  // Init beta and mu
  VectorXd beta_hat = beta_out;
  VectorXd mu_hat = calculate_mu<VectorXd>(model_matrix, beta_hat, exp_off);

  if (beta_hat.array().isNaN().any() || std::isnan(theta)) {
    beta_out.fill(NAN);
    dev_out = NAN;
    iters_out = 0;
    return;
  }

  // Init deviance
  double dev_old = cgp_sum(counts, mu_hat, theta, beta_hat);
  int iter = 0;
  for (; iter < max_iter; iter++) {
    // Find good direction to optimize beta
    const VectorXd step = fisher_step(model_matrix, counts, mu_hat, theta * mu_hat, beta_hat);
    // Find step size that actually decreases the deviance
    const double dev =
        decrease_deviance(beta_hat, mu_hat, step, model_matrix, cgp_sum, exp_off, counts, theta, dev_old, tolerance, max_rel_mu_change);

    const double conv_test = std::fabs(dev - dev_old) / (std::fabs(dev) + 0.1);
    dev_old = dev;
    if (std::isnan(conv_test)) {
      // This should not happen
      beta_hat.fill(NAN);
      iter = max_iter;
      break;
    }
    if (conv_test < tolerance) {
      break;
    }
  }

  if (try_recov_w_optim && iter == max_iter) {
    beta_hat = beta_out; // re-copy since from reference since values have been overwritten
    fitBeta_FS_optim_step(beta_hat, dev_old, iter, model_matrix, counts, exp_off, cgp_sum, theta, max_iter);
  }

  beta_out = beta_hat;
  dev_out = dev_old;
  iters_out = iter;
}

template <class D1, class D2> inline double fitBeta_NR_optim_step(const EMB<D1> &counts, const EMB<D2> &off, const double theta) {
  return optimize_fmax([&counts = std::as_const(counts), &off = std::as_const(off), &theta = std::as_const(theta)](double beta) -> double {
    double out = 0.0;
    // iterate through counts & offset simultaneously
    auto counts_begin = counts.cbegin();
    auto off_begin = off.cbegin();
    const auto counts_end = counts.cend();
    // we known that counts and off are of the same size, so we only check one vector
    for (; counts_begin != counts_end; ++counts_begin, ++off_begin) {
      out += dnbinom_impl(*counts_begin, beta, *off_begin, theta);
    }
    return out;
  });
}

template <class D1, class D2>
inline void fitBeta_NR_internal_step(
    // in-out parameters
    double &beta_out, double &dev_out, int &iters_out,
    // const params
    const EMB<D1> &counts, const EMB<D2> &off,
    // other params
    const double theta, const double tolerance, const int max_iter) {

  double beta = beta_out;
  if (std::isnan(beta) || std::isnan(theta)) {
    beta_out = NAN;
    dev_out = NAN;
    iters_out = 0;
    return;
  }
  if (counts.isApproxToConstant(0.0)) {
    beta_out = -INFINITY;
    dev_out = 0;
    iters_out = 0;
    return;
  }

  // Newton-Raphson
  int iter = 0;
  for (; iter < max_iter; iter++) {

    const ArrayXd mu = (off.array() + beta).exp();
    const ArrayXd denom = 1.0 + (mu * theta).array();

    const double dl = ((counts.array() - mu) / denom).sum();
    const double ddl = (mu * (1.0 + counts.array() * theta) / denom / denom).sum();

    const double step = dl / ddl;
    beta += step;
    if ((std::abs(step) < tolerance) || (std::isnan(beta))) {
      break;
    }
  }

  if (iter == max_iter || std::isnan(beta)) {
    beta = fitBeta_NR_optim_step(counts, off, theta);
  }

  beta_out = beta;
  dev_out = compute_gp_deviance_sum(counts, (off.array() + beta).exp().matrix(), theta);
  iters_out = iter;
}

#endif