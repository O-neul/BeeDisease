// generated with brms 2.23.0
functions {
 /* Return the log probability of a proper conditional autoregressive (CAR)
  * prior with a sparse representation for the adjacency matrix
  * Full credit to Max Joseph (https://github.com/mbjoseph/CARstan)
  * Args:
  *   phi: Vector containing the CAR parameters for each location
  *   car: Dependence (usually spatial) parameter for the CAR prior
  *   sdcar: Standard deviation parameter for the CAR prior
  *   Nloc: Number of locations
  *   Nedges: Number of edges (adjacency pairs)
  *   Nneigh: Number of neighbors for each location
  *   eigenW: Eigenvalues of D^(-1/2) * W * D^(-1/2)
  *   edges1, edges2: Sparse representation of adjacency matrix
  * Details:
  *   D = Diag(Nneigh)
  * Returns:
  *   Log probability density of CAR prior up to additive constant
  */
  real sparse_car_lpdf(vector phi, real car, real sdcar, int Nloc,
                       int Nedges, data vector Nneigh, data vector eigenW,
                       array[] int edges1, array[] int edges2) {
    real tau;  // precision parameter
    row_vector[Nloc] phit_D;  // phi' * D
    row_vector[Nloc] phit_W;  // phi' * W
    vector[Nloc] ldet;
    tau = inv_square(sdcar);
    phit_D = (phi .* Nneigh)';
    phit_W = rep_row_vector(0, Nloc);
    for (i in 1:Nedges) {
      phit_W[edges1[i]] = phit_W[edges1[i]] + phi[edges2[i]];
      phit_W[edges2[i]] = phit_W[edges2[i]] + phi[edges1[i]];
    }
    for (i in 1:Nloc) {
      ldet[i] = log1m(car * eigenW[i]);
    }
    return 0.5 * (Nloc * log(tau) + sum(ldet) -
           tau * (phit_D * phi - car * (phit_W * phi)));
  }
}
data {
  int<lower=1> N;  // total number of observations
  array[N] int Y;  // response variable
  array[N] int lb;  // lower truncation bounds;
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> Kc;  // number of population-level effects after centering
  // data for the CAR structure
  int<lower=1> Nloc;
  array[N] int<lower=1> Jloc;
  int<lower=0> Nedges;
  array[Nedges] int<lower=1> edges1;
  array[Nedges] int<lower=1> edges2;
  vector[Nloc] Nneigh;
  vector[Nloc] eigenMcar;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  // matrices for QR decomposition
  matrix[N, Kc] XQ;
  matrix[Kc, Kc] XR;
  matrix[Kc, Kc] XR_inv;
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
  // compute and scale QR decomposition
  XQ = qr_thin_Q(Xc) * sqrt(N - 1);
  XR = qr_thin_R(Xc) / sqrt(N - 1);
  XR_inv = inverse(XR);
}
parameters {
  vector[Kc] bQ;  // regression coefficients on QR scale
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> sdcar;  // SD of the CAR structure
  vector[Nloc] rcar;
  real<lower=0,upper=1> car;
  real<lower=0> shape;  // shape parameter
}
transformed parameters {
  // prior contributions to the log posterior
  real lprior = 0;
  lprior += normal_lpdf(bQ | 0, 1);
  lprior += normal_lpdf(Intercept | 1.386294361120, 1.5);
  lprior += student_t_lpdf(sdcar | 3, 0, 2.5)
    - 1 * student_t_lccdf(0 | 3, 0, 2.5);
  lprior += exponential_lpdf(shape | 1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    mu += Intercept + XQ * bQ;
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += rcar[Jloc[n]];
    }
    mu = exp(mu);
    for (n in 1:N) {
      target += neg_binomial_2_lpmf(Y[n] | mu[n], shape) -
        neg_binomial_2_lccdf(lb[n] - 1 | mu[n], shape);
    }
  }
  // priors including constants
  target += lprior;
  target += sparse_car_lpdf(
    rcar | car, sdcar, Nloc, Nedges, Nneigh, eigenMcar, edges1, edges2
  );
}
generated quantities {
  // obtain the actual coefficients
  vector[Kc] b = XR_inv * bQ;
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}

