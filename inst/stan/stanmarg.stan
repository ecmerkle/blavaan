/* This file is based on LERSIL.stan by Ben Goodrich.
   https://github.com/bgoodri/LERSIL */
functions { // you can use these in R following `rstan::expose_stan_functions("foo.stan")`
  /*
    Fills in the elements of a coefficient matrix containing some mix of 
    totally free, free subject to a sign constraint, and fixed elements
    
    @param free_elements vector of unconstrained elements
    @param skeleton matrix of the same dimensions as the output whose elements are
      positive_infinity(): if output element is totally free
      other: if output element is fixed to that number
    @return matrix of coefficients
  */
  matrix fill_matrix(vector free_elements, matrix skeleton, int[,] eq_skeleton, int pos_start, int spos_start) {
    int R = rows(skeleton);
    int C = cols(skeleton);
    matrix[R, C] out;

    int pos = spos_start; // position of eq_skeleton
    int freepos = pos_start; // position of free_elements
    int eqelem = 0;
    
    for (c in 1:C) for (r in 1:R) {
      real rc = skeleton[r, c];
      if (is_inf(rc)) { // free
	int eq = eq_skeleton[pos, 1];
	int wig = eq_skeleton[pos, 3];
	if (eq == 0 || wig == 1) {
	  out[r,c] = free_elements[freepos];
	  freepos += 1;
	} else {
	  eqelem = eq_skeleton[pos, 2];
	  out[r,c] = free_elements[eqelem];
	}
	pos += 1;
      } else out[r,c] = skeleton[r, c]; // fixed, so do not bump pos
    }
    return out;
  }
    
  vector fill_prior(vector free_elements, real[] pri_mean, int[,] eq_skeleton) {
    int R = dims(eq_skeleton)[1];
    int eqelem = 0;
    int pos = 1;
    vector[num_elements(pri_mean)] out;

    for (r in 1:R) {
      if (pos <= num_elements(pri_mean)) {
	int eq = eq_skeleton[r, 1];
	int wig = eq_skeleton[r, 3];

	if (eq == 0) {
	  out[pos] = pri_mean[pos];
	  pos += 1;
	} else if (wig == 1) {
	  eqelem = eq_skeleton[r, 2];
	  out[pos] = free_elements[eqelem];
	  pos += 1;
	}
      }
    }
    return out;
  }
  
  /*
   * This is a bug-free version of csr_to_dense_matrix and has the same arguments
   */
  matrix to_dense_matrix(int m, int n, vector w, int[] v, int[] u) {
    matrix[m, n] out = rep_matrix(0, m, n);
    int pos = 1;
    for (i in 1:m) {
      int start = u[i];
      int nnz = u[i + 1] - start;
      for (j in 1:nnz) {
        out[i, v[pos]] = w[pos];
        pos += 1;
      }
    }
    return out;
  }

  // sign function
  int sign(real x) {
    if (x > 0)
      return 1;
    else
      return -1;
  }

  // sign-constrain a vector of loadings
  vector sign_constrain_load(vector free_elements, int npar, int[,] sign_mat) {
    vector[npar] out;
    for (i in 1:npar) {
      if (sign_mat[i,1]) {
        int lookupval = sign_mat[i,2];
        if (free_elements[lookupval] < 0) {
	  out[i] = -free_elements[i];
	} else {
	  out[i] = free_elements[i];
	}
      } else {
        out[i] = free_elements[i];
      }
    }
    return out;
  }

  // sign-constrain a vector of regressions or covariances
  vector sign_constrain_reg(vector free_elements, int npar, int[,] sign_mat, vector load_par1, vector load_par2) {
    vector[npar] out;
    for (i in 1:npar) {
      if (sign_mat[i,1]) {
        int lookupval1 = sign_mat[i,2];
	int lookupval2 = sign_mat[i,3];
        if (sign(load_par1[lookupval1]) * sign(load_par2[lookupval2]) < 0) {
	  out[i] = -free_elements[i];
	} else {
	  out[i] = free_elements[i];
	}
      } else {
        out[i] = free_elements[i];
      }
    }
    return out;
  }

  // obtain covariance parameter vector for correlation/sd matrices
  vector cor2cov(matrix[] cormat, matrix[] sdmat, int num_free_elements, matrix[] matskel, int[,] wskel, int ngrp) {
    vector[num_free_elements] out;
    int R = rows(to_matrix(cormat[1]));
    int C = cols(to_matrix(cormat[1]));
    int pos = 1; // position of eq_skeleton
    int freepos = 1; // position of free_elements
    
    for (g in 1:ngrp) {
      for (c in 1:(R-1)) for (r in (c+1):R) {
        if (is_inf(matskel[g,r,c])) {
	  if (wskel[pos,1] == 0) {
	    out[freepos] = sdmat[g,r,r] * sdmat[g,c,c] * cormat[g,r,c];
	    freepos += 1;
	  }
	  pos += 1;
	}
      }
    }
    return out;
  }

  // E step of EM algorithm on latent continuous space
  matrix[] estep(vector[] YXstar, vector[] Mu, matrix[] Sigma, int[] Nobs, int[,] Obsvar, int[] startrow, int[] endrow, int[] grpnum, int Np, int Ng) {
    int p = dims(YXstar)[2];
    matrix[p, p + 1] out[Ng]; //mean vec + cov mat
    matrix[dims(YXstar)[1], p] YXfull; // columns consistenly ordered
    matrix[p, p] T2pat;
    int obsidx[p];
    int r1;
    int r2;
    int grpidx;
    int Nmis;

    for (g in 1:Ng) {
      out[g] = rep_matrix(0, p, p + 1);
    }

    for (mm in 1:Np) {
      obsidx = Obsvar[mm,];
      r1 = startrow[mm];
      r2 = endrow[mm];
      grpidx = grpnum[mm];
      Nmis = p - Nobs[mm];

      if (Nobs[mm] < p) {
	matrix[Nobs[mm], Nobs[mm]] Sig22 = Sigma[grpidx, obsidx[1:Nobs[mm]], obsidx[1:Nobs[mm]]];
	matrix[Nmis, Nmis] Sig11 = Sigma[grpidx, obsidx[(Nobs[mm] + 1):p], obsidx[(Nobs[mm] + 1):p]];
	matrix[Nmis, Nobs[mm]] Sig12 = Sigma[grpidx, obsidx[(Nobs[mm] + 1):p], obsidx[1:Nobs[mm]]];
	matrix[Nobs[mm], Nobs[mm]] S22inv = inverse_spd(Sig22);
	matrix[Nmis, Nmis] T2p11 = Sig11 - (Sig12 * S22inv * Sig12');
	
        // partition into observed/missing, compute Sigmas, add to out
	for (jj in r1:r2) {
	  vector[Nmis] ymis;
	  ymis = Mu[grpidx, obsidx[(Nobs[mm] + 1):p]] + (Sig12 * S22inv * (YXstar[jj, 1:Nobs[mm]] - Mu[grpidx, obsidx[1:Nobs[mm]]]));
	  for (kk in 1:Nobs[mm]) {
	    YXfull[jj, obsidx[kk]] = YXstar[jj, kk];
	  }
	  for (kk in (Nobs[mm] + 1):p) {
	    YXfull[jj, obsidx[kk]] = ymis[kk - Nobs[mm]];
	  }
	}
	T2pat = crossprod(YXfull[r1:r2,]);
	// correction for missing cells/conditional covariances
	for (jj in 1:Nmis) {
	  for (kk in jj:Nmis) {
	    T2pat[obsidx[Nobs[mm] + jj], obsidx[Nobs[mm] + kk]] = T2pat[obsidx[Nobs[mm] + jj], obsidx[Nobs[mm] + kk]] + (r2 - r1 + 1) * T2p11[jj, kk];
	    if (kk > jj) {
	      T2pat[obsidx[Nobs[mm] + kk], obsidx[Nobs[mm] + jj]] = T2pat[obsidx[Nobs[mm] + jj], obsidx[Nobs[mm] + kk]];
	    }
	  }
	}
      } else {
	// complete data
	for (jj in r1:r2) {
	  for (kk in 1:Nobs[mm]) {
	    YXfull[jj, obsidx[kk]] = YXstar[jj, kk];
	  }
	}
	T2pat = crossprod(YXfull[r1:r2,]);
      }
      for (i in 1:p) {
	out[grpidx,i,1] += sum(YXfull[r1:r2,i]);
      }
      out[grpidx,,2:(p+1)] += T2pat;
    }
    
    return out;
  }

  matrix sig_inv_update(matrix Sigmainv, int[] obsidx, int Nobs, int np, real logdet) {
    matrix[Nobs + 1, Nobs + 1] out = rep_matrix(0, Nobs + 1, Nobs + 1);
    int nrm = np - Nobs;
    matrix[nrm, nrm] H;
    matrix[nrm, Nobs] A;

    if (nrm == 0) {
      out[1:Nobs, 1:Nobs] = Sigmainv;
      out[Nobs + 1, Nobs + 1] = logdet;
    } else {
      H = Sigmainv[obsidx[(Nobs + 1):np], obsidx[(Nobs + 1):np]];
      A = Sigmainv[obsidx[(Nobs + 1):np], obsidx[1:Nobs]];

      out[1:Nobs, 1:Nobs] = Sigmainv[obsidx[1:Nobs], obsidx[1:Nobs]] - A' * mdivide_left_spd(H, A);
      out[Nobs + 1, Nobs + 1] = logdet + log_determinant(H);
    }

    return out;
  }
  
  real multi_normal_suff(vector xbar, matrix S, vector Mu, matrix Supdate, int N) {
    int Nobs = dims(S)[1];
    real out;

    // using elementwise multiplication + sum here for efficiency
    out = -.5 * N * ( sum(Supdate[1:Nobs, 1:Nobs] .* (S + (xbar - Mu) * (xbar - Mu)')) + Supdate[Nobs + 1, Nobs + 1] + Nobs * log(2 * pi()) );

    return(out);
  }
}
data {
  // see https://books.google.com/books?id=9AC-s50RjacC&lpg=PP1&dq=LISREL&pg=PA2#v=onepage&q=LISREL&f=false
  int<lower=0> p; // number of manifest response variables
  int<lower=0> q; // number of manifest predictors
  int<lower=0> m; // number of latent endogenous variables
  int<lower=0> n; // number of latent exogenous variables
  int<lower=1> Ng; // number of groups
  int<lower=0, upper=1> missing; // are there missing values?
  int<lower=0, upper=1> save_lvs; // should we save lvs?
  int<lower=1> Np; // number of group-by-missing patterns combos
  int<lower=1> N[Ng]; // number of observations per group
  int<lower=1> Nobs[Np]; // number of observed variables in each missing pattern
  int<lower=0> Nordobs[Np]; // number of ordinal observed variables in each missing pattern
  int<lower=0> Obsvar[Np, p + q]; // indexing of observed variables
  int<lower=1> Ntot; // number of observations across all groups
  int<lower=1> startrow[Np]; // starting row for each missing pattern
  int<lower=1,upper=Ntot> endrow[Np]; // ending row for each missing pattern
  int<lower=1,upper=Ng> grpnum[Np]; // group number for each row of data
  int<lower=0,upper=1> wigind; // do any parameters have approx equality constraint ('wiggle')?
  int<lower=0, upper=1> has_data; // are the raw data on y and x available?
  int<lower=0, upper=1> ord; // are there any ordinal variables?
  int<lower=0> Nord; // how many ordinal variables?
  int<lower=0> ordidx[Nord]; // indexing of ordinal variables
  int<lower=0> OrdObsvar[Np, Nord]; // indexing of observed ordinal variables in YXo
  int<lower=0> Noent; // how many observed entries of ordinal variables (for data augmentation)
  int<lower=0> contidx[p + q - Nord]; // indexing of continuous variables
  int<lower=1> nlevs[Nord]; // how many levels does each ordinal variable have
  vector[p + q - Nord] YX[Ntot]; // continuous data
  int YXo[Ntot, Nord]; // ordinal data
  int<lower=0> Nx[Np]; // number of fixed.x variables
  int<lower=0> Xvar[Np, p + q]; // indexing of fixed.x variables
  int<lower=0> Xdatvar[Np, p + q]; // indexing of fixed.x in data (differs from Xvar when missing)
  int<lower=0, upper=1> use_cov;
  int<lower=0, upper=1> pri_only; // sample only from the prior?
  int<lower=0> emiter; // number of em iterations for saturated model in ppp (missing data only)
  int<lower=0, upper=1> use_suff; // should we compute likelihood via mvn sufficient stats?
  int<lower=0, upper=1> do_test; // should we do everything in generated quantities?
  vector[p + q - Nord] YXbar[Np]; // sample means of continuous manifest variables
  matrix[p + q - Nord + 1, p + q - Nord + 1] S[Np];     // sample covariance matrix among all continuous manifest variables NB!! multiply by (N-1) to use wishart lpdf!!

  
  /* sparse matrix representations of skeletons of coefficient matrices, 
     which is not that interesting but necessary because you cannot pass
     missing values into the data block of a Stan program from R */
  int<lower=0> len_w1;        // max number of free elements in Lambda_y per grp
  int<lower=0> wg1[Ng];           // number of free elements in Lambda_y per grp
  vector[len_w1] w1[Ng];          // values of free elements in Lambda_y
  int<lower=1> v1[Ng, len_w1];    // index  of free elements in Lambda_y
  int<lower=1> u1[Ng, p + 1];     // index  of free elements in Lambda_y
  int<lower=0> w1skel[sum(wg1), 3];
  int<lower=0> lam_y_sign[sum(wg1), 2];
  int<lower=0> len_lam_y;     // number of free elements minus equality constraints
  real lambda_y_mn[len_lam_y];           // prior
  real<lower=0> lambda_y_sd[len_lam_y];

  // same things but for Gamma
  int<lower=0> len_w3;
  int<lower=0> wg3[Ng];
  vector[len_w3] w3[Ng];
  int<lower=1> v3[Ng, len_w3];
  int<lower=1> u3[Ng, m + 1];
  int<lower=0> w3skel[sum(wg3), 3];
  int<lower=0> gam_sign[sum(wg3), 3];
  int<lower=0> len_gam;
  real gamma_mn[len_gam];
  real<lower=0> gamma_sd[len_gam];
  
  // same things but for B
  int<lower=0> len_w4;
  int<lower=0> wg4[Ng];
  vector[len_w4] w4[Ng];
  int<lower=1> v4[Ng, len_w4];
  int<lower=1> u4[Ng, m + 1];
  int<lower=0> w4skel[sum(wg4), 3];
  int<lower=0> b_sign[sum(wg4), 3];
  int<lower=0> len_b;
  real b_mn[len_b];
  real<lower=0> b_sd[len_b];
  
  // same things but for diag(Theta)
  int<lower=0> len_w5;
  int<lower=0> wg5[Ng];
  vector[len_w5] w5[Ng];
  int<lower=1> v5[Ng, len_w5];
  int<lower=1> u5[Ng, p + 1];
  int<lower=0> w5skel[sum(wg5), 3];
  int<lower=0> len_thet_sd;
  real<lower=0> theta_sd_shape[len_thet_sd];
  real<lower=0> theta_sd_rate[len_thet_sd];
  int<lower=-2, upper=2> theta_pow;
  int<lower=0> w5use;
  int<lower=1> usethet[w5use];

  // same things but for Theta_r
  int<lower=0> len_w7;
  int<lower=0> wg7[Ng];
  vector[len_w7] w7[Ng];
  int<lower=1> v7[Ng, len_w7];
  int<lower=1> u7[Ng, p + 1];
  int<lower=0> w7skel[sum(wg7), 3];
  int<lower=0> len_thet_r;
  real<lower=0> theta_r_alpha[len_thet_r];
  real<lower=0> theta_r_beta[len_thet_r];
  
  // same things but for Theta_r_x
  int<lower=0> len_w8;
  int<lower=0> wg8[Ng];
  vector[len_w8] w8[Ng];
  int<lower=1> v8[Ng, len_w8];
  int<lower=1> u8[Ng, q + 1];
  int<lower=0> w8skel[sum(wg8), 3];
  int<lower=0> len_thet_x_r;
  real<lower=0> theta_x_r_alpha[len_thet_x_r];
  real<lower=0> theta_x_r_beta[len_thet_x_r];
  
  // same things but for Psi
  int<lower=0> len_w9;
  int<lower=0> wg9[Ng];
  vector[len_w9] w9[Ng];
  int<lower=1> v9[Ng, len_w9];
  int<lower=1> u9[Ng, m + 1];
  int<lower=0> w9skel[sum(wg9), 3];
  int<lower=0> len_psi_sd;
  real<lower=0> psi_sd_shape[len_psi_sd];
  real<lower=0> psi_sd_rate[len_psi_sd];
  int<lower=-2,upper=2> psi_pow;
  int<lower=0> w9use;
  int<lower=1> usepsi[w9use];
  int<lower=0> w9no;
  int<lower=1> nopsi[w9no];
  
  // same things but for Psi_r
  int<lower=0> len_w10;
  int<lower=0> wg10[Ng];
  vector[len_w10] w10[Ng];
  int<lower=1> v10[Ng, len_w10];
  int<lower=1> u10[Ng, m + 1];
  int<lower=0> w10skel[sum(wg10), 3];
  int<lower=0> psi_r_sign[sum(wg10), 3];
  int<lower=0> len_psi_r;
  real<lower=0> psi_r_alpha[len_psi_r];
  real<lower=0> psi_r_beta[len_psi_r];
  int<lower=0,upper=1> fullpsi;
    
  // same things but for Nu
  int<lower=0> len_w13;
  int<lower=0> wg13[Ng];
  vector[len_w13] w13[Ng];
  int<lower=1> v13[Ng, len_w13];
  int<lower=1> u13[Ng, use_cov ? 1 : p + q + 1];
  int<lower=0> w13skel[sum(wg13), 3];
  int<lower=0> len_nu;
  real nu_mn[len_nu];
  real<lower=0> nu_sd[len_nu];
  
  // same things but for Alpha
  int<lower=0> len_w14;
  int<lower=0> wg14[Ng];
  vector[len_w14] w14[Ng];
  int<lower=0> v14[Ng, len_w14];
  int<lower=1> u14[Ng, use_cov ? 1 : m + n + 1];
  int<lower=0> w14skel[sum(wg14), 3];
  int<lower=0> len_alph;
  real alpha_mn[len_alph];
  real<lower=0> alpha_sd[len_alph];

  // same things but for Tau
  int<lower=0> len_w15;
  int<lower=0> wg15[Ng];
  vector[len_w15] w15[Ng];
  int<lower=0> v15[Ng, len_w15];
  int<lower=1> u15[Ng, sum(nlevs) - Nord + 1];
  int<lower=0> w15skel[sum(wg15), 3];
  int<lower=0> len_tau;
  real tau_mn[len_tau];
  real<lower=0> tau_sd[len_tau];
}
transformed data { // (re)construct skeleton matrices in Stan (not that interesting)
  matrix[p, m] Lambda_y_skeleton[Ng];
  //matrix[q, n] Lambda_x_skeleton[Ng];
  matrix[m, n] Gamma_skeleton[Ng];
  matrix[m, m] B_skeleton[Ng];
  matrix[p, p] Theta_skeleton[Ng];
  matrix[p, p] Theta_r_skeleton[Ng];
  matrix[m, m] Psi_skeleton[Ng];
  matrix[m, m] Psi_r_skeleton[Ng];
  //matrix[n, n] Phi_skeleton[Ng];
  //matrix[n, n] Phi_r_skeleton[Ng];
  matrix[(p + q), 1] Nu_skeleton[Ng];
  matrix[(m + n), 1] Alpha_skeleton[Ng];
  matrix[sum(nlevs) - Nord, 1] Tau_skeleton[Ng];
  vector[ord ? 0 : (p + q)] YXbarstar[Np];
  matrix[ord ? 0 : (p + q), ord ? 0 : (p + q)] Sstar[Np];

  matrix[m, m] I = diag_matrix(rep_vector(1, m));

  int Ncont = p + q - Nord;
  
  int g_start1[Ng];
  int g_start2[Ng];
  int g_start3[Ng];
  int g_start4[Ng];
  int g_start5[Ng];
  int g_start6[Ng];
  int g_start7[Ng];
  int g_start8[Ng];
  int g_start9[Ng];
  int g_start10[Ng];
  int g_start11[Ng];
  int g_start12[Ng];
  int g_start13[Ng];
  int g_start14[Ng];
  int g_start15[Ng];
  
  int f_start1[Ng];
  int f_start2[Ng];
  int f_start3[Ng];
  int f_start4[Ng];
  int f_start5[Ng];
  int f_start6[Ng];
  int f_start7[Ng];
  int f_start8[Ng];
  int f_start9[Ng];
  int f_start10[Ng];
  int f_start11[Ng];
  int f_start12[Ng];
  int f_start13[Ng];
  int f_start14[Ng];
  int f_start15[Ng];
  
  int len_free[15];
  int pos[15];
  
  for (i in 1:15) {
    len_free[i] = 0;
    pos[i] = 1;
  }
  
  for (g in 1:Ng) {
    Lambda_y_skeleton[g] = to_dense_matrix(p, m, w1[g], v1[g,], u1[g,]);
    Gamma_skeleton[g] = to_dense_matrix(m, n, w3[g], v3[g,], u3[g,]);
    B_skeleton[g] = to_dense_matrix(m, m, w4[g], v4[g,], u4[g,]);
    Theta_skeleton[g] = to_dense_matrix(p, p, w5[g], v5[g,], u5[g,]);
    Theta_r_skeleton[g] = to_dense_matrix(p, p, w7[g], v7[g,], u7[g,]);
    Psi_skeleton[g] = to_dense_matrix(m, m, w9[g], v9[g,], u9[g,]);
    Psi_r_skeleton[g] = to_dense_matrix(m, m, w10[g], v10[g,], u10[g,]);
    if (!use_cov) {
      Nu_skeleton[g] = to_dense_matrix((p + q), 1, w13[g], v13[g,], u13[g,]);
      Alpha_skeleton[g] = to_dense_matrix((m + n), 1, w14[g], v14[g,], u14[g,]);
    }
    Tau_skeleton[g] = to_dense_matrix(sum(nlevs) - Nord, 1, w15[g], v15[g,], u15[g,]);
    
    // count free elements in Lambda_y_skeleton
    g_start1[g] = len_free[1] + 1;
    f_start1[g] = pos[1];
    for (i in 1:p) {
      for (j in 1:m) {
        if (is_inf(Lambda_y_skeleton[g,i,j])) {
	  if (w1skel[pos[1],2] == 0 || w1skel[pos[1],3] == 1) len_free[1] += 1;
	  pos[1] += 1;
        }
      }
    }

    // same thing but for Gamma_skeleton
    g_start3[g] = len_free[3] + 1;
    f_start3[g] = pos[3];
    for (i in 1:m) {
      for (j in 1:n) {
	if (is_inf(Gamma_skeleton[g,i,j])) {
	  if (w3skel[pos[3],2] == 0 || w3skel[pos[3],3] == 1) len_free[3] += 1;
	  pos[3] += 1;
	}
      }
    }

    // same thing but for B_skeleton
    g_start4[g] = len_free[4] + 1;
    f_start4[g] = pos[4];
    for (i in 1:m) {
      for (j in 1:m) {
	if (is_inf(B_skeleton[g,i,j])) {
	  if (w4skel[pos[4],2] == 0 || w4skel[pos[4],3] == 1) len_free[4] += 1;
	  pos[4] += 1;
	}
      }
    }
    
    // same thing but for Theta_skeleton
    g_start5[g] = len_free[5] + 1;
    f_start5[g] = pos[5];
    for (i in 1:p) {
      if (is_inf(Theta_skeleton[g,i,i])) {
	if (w5skel[pos[5],2] == 0 || w5skel[pos[5],3] == 1) len_free[5] += 1;
	pos[5] += 1;
      }
    }

    // same thing but for Theta_r_skeleton
    g_start7[g] = len_free[7] + 1;
    f_start7[g] = pos[7];
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
	if (is_inf(Theta_r_skeleton[g,j,i])) {
	  if (w7skel[pos[7],2] == 0 || w7skel[pos[7],3] == 1) len_free[7] += 1;
	  pos[7] += 1;
	}
      }
    }

    // same thing but for Psi_skeleton
    g_start9[g] = len_free[9] + 1;
    f_start9[g] = pos[9];
    for (i in 1:m) {
      if (is_inf(Psi_skeleton[g,i,i])) {
	if (w9skel[pos[9],2] == 0 || w9skel[pos[9],3] == 1) len_free[9] += 1;
	pos[9] += 1;
      }
    }

    // same thing but for Psi_r_skeleton
    g_start10[g] = len_free[10] + 1;
    f_start10[g] = pos[10];
    for (i in 1:(m-1)) {
      for (j in (i+1):m) {
	if (is_inf(Psi_r_skeleton[g,j,i])) {
	  if (w10skel[pos[10],2] == 0 || w10skel[pos[10],3] == 1) len_free[10] += 1;
	  pos[10] += 1;
	}
      }
    }

    if (!use_cov) {
      // same thing but for Nu_skeleton
      // pos = len_free13 + 1;
      g_start13[g] = len_free[13] + 1;
      f_start13[g] = pos[13];
      for (i in 1:(p+q)) {
	if (is_inf(Nu_skeleton[g,i,1])) {
	  if (w13skel[pos[13],2] == 0 || w13skel[pos[13],3] == 1) len_free[13] += 1;
	  pos[13] += 1;
	}
      }

      // same thing but for Alpha_skeleton
      g_start14[g] = len_free[14] + 1;
      f_start14[g] = pos[14];
      for (i in 1:(m+n)) {
	if (is_inf(Alpha_skeleton[g,i,1])) {
	  if (w14skel[pos[14],2] == 0 || w14skel[pos[14],3] == 1) len_free[14] += 1;
	  pos[14] += 1;
	}
      }
    }

    // same thing but for Tau_skeleton
    g_start15[g] = len_free[15] + 1;
    f_start15[g] = pos[15];
    for (i in 1:(sum(nlevs) - Nord)) {
      if (is_inf(Tau_skeleton[g,i,1])) {
	if (w15skel[pos[15],2] == 0 || w15skel[pos[15],3] == 1) len_free[15] += 1;
	pos[15] += 1;
      }
    }    
  }

  if (!ord) {
    // sufficient stat matrices by pattern, moved to left for missing
    for (patt in 1:Np) {
      Sstar[patt] = rep_matrix(0, p + q, p + q);
      Sstar[patt, 1:Nobs[patt], 1:Nobs[patt]] = S[patt, Obsvar[patt, 1:Nobs[patt]], Obsvar[patt, 1:Nobs[patt]]];

      for (j in 1:Nobs[patt]) {
	YXbarstar[patt,j] = YXbar[patt, Obsvar[patt,j]];
      }
    }
  }
}
parameters {
  // free elements (possibly with inequality constraints) for coefficient matrices
  vector[len_free[1]] Lambda_y_free;
  //vector[len_free[2]] Lambda_x_free;
  vector[len_free[3]] Gamma_free;
  vector[len_free[4]] B_free;
  vector<lower=0>[len_free[5]] Theta_sd_free;
  vector<lower=-1,upper=1>[len_free[7]] Theta_r_free; // to use beta prior
  vector<lower=0>[len_free[9]] Psi_sd_free;
  corr_matrix[m] Psi_r_mat[Ng * fullpsi];
  vector<lower=-1,upper=1>[fullpsi ? 0 : len_free[10]] Psi_r_free;
  vector[len_free[13]] Nu_free;
  vector[len_free[14]] Alpha_free;
  vector[len_free[15]] Tau_ufree;

  vector<lower=0,upper=1>[Noent] z_aug; //augmented ordinal data
}
transformed parameters {
  matrix[p, m] Lambda_y[Ng];
  matrix[m, n] Gamma[Ng];
  matrix[m, m] B[Ng];
  matrix[p, p] Theta_sd[Ng];
  matrix[p, p] T_r_lower[Ng];
  matrix[p, p] Theta_r[Ng];
  matrix[p + q, 1] Nu[Ng];
  matrix[m + n, 1] Alpha[Ng];

  matrix[sum(nlevs) - Nord, 1] Tau_un[Ng];
  matrix[sum(nlevs) - Nord, 1] Tau[Ng];
  vector[len_free[15]] Tau_free;
  real tau_jacobian;
  
  matrix[m, m] Psi[Ng];
  
  matrix[m, m] Psi_sd[Ng];
  matrix[m, m] Psi_r_lower[Ng];
  matrix[m, m] Psi_r[Ng];

  vector[len_free[1]] lambda_y_primn;
  vector[len_free[4]] b_primn;
  vector[len_free[13]] nu_primn;
  vector[len_free[14]] alpha_primn;
  vector[len_free[15]] tau_primn;

  matrix[p, m] Lambda_y_A[Ng];     // = Lambda_y * (I - B)^{-1}

  vector[p + q] Mu[Ng];
  matrix[p + q, p + q] Sigma[Ng];  // model covariance matrix
  matrix[p + q, p + q] Sigmainv_grp[Ng];
  real logdetSigma_grp[Ng];
  matrix[p + q + 1, p + q + 1] Sigmainv[Np];  // for updating S^-1 by missing data pattern
  
  matrix[p, q] top_right[Ng]; // top right block of Sigma
  vector[p + q] YXstar[Ntot];
  vector[Nord] YXostar[Ntot]; // ordinal data

  for (g in 1:Ng) {
    // model matrices
    Lambda_y[g] = fill_matrix(Lambda_y_free, Lambda_y_skeleton[g], w1skel, g_start1[g], f_start1[g]);
    Gamma[g] = fill_matrix(Gamma_free, Gamma_skeleton[g], w3skel, g_start3[g], f_start3[g]);
    B[g] = fill_matrix(B_free, B_skeleton[g], w4skel, g_start4[g], f_start4[g]);
    Theta_sd[g] = fill_matrix(Theta_sd_free, Theta_skeleton[g], w5skel, g_start5[g], f_start5[g]);
    T_r_lower[g] = fill_matrix(Theta_r_free, Theta_r_skeleton[g], w7skel, g_start7[g], f_start7[g]);
    Theta_r[g] = T_r_lower[g] + transpose(T_r_lower[g]) - diag_matrix(rep_vector(1, p));
    if (!use_cov) {
      Nu[g] = fill_matrix(Nu_free, Nu_skeleton[g], w13skel, g_start13[g], f_start13[g]);
      Alpha[g] = fill_matrix(Alpha_free, Alpha_skeleton[g], w14skel, g_start14[g], f_start14[g]);
    }

    Psi[g] = diag_matrix(rep_vector(0, m));
  
    if (m > 0) {
      Psi_sd[g] = fill_matrix(Psi_sd_free, Psi_skeleton[g], w9skel, g_start9[g], f_start9[g]);
      if (fullpsi) {
	Psi_r[g] = Psi_r_mat[g];
      } else {
        Psi_r_lower[g] = fill_matrix(Psi_r_free, Psi_r_skeleton[g], w10skel, g_start10[g], f_start10[g]);
        Psi_r[g] = Psi_r_lower[g] + transpose(Psi_r_lower[g]) - diag_matrix(rep_vector(1, m));
      }
      Psi[g] = quad_form_sym(Psi_r[g], Psi_sd[g]);
    }
  }

  // see https://books.google.com/books?id=9AC-s50RjacC&lpg=PP1&dq=LISREL&pg=PA3#v=onepage&q=LISREL&f=false
  for (g in 1:Ng) {
    if (m > 0) {
      Lambda_y_A[g] = mdivide_right(Lambda_y[g], I - B[g]);     // = Lambda_y * (I - B)^{-1}
    }

    if (!use_cov) {
      Mu[g] = to_vector(Nu[g]);
    } else if(has_data) {
      Mu[g] = YXbar[g]; // doesn't enter in likelihood, just for lppd + loo
    }
      
    if (p > 0) {
      Sigma[g, 1:p, 1:p] = quad_form_sym(Theta_r[g], Theta_sd[g]);
      if (m > 0) {
        Sigma[g, 1:p, 1:p] += quad_form_sym(Psi[g], transpose(Lambda_y_A[g]));
	if (!use_cov) Mu[g, 1:p] += to_vector(Lambda_y_A[g] * Alpha[g, 1:m, 1]);
      }
    }
  }

  // obtain ordered thresholds
  if (ord) {
    int opos = 1;
    int ofreepos = 1;
    tau_jacobian = 0;
    for (g in 1:Ng) {
      int vecpos = 1;
      Tau_un[g] = fill_matrix(Tau_ufree, Tau_skeleton[g], w15skel, g_start15[g], f_start15[g]);
      for (i in 1:Nord) {
	for (j in 1:(nlevs[i] - 1)) {
	  real rc = Tau_skeleton[g, vecpos, 1];
	  int eq = w15skel[opos, 1];
	  int wig = w15skel[opos, 3];

	  if (is_inf(rc)) {
	    if (eq == 0 || wig == 1) {
	      if (j == 1) {
		Tau[g, vecpos, 1] = Tau_un[g, vecpos, 1];
	      } else {
		Tau[g, vecpos, 1] = Tau[g, (vecpos - 1), 1] + exp(Tau_un[g, vecpos, 1]);
	      }

	      Tau_free[ofreepos] = Tau[g, vecpos, 1];	      
	      // this is used if a prior goes on Tau_free, instead of Tau_ufree:
	      //if (j > 1) {
	      //  tau_jacobian += Tau_un[g, vecpos, 1]; // see https://mc-stan.org/docs/2_24/reference-manual/ordered-vector.html
	      // }
	      ofreepos += 1;
	    } else if (eq == 1) {
	      int eqent = w15skel[opos, 2];
	      Tau[g, vecpos, 1] = Tau_free[eqent];
	    }	    
	    opos += 1;
	  } else {
	    // fixed value
	    Tau[g, vecpos, 1] = Tau_un[g, vecpos, 1];
	  }	  
	  vecpos += 1;
	}
      }
    }
  }

  // prior vectors
  if (wigind) {
    lambda_y_primn = fill_prior(Lambda_y_free, lambda_y_mn, w1skel);
    b_primn = fill_prior(B_free, b_mn, w4skel);
    nu_primn = fill_prior(Nu_free, nu_mn, w13skel);
    alpha_primn = fill_prior(Alpha_free, alpha_mn, w14skel);
    tau_primn = fill_prior(Tau_ufree, tau_mn, w15skel);
  } else {
    lambda_y_primn = to_vector(lambda_y_mn);
    b_primn = to_vector(b_mn);
    nu_primn = to_vector(nu_mn);
    alpha_primn = to_vector(alpha_mn);
    tau_primn = to_vector(tau_mn);
  }

  // continuous responses underlying ordinal data
  if (ord) {
    int idxvec = 0;
    for (patt in 1:Np) {
      for (i in startrow[patt]:endrow[patt]) {
	for (j in 1:Nordobs[patt]) {
	  int obspos = OrdObsvar[patt,j];
	  int vecpos = YXo[i,obspos] - 1;
	  idxvec += 1;
	  if (obspos > 1) vecpos += sum(nlevs[1:(obspos - 1)]) - (obspos - 1);
	  if (YXo[i,obspos] == 1) {
	    YXostar[i,obspos] = -10 + (Tau[grpnum[patt], (vecpos + 1), 1] + 10) .* z_aug[idxvec];
	    tau_jacobian += log(fabs(Tau[grpnum[patt], (vecpos + 1), 1] + 10));  // must add log(U) to tau_jacobian
	  } else if (YXo[i,obspos] == nlevs[obspos]) {
	    YXostar[i,obspos] = Tau[grpnum[patt], vecpos, 1] + (10 - Tau[grpnum[patt], vecpos, 1]) .* z_aug[idxvec];
	    tau_jacobian += log(fabs(10 - Tau[grpnum[patt], vecpos, 1]));
	  } else {
	    YXostar[i,obspos] = Tau[grpnum[patt], vecpos, 1] + (Tau[grpnum[patt], (vecpos + 1), 1] - Tau[grpnum[patt], vecpos, 1]) .* z_aug[idxvec];
	    tau_jacobian += Tau_un[grpnum[patt], (vecpos + 1), 1]; // jacobian is log(exp(Tau_un))
	  }
	  YXstar[i, ordidx[obspos]] = YXostar[i, obspos];
	}
      }
    }
  }

  if (Ncont > 0) {
    for (patt in 1:Np) {
      for (i in startrow[patt]:endrow[patt]) {
	for (j in 1:Ncont) {
	  YXstar[i, contidx[j]] = YX[i,j];
	}
      }
    }
  }

  // move observations to the left
  if (missing) {
    for (patt in 1:Np) {
      for (i in startrow[patt]:endrow[patt]) {
	for (j in 1:Nobs[patt]) {
	  YXstar[i,j] = YXstar[i, Obsvar[patt,j]];
	}
      }
    }
  }

  // for computing mvn with sufficient stats
  for (g in 1:Ng) {
    Sigmainv_grp[g] = inverse_spd(Sigma[g]);
    logdetSigma_grp[g] = log_determinant(Sigma[g]);
  }
  for (patt in 1:Np) {    
    Sigmainv[patt, 1:(Nobs[patt] + 1), 1:(Nobs[patt] + 1)] = sig_inv_update(Sigmainv_grp[grpnum[patt]], Obsvar[patt,], Nobs[patt], p + q, logdetSigma_grp[grpnum[patt]]);
  }
}
model { // N.B.: things declared in the model block do not get saved in the output, which is okay here

  /* transformed sd parameters for priors */
  vector[len_free[5]] Theta_pri;
  vector[len_free[9]] Psi_pri;
  
  /* log-likelihood */
  if (use_cov && !pri_only) {
    for (g in 1:Ng) {
      target += wishart_lpdf((N[g] - 1) * Sstar[g] | N[g] - 1, Sigma[g]);
      if (Nx[g] > 0) {
	int xvars[Nx[g]] = Xdatvar[g, 1:Nx[g]];
	target += -wishart_lpdf((N[g] - 1) * Sstar[g, xvars, xvars] | N[g] - 1, Sigma[g, xvars, xvars]);
      }
    }
  } else if (has_data && !pri_only) {
    int obsidx[p + q];
    int xidx[p + q];
    int xdatidx[p + q];
    int grpidx;
    int r1;
    int r2;
        
    for (mm in 1:Np) {
      obsidx = Obsvar[mm,];
      xidx = Xvar[mm,];
      xdatidx = Xdatvar[mm,];
      grpidx = grpnum[mm];
      r1 = startrow[mm];
      r2 = endrow[mm];

      if (!use_suff) {
	target += multi_normal_lpdf(YXstar[r1:r2,1:Nobs[mm]] | Mu[grpidx, obsidx[1:Nobs[mm]]], Sigma[grpidx, obsidx[1:Nobs[mm]], obsidx[1:Nobs[mm]]]);

	if (Nx[mm] > 0) {
	  target += -multi_normal_lpdf(YXstar[r1:r2,xdatidx[1:Nx[mm]]] | Mu[grpidx, xidx[1:Nx[mm]]], Sigma[grpidx, xidx[1:Nx[mm]], xidx[1:Nx[mm]]]);
	}
      } else {
	// sufficient stats
	target += multi_normal_suff(YXbarstar[mm, 1:Nobs[mm]], Sstar[mm, 1:Nobs[mm], 1:Nobs[mm]], Mu[grpidx, obsidx[1:Nobs[mm]]], Sigmainv[mm, 1:(Nobs[mm] + 1), 1:(Nobs[mm] + 1)], r2 - r1 + 1);
      
	if (Nx[mm] > 0) {
	  target += -multi_normal_suff(YXbarstar[mm, xdatidx[1:Nx[mm]]], Sstar[mm, xdatidx[1:Nx[mm]], xdatidx[1:Nx[mm]]], Mu[grpidx, xidx[1:Nx[mm]]], sig_inv_update(Sigmainv[grpidx], xidx, Nx[mm], p + q, logdetSigma_grp[grpidx]), r2 - r1 + 1);
	}
      }
    }
    if (ord) {
      target += tau_jacobian;
    }
  }
  
  
  /* prior densities in log-units */
  target += normal_lpdf(Lambda_y_free | lambda_y_primn, lambda_y_sd);
  target += normal_lpdf(Gamma_free    | gamma_mn, gamma_sd);
  target += normal_lpdf(B_free        | b_primn, b_sd);

  target += normal_lpdf(Nu_free       | nu_primn, nu_sd);
  target += normal_lpdf(Alpha_free    | alpha_primn, alpha_sd);
  target += normal_lpdf(Tau_ufree      | tau_primn, tau_sd);

  /* transform sd parameters to var or prec, depending on
     what the user wants. */
  Theta_pri = Theta_sd_free;
  if (len_free[5] > 0 && theta_pow != 1) {
    for (i in 1:len_free[5]) {
      Theta_pri[i] = Theta_sd_free[i]^(theta_pow);
      target += log(fabs(theta_pow)) + (theta_pow - 1)*log(Theta_sd_free[i]);
    }
  }
  Psi_pri = Psi_sd_free;
  if (len_free[9] > 0 && psi_pow != 1) {
    for (i in 1:len_free[9]) {
      Psi_pri[i] = Psi_sd_free[i]^(psi_pow);
      target += log(fabs(psi_pow)) + (psi_pow - 1)*log(Psi_sd_free[i]);
    }
  }

  target += gamma_lpdf(Theta_pri | theta_sd_shape, theta_sd_rate);
  target += gamma_lpdf(Psi_pri | psi_sd_shape, psi_sd_rate);

  target += beta_lpdf(.5 * (1 + Theta_r_free) | theta_r_alpha, theta_r_beta) + log(.5) * len_free[7]; // the latter term is the jacobian moving from (-1,1) to (0,1), because beta_lpdf is defined on (0,1)
  if (fullpsi) {
    for (g in 1:Ng) {
      target += lkj_corr_lpdf(Psi_r_mat[g] | psi_r_alpha[1]);
    }
  } else if (len_free[10] > 0) {
    target += beta_lpdf(.5 * (1 + Psi_r_free) | psi_r_alpha, psi_r_beta) + log(.5) * len_free[10];
  }
}
generated quantities { // these matrices are saved in the output but do not figure into the likelihood
  // see https://books.google.com/books?id=9AC-s50RjacC&lpg=PP1&dq=LISREL&pg=PA34#v=onepage&q=LISREL&f=false

  // sign constraints and correlations
  vector[len_free[1]] ly_sign;
  matrix[p, m] L_Y[Ng];
  //vector[len_free[2]] lx_sign;
  //matrix[q, n] L_X[Ng];
  vector[len_free[3]] g_sign;
  matrix[m, n] Gam[Ng];
  vector[len_free[4]] bet_sign;
  matrix[m, m] Bet[Ng];
  matrix[p, p] Theta[Ng];
  matrix[m, m] PSmat[Ng];
  matrix[m, m] PS[Ng];
  vector[len_free[7]] Theta_cov;
  vector[len_free[5]] Theta_var;
  vector[len_free[10]] P_r;
  vector[len_free[10]] Psi_cov;
  vector[len_free[9]] Psi_var;

  vector[use_cov ? Ng : Ntot] log_lik; // for loo, etc
  vector[use_cov ? Ng : Ntot] log_lik_sat; // for ppp
  vector[p + q] YXstar_rep[Ntot]; // artificial data
  vector[use_cov ? Ng : Ntot] log_lik_rep; // for loo, etc
  vector[use_cov ? Ng : Ntot] log_lik_rep_sat; // for ppp
  matrix[p + q, p + q + 1] satout[Ng];
  matrix[p + q, p + q + 1] satrep_out[Ng];
  vector[p + q] Mu_sat[Ng];
  matrix[p + q, p + q] Sigma_sat[Ng];
  matrix[p + q, p + q] Sigma_sat_inv_grp[Ng];
  real logdetS_sat_grp[Ng];
  matrix[p + q + 1, p + q + 1] Sigma_sat_inv[Np];
  vector[p + q] Mu_rep_sat[Ng];
  matrix[p + q, p + q] Sigma_rep_sat[Ng];
  matrix[p + q, p + q] Sigma_rep_sat_inv_grp[Ng];
  matrix[p + q + 1, p + q + 1] Sigma_rep_sat_inv[Np];
  real logdetS_rep_sat_grp[Ng];
  matrix[p + q, p + q] zmat;
  real<lower=0, upper=1> ppp;
  
  // first deal with sign constraints:
  ly_sign = sign_constrain_load(Lambda_y_free, len_free[1], lam_y_sign);
  bet_sign = sign_constrain_reg(B_free, len_free[4], b_sign, Lambda_y_free, Lambda_y_free);
  if (fullpsi == 0) {
    P_r = sign_constrain_reg(Psi_r_free, len_free[10], psi_r_sign, Lambda_y_free, Lambda_y_free);
  }
  
  for (g in 1:Ng) {
    L_Y[g] = fill_matrix(ly_sign, Lambda_y_skeleton[g], w1skel, g_start1[g], f_start1[g]);

    Gam[g] = fill_matrix(g_sign, Gamma_skeleton[g], w3skel, g_start3[g], f_start3[g]);

    Bet[g] = fill_matrix(bet_sign, B_skeleton[g], w4skel, g_start4[g], f_start4[g]);

    Theta[g] = quad_form_sym(Theta_r[g], Theta_sd[g]);

    if (m > 0) {
      if (fullpsi) {
	PSmat[g] = Psi_r_mat[g];
	PS[g] = quad_form_sym(PSmat[g], Psi_sd[g]);
      } else {
	PSmat[g] = fill_matrix(P_r, Psi_r_skeleton[g], w10skel, g_start10[g], f_start10[g]);
	PS[g] = quad_form_sym(PSmat[g] + transpose(PSmat[g]) - diag_matrix(rep_vector(1, m)), Psi_sd[g]);
      }

    }
    
  }

  // off-diagonal covariance parameter vectors, from cor/sd matrices:
  Theta_cov = cor2cov(Theta_r, Theta_sd, num_elements(Theta_r_free), Theta_r_skeleton, w7skel, Ng);
  Theta_var = Theta_sd_free .* Theta_sd_free;
  if (m > 0 && len_free[10] > 0) {
    /* iden is created so that we can re-use cor2cov, even though
       we don't need to multiply to get covariances */
    matrix[m, m] iden[Ng];
    for (g in 1:Ng) {
      iden[g] = diag_matrix(rep_vector(1, m));
    }
    Psi_cov = cor2cov(PS, iden, len_free[10], Psi_r_skeleton, w10skel, Ng);
  } else {
    Psi_cov = P_r;
  }
  Psi_var = Psi_sd_free .* Psi_sd_free;

  { // log-likelihood
    int obsidx[p + q];
    int xidx[p + q];
    int xdatidx[p + q];
    int r1;
    int r2;
    int grpidx;

    if (do_test && use_cov) {
      for (g in 1:Ng) {
	Sigma_rep_sat[g] = wishart_rng(N[g] - 1, Sigma[g]);
      }
    } else if (do_test && has_data) {
      for (mm in 1:Np) {
	obsidx = Obsvar[mm,];
	xidx = Xvar[mm,];
	xdatidx = Xdatvar[mm,];
	r1 = startrow[mm];
	r2 = endrow[mm];
	grpidx = grpnum[mm];
	for (jj in r1:r2) {
	  YXstar_rep[jj, 1:Nobs[mm]] = multi_normal_rng(Mu[grpidx, obsidx[1:Nobs[mm]]], Sigma[grpidx, obsidx[1:Nobs[mm]], obsidx[1:Nobs[mm]]]);
	}
      }

      if (missing) {
	// start values for Mu and Sigma
	for (g in 1:Ng) {
	  Mu_sat[g] = rep_vector(0, p + q);
	  Mu_rep_sat[g] = Mu_sat[g];
	  Sigma_sat[g] = diag_matrix(rep_vector(1, p + q));
	  Sigma_rep_sat[g] = Sigma_sat[g];
	}

	for (jj in 1:emiter) {
	  satout = estep(YXstar, Mu_sat, Sigma_sat, Nobs, Obsvar, startrow, endrow, grpnum, Np, Ng);
	  satrep_out = estep(YXstar_rep, Mu_rep_sat, Sigma_rep_sat, Nobs, Obsvar, startrow, endrow, grpnum, Np, Ng);

	  // M step
	  for (g in 1:Ng) {
	    Mu_sat[g] = satout[g,,1]/N[g];
	    Sigma_sat[g] = satout[g,,2:(p + q + 1)]/N[g] - Mu_sat[g] * Mu_sat[g]';
	    Mu_rep_sat[g] = satrep_out[g,,1]/N[g];
	    Sigma_rep_sat[g] = satrep_out[g,,2:(p + q + 1)]/N[g] - Mu_rep_sat[g] * Mu_rep_sat[g]';
	  }
	}
      } else {
	// complete data; Np patterns must only correspond to groups
	for (mm in 1:Np) {
	  int arr_dims[3] = dims(YXstar);
	  matrix[endrow[mm] - startrow[mm] + 1, arr_dims[2]] YXsmat; // crossprod needs matrix
	  matrix[endrow[mm] - startrow[mm] + 1, arr_dims[2]] YXsrepmat;
	  r1 = startrow[mm];
	  r2 = endrow[mm];
	  grpidx = grpnum[mm];
	  for (jj in 1:(p + q)) {
	    Mu_sat[grpidx,jj] = mean(YXstar[r1:r2,jj]);
	    Mu_rep_sat[grpidx,jj] = mean(YXstar_rep[r1:r2,jj]);
	  }
	  for (jj in r1:r2) {
	    YXsmat[jj - r1 + 1] = (YXstar[jj] - Mu_sat[grpidx])';
	    YXsrepmat[jj - r1 + 1] = (YXstar_rep[jj] - Mu_rep_sat[grpidx])';
	  }
	  Sigma_sat[grpidx] = crossprod(YXsmat)/N[grpidx];
	  Sigma_rep_sat[grpidx] = crossprod(YXsrepmat)/N[grpidx];
	  // FIXME? Sigma_sat[grpidx] = tcrossprod(YXsmat); does not throw an error??
	}
      }

      for (g in 1:Ng) {
	Sigma_sat_inv_grp[g] = inverse_spd(Sigma_sat[g]);
	logdetS_sat_grp[g] = log_determinant(Sigma_sat[g]);

	Sigma_rep_sat_inv_grp[g] = inverse_spd(Sigma_rep_sat[g]);
	logdetS_rep_sat_grp[g] = log_determinant(Sigma_rep_sat[g]);
      }

      for (mm in 1:Np) {
	Sigma_sat_inv[mm, 1:(Nobs[mm] + 1), 1:(Nobs[mm] + 1)] = sig_inv_update(Sigma_sat_inv_grp[grpnum[mm]], Obsvar[mm,], Nobs[mm], p + q, logdetS_sat_grp[grpnum[mm]]);
	Sigma_rep_sat_inv[mm, 1:(Nobs[mm] + 1), 1:(Nobs[mm] + 1)] = sig_inv_update(Sigma_rep_sat_inv_grp[grpnum[mm]], Obsvar[mm,], Nobs[mm], p + q, logdetS_rep_sat_grp[grpnum[mm]]);
      }
    }

    // compute log-likelihoods
    if (do_test) {
      zmat = rep_matrix(0, p + q, p + q);
      for (mm in 1:Np) {
	obsidx = Obsvar[mm,];
	xidx = Xvar[mm,];
	xdatidx = Xdatvar[mm,];
	r1 = startrow[mm];
	r2 = endrow[mm];
	grpidx = grpnum[mm];

	if (use_cov) {
	  log_lik[mm] = wishart_lpdf((N[mm] - 1) * Sstar[mm] | N[mm] - 1, Sigma[mm]);
	  log_lik_sat[mm] = -log_lik[mm] + wishart_lpdf((N[mm] - 1) * Sstar[mm] | N[mm] - 1, Sstar[mm]);
	  log_lik_rep[mm] = wishart_lpdf(Sigma_rep_sat[mm] | N[mm] - 1, Sigma[mm]);
	  log_lik_rep_sat[mm] = wishart_lpdf(Sigma_rep_sat[mm] | N[mm] - 1, pow(N[mm] - 1, -1) * Sigma_rep_sat[mm]);

	  if (Nx[mm] > 0) {
	    int xvars[Nx[mm]] = xdatidx[1:Nx[mm]];
	    log_lik[mm] += -wishart_lpdf((N[mm] - 1) * Sstar[mm, xvars, xvars] | N[mm] - 1, Sigma[mm, xvars, xvars]);
	    log_lik_sat[mm] += wishart_lpdf((N[mm] - 1) * Sstar[mm, xvars, xvars] | N[mm] - 1, Sigma[mm, xvars, xvars]);
	    log_lik_sat[mm] += -wishart_lpdf((N[mm] - 1) * Sstar[mm, xvars, xvars] | N[mm] - 1, Sstar[mm, xvars, xvars]);
	    log_lik_rep[mm] += -wishart_lpdf(Sigma_rep_sat[mm, xvars, xvars] | N[mm] - 1, Sigma[mm, xvars, xvars]);
	    log_lik_rep_sat[mm] += -wishart_lpdf(Sigma_rep_sat[mm, xvars, xvars] | N[mm] - 1, pow(N[mm] - 1, -1) * Sigma_rep_sat[mm, xvars, xvars]);
	  }
	} else if (has_data) {
	  for (jj in r1:r2) {
	    log_lik[jj] = multi_normal_suff(YXstar[jj, 1:Nobs[mm]], zmat[1:Nobs[mm], 1:Nobs[mm]], Mu[grpidx, obsidx[1:Nobs[mm]]], Sigmainv[mm], 1);
	    // we add loglik[jj] here so that _sat always varies and does not lead to
	    // problems with rhat and neff computations
	    log_lik_sat[jj] = -log_lik[jj] + multi_normal_suff(YXstar[jj, 1:Nobs[mm]], zmat[1:Nobs[mm], 1:Nobs[mm]], Mu_sat[grpidx, obsidx[1:Nobs[mm]]], Sigma_sat_inv[mm], 1);
	
	    log_lik_rep[jj] = multi_normal_suff(YXstar_rep[jj, 1:Nobs[mm]], zmat[1:Nobs[mm], 1:Nobs[mm]], Mu[grpidx, obsidx[1:Nobs[mm]]], Sigmainv[mm], 1);
	    log_lik_rep_sat[jj] = multi_normal_suff(YXstar_rep[jj, 1:Nobs[mm]], zmat[1:Nobs[mm], 1:Nobs[mm]], Mu_rep_sat[grpidx, obsidx[1:Nobs[mm]]], Sigma_rep_sat_inv[mm], 1);

	    // log_lik_sat, log_lik_sat_rep
	    if (Nx[mm] > 0) {
	      log_lik[jj] += -multi_normal_suff(YXstar[jj, xdatidx[1:Nx[mm]]], zmat[1:Nx[mm], 1:Nx[mm]], Mu[grpidx, xidx[1:Nx[mm]]], sig_inv_update(Sigmainv[grpidx], xidx, Nx[mm], p + q, logdetSigma_grp[grpidx]), 1);
	      log_lik_sat[jj] += multi_normal_suff(YXstar[jj, xdatidx[1:Nx[mm]]], zmat[1:Nx[mm], 1:Nx[mm]], Mu[grpidx, xidx[1:Nx[mm]]], sig_inv_update(Sigmainv[grpidx], xidx, Nx[mm], p + q, logdetSigma_grp[grpidx]), 1); // "remove" the fixed x part of log_lik[jj] that we added above
	      log_lik_sat[jj] += -multi_normal_suff(YXstar[jj, xdatidx[1:Nx[mm]]], zmat[1:Nx[mm], 1:Nx[mm]], Mu_sat[grpidx, xidx[1:Nx[mm]]], sig_inv_update(Sigma_sat_inv[grpidx], xidx, Nx[mm], p + q, logdetS_sat_grp[grpidx]), 1);
	  
	      log_lik_rep[jj] += -multi_normal_suff(YXstar_rep[jj, xdatidx[1:Nx[mm]]], zmat[1:Nx[mm], 1:Nx[mm]], Mu[grpidx, xidx[1:Nx[mm]]], sig_inv_update(Sigmainv[grpidx], xidx, Nx[mm], p + q, logdetSigma_grp[grpidx]), 1);

	      log_lik_rep_sat[jj] += -multi_normal_suff(YXstar_rep[jj, xdatidx[1:Nx[mm]]], zmat[1:Nx[mm], 1:Nx[mm]], Mu_rep_sat[grpidx, xidx[1:Nx[mm]]], sig_inv_update(Sigma_rep_sat_inv[grpidx], xidx, Nx[mm], p + q, logdetS_rep_sat_grp[grpidx]), 1);
	    }
	  }
	}
      }
      ppp = step((-sum(log_lik_rep) + sum(log_lik_rep_sat)) - (sum(log_lik_sat)));
    } else {
      ppp = 0;
    }
  }
  
} // end with a completely blank line (not even whitespace)
