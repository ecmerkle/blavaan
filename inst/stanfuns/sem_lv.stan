  real sem_lv_lpdf(matrix x, real[,,] alpha, real[,,] B, real[,,] psi, int[] g, int[] regind, int[] exoind, int k, int N, int Ng, int diagpsi, int fullbeta, int nlv, int nexo){
    real ldetcomp[Ng];
    matrix[k,k] iden;
    vector[k] psivecinv[Ng];
    matrix[k,k] psimatinv[Ng];
    matrix[k,k] siginv[Ng];
    vector[k] xvec;
    vector[k] evlv[Ng];
    real xvectm;
    real ldetsum;

    iden = diag_matrix(rep_vector(1.0, k));

    for(j in 1:Ng){
      if(nexo == 0){
        evlv[j] = inverse(iden - to_matrix(B[,,j])) * (to_vector(alpha[,1,j]) + to_matrix(B[,,j]) * to_vector(alpha[,1,j]));
      } else {
        evlv[j,regind] = inverse(iden[regind,regind] - to_matrix(B[regind,regind,j])) * (to_vector(alpha[regind,1,j]) + to_matrix(B[regind,exoind,j]) * to_vector(alpha[exoind,1,j]));
        evlv[j,exoind] = to_vector(alpha[exoind,1,j]);
      }
    }

    if(diagpsi){
      for(j in 1:Ng){
        for(i in 1:k){
          psivecinv[j,i] = 1/psi[i,i,j];
        }
        psimatinv[j] = diag_matrix(psivecinv[j]);

        siginv[j] = (iden - to_matrix(B[,,j])') * psimatinv[j] * (iden - to_matrix(B[,,j]));

	if(fullbeta){
	  ldetcomp[j] = log_determinant(iden - to_matrix(B[,,j]));
	  ldetcomp[j] = ldetcomp[j] + ldetcomp[j] + sum(log(diagonal(to_matrix(psi[,,j]))));
	} else {
          ldetcomp[j] = sum(log(diagonal(to_matrix(psi[,,j]))));
  	}
      }
    } else {
      for(j in 1:Ng){
        psimatinv[j] = to_matrix(psi[,,j]);
	psimatinv[j] = psimatinv[j] + psimatinv[j]' - diag_matrix(diagonal(psimatinv[j]));

	ldetcomp[j] = log_determinant(psimatinv[j]);
	if(fullbeta){
	  ldetcomp[j] = ldetcomp[j] + 2 * log_determinant(iden - to_matrix(B[,,j]));
	}

	psimatinv[j] = inverse_spd(psimatinv[j]);
        siginv[j] = (iden - to_matrix(B[,,j])') * psimatinv[j] * (iden - to_matrix(B[,,j]));
      }
    }

    xvectm = 0;
    ldetsum = 0;
    for(i in 1:N){
      xvec = x[i,]';
      xvectm = xvectm + (xvec - evlv[g[i],])' * siginv[g[i]] * (xvec - evlv[g[i],]);
      ldetsum = ldetsum + ldetcomp[g[i]];
    }

    return -0.5 * (ldetsum + xvectm);
  }
