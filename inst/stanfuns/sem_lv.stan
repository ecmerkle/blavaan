  real sem_lv_lpdf(matrix x, real[,,] alpha, real[,,] B, real[,,] psi, int[] g, int k, int N, int Ng, int diagpsi){
    real ldetcomp[Ng];
    matrix[k,k] iden;
    vector[k] psivecinv[Ng];
    matrix[k,k] psimatinv[Ng];
    matrix[k,k] siginv[Ng];
    vector[k] xvec;
    real xvectm;
    real ldetsum;

    iden = diag_matrix(rep_vector(1.0, k));

    if(diagpsi){
      for(j in 1:Ng){
        for(i in 1:k){
          psivecinv[j,i] = 1/psi[i,i,j];
        }
        psimatinv[j] = diag_matrix(psivecinv[j]);

        siginv[j] = (iden - to_matrix(B[,,j])') * psimatinv[j] * (iden - to_matrix(B[,,j]));
        ldetcomp[j] = sum(log(diagonal(to_matrix(psi[,,j]))));
      }
    } else {
      for(j in 1:Ng){
        psimatinv[j] = to_matrix(psi[,,j]);
	psimatinv[j] = psimatinv[j] + psimatinv[j]' - diag_matrix(diagonal(psimatinv[j]));
	ldetcomp[j] = log_determinant(psimatinv[j]);
	psimatinv[j] = inverse_spd(psimatinv[j]);
        siginv[j] = (iden - to_matrix(B[,,j])') * psimatinv[j] * (iden - to_matrix(B[,,j]));
      }
    }

    xvectm = 0;
    ldetsum = 0;
    for(i in 1:N){
      xvec = x[i,]';
      xvectm = xvectm + (xvec - to_vector(alpha[,1,g[i]]))' * siginv[g[i]] * (xvec - to_vector(alpha[,1,g[i]]));
      ldetsum = ldetsum + ldetcomp[g[i]];
    }

    return -0.5 * (ldetsum + xvectm);
  }
