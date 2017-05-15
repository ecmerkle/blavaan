  real sem_lv_missing_lpdf(matrix x, real[,,] alpha, real[,,] B, real[,,] psi, int[] g, int[] regind, int[] exoind, int k, int N, int Ng, int diagpsi, int fullbeta, int nlv, int nexo, int[] nseen, int[,] obsvar){
    real ldetcomp[Ng];
    matrix[k,k] iden;
    vector[k] psivecinv[Ng];
    matrix[k,k] psimatinv[Ng];
    matrix[k,k] siginv[Ng];
    vector[k] xvec;
    int tmpobs[k];
    real xvectm;
    real ldetsum;

    iden = diag_matrix(rep_vector(1.0, k));

    // TODO speed up by using missingness patterns
    xvectm = 0;
    ldetsum = 0;
    for(i in 1:N){
      if(nlv > 0){
        for(j in 1:nlv){
	  tmpobs[j] = j;
	}
	for(j in (nlv + 1):(nlv + nseen[i])){
	  tmpobs[j] = nlv + obsvar[i,(j - nlv)];
	}
      } else {
        for(j in 1:nseen[i]){
	  tmpobs[j] = obsvar[i,j];
	}
      }

      if(diagpsi){
        for(j in 1:k){
          psivecinv[g[i],j] = 1/psi[j,j,g[i]];
        }
        psimatinv[g[i]] = diag_matrix(psivecinv[g[i]]);

        siginv[g[i],1:(nlv+nseen[i]),1:(nlv+nseen[i])] = (iden[tmpobs[1:(nlv+nseen[i])],tmpobs[1:(nlv+nseen[i])]] - to_matrix(B[tmpobs[1:(nlv+nseen[i])],tmpobs[1:(nlv+nseen[i])],g[i]])') * psimatinv[g[i],tmpobs[1:(nlv+nseen[i])],tmpobs[1:(nlv+nseen[i])]] * (iden[tmpobs[1:(nlv+nseen[i])],tmpobs[1:(nlv+nseen[i])]] - to_matrix(B[tmpobs[1:(nlv+nseen[i])],tmpobs[1:(nlv+nseen[i])],g[i]]));

	if(fullbeta){
	  ldetcomp[g[i]] = log_determinant(iden[tmpobs[1:(nlv+nseen[i])],tmpobs[1:(nlv+nseen[i])]] - to_matrix(B[tmpobs[1:(nlv+nseen[i])],tmpobs[1:(nlv+nseen[i])],g[i]]));
	  ldetcomp[g[i]] = ldetcomp[g[i]] + ldetcomp[g[i]] + sum(log(diagonal(to_matrix(psi[tmpobs[1:(nlv+nseen[i])],tmpobs[1:(nlv+nseen[i])],g[i]]))));
	} else {
          ldetcomp[g[i]] = sum(log(diagonal(to_matrix(psi[tmpobs[1:(nlv+nseen[i])],tmpobs[1:(nlv+nseen[i])],g[i]]))));
  	}
      } else {
        psimatinv[g[i]] = to_matrix(psi[,,g[i]]);
	psimatinv[g[i]] = psimatinv[g[i]] + psimatinv[g[i]]' - diag_matrix(diagonal(psimatinv[g[i]]));

	ldetcomp[g[i]] = log_determinant(psimatinv[g[i],tmpobs[1:(nlv+nseen[i])],tmpobs[1:(nlv+nseen[i])]]);
	if(fullbeta){
	  ldetcomp[g[i]] = ldetcomp[g[i]] + 2 * log_determinant(iden[tmpobs[1:(nlv+nseen[i])],tmpobs[1:(nlv+nseen[i])]] - to_matrix(B[tmpobs[1:(nlv+nseen[i])],tmpobs[1:(nlv+nseen[i])],g[i]]));
	}

	psimatinv[g[i]] = inverse_spd(psimatinv[g[i]]);
        siginv[g[i],1:(nlv+nseen[i]),1:(nlv+nseen[i])] = (iden[tmpobs[1:(nlv+nseen[i])],tmpobs[1:(nlv+nseen[i])]] - to_matrix(B[tmpobs[1:(nlv+nseen[i])],tmpobs[1:(nlv+nseen[i])],g[i]])') * inverse_spd(psimatinv[g[i],tmpobs[1:(nlv+nseen[i])],tmpobs[1:(nlv+nseen[i])]]) * (iden[tmpobs[1:(nlv+nseen[i])],tmpobs[1:(nlv+nseen[i])]] - to_matrix(B[tmpobs[1:(nlv+nseen[i])],tmpobs[1:(nlv+nseen[i])],g[i]]));
      }
    
      xvec[1:(nlv+nseen[i])] = x[i,1:(nlv+nseen[i])]';
      xvectm = xvectm + (xvec[1:(nlv+nseen[i])] - to_vector(alpha[tmpobs[1:(nlv+nseen[i])],1,g[i]]))' * siginv[g[i],1:(nlv+nseen[i]),1:(nlv+nseen[i])] * (xvec[1:(nlv+nseen[i])] - to_vector(alpha[tmpobs[1:(nlv+nseen[i])],1,g[i]]));
      ldetsum = ldetsum + ldetcomp[g[i]];
    }

    return -0.5 * (ldetsum + xvectm);
  }
