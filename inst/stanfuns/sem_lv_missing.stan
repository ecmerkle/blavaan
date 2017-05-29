  real sem_lv_missing_lpdf(matrix x, real[,,] alpha, real[,,] B, real[,,] psi, int[] g, int[] regind, int[] exoind, int k, int N, int Ng, int diagpsi, int fullbeta, int nlv, int nexo, int[] nseen, int[,] obsvar, int[] nseenexo, int[,] obsexo){
    real ldetcomp[Ng];
    matrix[k,k] iden;
    vector[k] psivecinv[Ng];
    matrix[k,k] psimatinv[Ng];
    matrix[k,k] siginv[Ng];
    vector[k] xvec;
    vector[k] evlv[Ng];
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

      if(nexo == 0){
        evlv[g[i]] = inverse(iden - to_matrix(B[,,g[i]])) * (to_vector(alpha[,1,g[i]]) + to_matrix(B[,,g[i]]) * to_vector(alpha[,1,g[i]]));
        evlv[g[i],1:(nlv+nseen[i])] = evlv[g[i],tmpobs[1:(nlv+nseen[i])]];
      } else {
        if(nseen[i] - nseenexo[i] > 0){
          evlv[g[i],1:(nlv+nseen[i]-nseenexo[i])] = inverse(iden[tmpobs[1:(nlv+nseen[i]-nseenexo[i])],tmpobs[1:(nlv+nseen[i]-nseenexo[i])]] - to_matrix(B[tmpobs[1:(nlv+nseen[i]-nseenexo[i])],tmpobs[1:(nlv+nseen[i]-nseenexo[i])],g[i]])) * (to_vector(alpha[tmpobs[1:(nlv+nseen[i]-nseenexo[i])],1,g[i]]) + to_matrix(B[tmpobs[1:(nlv+nseen[i]-nseenexo[i])],tmpobs[1:(nlv+nseen[i]-nseenexo[i])],g[i]]) * to_vector(alpha[tmpobs[1:(nlv+nseen[i]-nseenexo[i])],1,g[i]]));
	} else if(nlv > 0){
          evlv[g[i],1:nlv] = inverse(iden[tmpobs[1:nlv],tmpobs[1:nlv]] - to_matrix(B[tmpobs[1:nlv],tmpobs[1:nlv],g[i]])) * (to_vector(alpha[tmpobs[1:nlv],1,g[i]]) + to_matrix(B[tmpobs[1:nlv],tmpobs[1:nlv],g[i]]) * to_vector(alpha[tmpobs[1:nlv],1,g[i]]));
	}

        if(nseenexo[i] > 0){
          // assume exo is always at the end
          evlv[g[i],(nlv + nseen[i] - nseenexo[i] + 1):(nlv + nseen[i])] = to_vector(alpha[obsexo[i,1:nseenexo[i]],1,g[i]]);
        }
      }

      // remove after testing this is needed for reg model with 
      // 0 constraint, but not coded as exo?
      evlv[1,1] = alpha[1,1,1] + to_vector(B[1,2:3,1])' * to_vector(alpha[2:3,1,1]);
      if(nseen[i] > 1){
        evlv[g[i],2:nseen[i]] = to_vector(alpha[tmpobs[2:nseen[i]],1,g[i]]);
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
        siginv[g[i],1:(nlv+nseen[i]),1:(nlv+nseen[i])] = (iden[tmpobs[1:(nlv+nseen[i])],tmpobs[1:(nlv+nseen[i])]] - to_matrix(B[tmpobs[1:(nlv+nseen[i])],tmpobs[1:(nlv+nseen[i])],g[i]])') * psimatinv[g[i],tmpobs[1:(nlv+nseen[i])],tmpobs[1:(nlv+nseen[i])]] * (iden[tmpobs[1:(nlv+nseen[i])],tmpobs[1:(nlv+nseen[i])]] - to_matrix(B[tmpobs[1:(nlv+nseen[i])],tmpobs[1:(nlv+nseen[i])],g[i]]));
      }

      xvec[1:(nlv+nseen[i])] = x[i,1:(nlv+nseen[i])]';
      xvectm = xvectm + (xvec[1:(nlv+nseen[i])] - evlv[g[i],1:(nlv+nseen[i])])' * siginv[g[i],1:(nlv+nseen[i]),1:(nlv+nseen[i])] * (xvec[1:(nlv+nseen[i])] - evlv[g[i],1:(nlv+nseen[i])]);
      ldetsum = ldetsum + ldetcomp[g[i]];
    }

    return -0.5 * (ldetsum + xvectm);
  }
