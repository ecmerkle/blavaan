  real sem_lv_missing_lpdf(matrix x, real[,,] alpha, real[,,] B, real[,,] psi, int[] g, int[] regind, int[] exoind, int k, int N, int Ng, int diagpsi, int fullbeta, int nlv, int nexo, int[] lvind, int nlvno0, int[] nseen, int[,] obsvar, int[] nseenexo, int[,] obsexo){
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

    evlv = sem_mean(alpha, B, g, regind, exoind, k, Ng, nexo);

    // TODO speed up by using missingness patterns
    xvectm = 0;
    ldetsum = 0;
    for(i in 1:N){
      if(nlvno0 > 0){
        for(j in 1:nlvno0){
	  tmpobs[j] = lvind[j];
	}
	for(j in (nlvno0 + 1):(nlvno0 + nseen[i])){
	  tmpobs[j] = nlv + obsvar[i,(j - nlv)];
	}
      } else {
        for(j in 1:nseen[i]){
	  tmpobs[j] = nlv + obsvar[i,j];
	}
      }

      // remove after testing this is needed for reg model with 
      // 0 constraint, but not coded as exo?
      //evlv[1,1] = alpha[1,1,1] + to_vector(B[1,2:3,1])' * to_vector(alpha[2:3,1,1]);
      //if(nseen[i] > 1){
      //  evlv[g[i],2:nseen[i]] = to_vector(alpha[tmpobs[2:nseen[i]],1,g[i]]);
      //}

      if(diagpsi){
        for(j in 1:(nlvno0 + nseen[i])){
          psivecinv[g[i],tmpobs[j]] = 1/psi[tmpobs[j],tmpobs[j],g[i]];
        }
        psimatinv[g[i]] = diag_matrix(psivecinv[g[i]]);

        siginv[g[i],1:(nlvno0+nseen[i]),1:(nlvno0+nseen[i])] = (iden[tmpobs[1:(nlvno0+nseen[i])],tmpobs[1:(nlvno0+nseen[i])]] - to_matrix(B[tmpobs[1:(nlvno0+nseen[i])],tmpobs[1:(nlvno0+nseen[i])],g[i]])') * psimatinv[g[i],tmpobs[1:(nlvno0+nseen[i])],tmpobs[1:(nlvno0+nseen[i])]] * (iden[tmpobs[1:(nlvno0+nseen[i])],tmpobs[1:(nlvno0+nseen[i])]] - to_matrix(B[tmpobs[1:(nlvno0+nseen[i])],tmpobs[1:(nlvno0+nseen[i])],g[i]]));

	if(fullbeta){
	  ldetcomp[g[i]] = log_determinant(iden[tmpobs[1:(nlvno0+nseen[i])],tmpobs[1:(nlvno0+nseen[i])]] - to_matrix(B[tmpobs[1:(nlvno0+nseen[i])],tmpobs[1:(nlvno0+nseen[i])],g[i]]));
	  ldetcomp[g[i]] = ldetcomp[g[i]] + ldetcomp[g[i]] + sum(log(diagonal(to_matrix(psi[tmpobs[1:(nlvno0+nseen[i])],tmpobs[1:(nlvno0+nseen[i])],g[i]]))));
	} else {
          ldetcomp[g[i]] = sum(log(diagonal(to_matrix(psi[tmpobs[1:(nlvno0+nseen[i])],tmpobs[1:(nlvno0+nseen[i])],g[i]]))));
  	}
      } else {
        psimatinv[g[i]] = to_matrix(psi[,,g[i]]);
	psimatinv[g[i]] = psimatinv[g[i]] + psimatinv[g[i]]' - diag_matrix(diagonal(psimatinv[g[i]]));

	ldetcomp[g[i]] = log_determinant(psimatinv[g[i],tmpobs[1:(nlvno0+nseen[i])],tmpobs[1:(nlvno0+nseen[i])]]);
	if(fullbeta){
	  ldetcomp[g[i]] = ldetcomp[g[i]] + 2 * log_determinant(iden[tmpobs[1:(nlvno0+nseen[i])],tmpobs[1:(nlvno0+nseen[i])]] - to_matrix(B[tmpobs[1:(nlvno0+nseen[i])],tmpobs[1:(nlvno0+nseen[i])],g[i]]));
	}

	psimatinv[g[i]] = inverse_spd(psimatinv[g[i]]);
        siginv[g[i],1:(nlvno0+nseen[i]),1:(nlvno0+nseen[i])] = (iden[tmpobs[1:(nlvno0+nseen[i])],tmpobs[1:(nlvno0+nseen[i])]] - to_matrix(B[tmpobs[1:(nlvno0+nseen[i])],tmpobs[1:(nlvno0+nseen[i])],g[i]])') * psimatinv[g[i],tmpobs[1:(nlvno0+nseen[i])],tmpobs[1:(nlvno0+nseen[i])]] * (iden[tmpobs[1:(nlvno0+nseen[i])],tmpobs[1:(nlvno0+nseen[i])]] - to_matrix(B[tmpobs[1:(nlvno0+nseen[i])],tmpobs[1:(nlvno0+nseen[i])],g[i]]));
      }

      xvec[1:(nlvno0+nseen[i])] = x[i,1:(nlvno0+nseen[i])]';
      xvectm = xvectm + (xvec[1:(nlvno0+nseen[i])] - evlv[g[i],tmpobs[1:(nlvno0+nseen[i])]])' * siginv[g[i],1:(nlvno0+nseen[i]),1:(nlvno0+nseen[i])] * (xvec[1:(nlvno0+nseen[i])] - evlv[g[i],tmpobs[1:(nlvno0+nseen[i])]]);
      ldetsum = ldetsum + ldetcomp[g[i]];
    }

    return -0.5 * (ldetsum + xvectm);
  }
