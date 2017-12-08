  real sem_lv_lpdf(matrix x, real[,,] alpha, real[,,] B, real[,,] psi, real[,,] gamma, int gamind, real[,] meanx, int[] g, int k, int N, int Ng, int diagpsi, int fullbeta, int nlv, int[] lvind, int nlvno0){
    real ldetcomp[Ng];
    matrix[k,k] iden;
    vector[k] alpha2[Ng];
    vector[k] psivecinv[Ng];
    matrix[k,k] psimatinv[Ng];
    matrix[k,k] psimat[Ng];
    matrix[k,k] siginv[Ng];
    vector[k] xvec;
    vector[k] evlv[Ng];
    int idx[(k-nlv+nlvno0)];
    real xvectm;
    real ldetsum;
    int nov;
    int nidx;

    nov = k - nlv;
    nidx = nov + nlvno0;

    iden = diag_matrix(rep_vector(1.0, k));

    if(nlvno0 > 0){
      idx[1:nlvno0] = lvind;
    }
    if(nov > 0){
      for(j in 1:nov){
        idx[nlvno0+j] = nlv + j; //nlvno0 + j?
      }
    }

    for(j in 1:Ng){
      alpha2[j] = to_vector(alpha[,1,j]);
    }

    evlv = sem_mean(alpha2, B, gamma, g, k, Ng, gamind, meanx);

    if(diagpsi){
      for(j in 1:Ng){
        for(i in 1:nidx){
          psivecinv[j,idx[i]] = 1/psi[idx[i],idx[i],j];
        }
        psimatinv[j] = diag_matrix(psivecinv[j]);

        siginv[j,1:nidx,1:nidx] = (iden[idx,idx] - to_matrix(B[idx,idx,j])') * psimatinv[j,idx,idx] * (iden[idx,idx] - to_matrix(B[idx,idx,j]));

	if(fullbeta){
	  ldetcomp[j] = log_determinant(iden[idx,idx] - to_matrix(B[idx,idx,j]));
	  ldetcomp[j] = -2 * ldetcomp[j] + sum(log(diagonal(to_matrix(psi[idx,idx,j]))));
	} else {
          ldetcomp[j] = sum(log(diagonal(to_matrix(psi[idx,idx,j]))));
  	}
      }
    } else {
      for(j in 1:Ng){
	psimat[j] = to_matrix(psi[,,j]) + to_matrix(psi[,,j])' - diag_matrix(diagonal(to_matrix(psi[,,j])));

	ldetcomp[j] = log_determinant(psimat[j,idx,idx]);
	if(fullbeta){
	  ldetcomp[j] = ldetcomp[j] - 2 * log_determinant(iden[idx,idx] - to_matrix(B[idx,idx,j]));
	}

	psimatinv[j] = psimat[j];
	psimatinv[j,1:nidx,1:nidx] = inverse_spd(psimat[j,idx,idx]);
        siginv[j,1:nidx,1:nidx] = (iden[idx,idx] - to_matrix(B[idx,idx,j])') * psimatinv[j,1:nidx,1:nidx] * (iden[idx,idx] - to_matrix(B[idx,idx,j]));
      }
    }

    xvectm = 0;
    ldetsum = 0;
    for(i in 1:N){
      xvec = x[i,]';
      xvectm = xvectm + (xvec[idx] - evlv[g[i],idx])' * siginv[g[i],1:nidx,1:nidx] * (xvec[idx] - evlv[g[i],idx]);
      ldetsum = ldetsum + ldetcomp[g[i]];
    }

    return -0.5 * (ldetsum + xvectm);
  }
