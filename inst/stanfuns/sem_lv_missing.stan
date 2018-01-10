  real sem_lv_missing_lpdf(matrix x, real[,,] alpha, real[,,] B, real[,,] psi, real[,,] gamma, int gamind, real[,] meanx, int[] g, int k, int N, int Ng, int diagpsi, int fullbeta, int nlv, int[] lvind, int nlvno0, int[,] nseen, int[,,] obsvar, int[] obspatt, int[] gpatt){
    real ldetcomp[Ng,max(gpatt)];
    matrix[k,k] iden;
    vector[k] alpha2[Ng];
    vector[k] psivecinv[Ng];
    matrix[k,k] psimatinv[Ng];
    matrix[k,k] psimat[Ng];
    matrix[k,k] siginv[Ng,max(gpatt)];
    vector[k] xvec;
    vector[k] evlv[Ng];
    int idx[(k-nlv+nlvno0)];
    int tmpobs[k];
    real xvectm;
    real ldetsum;
    int nov;
    int nidx;

    nov = k - nlv;

    iden = diag_matrix(rep_vector(1.0, k));

    for(j in 1:Ng){
      alpha2[j] = to_vector(alpha[,1,j]);
    }

    evlv = sem_mean(alpha2, B, gamma, g, k, Ng, gamind, meanx);

    //     compute siginv, ldetcomp by missingness pattern
    //     siginv: matrix[k,k] siginv[Ng,max(gpatt)]
    //     ldetcomp: vector[max(gpatt)] ldetcomp[Ng]
    for(gg in 1:Ng){
      for(m in 1:gpatt[gg]){
        if(nlvno0 > 0){
	  idx[1:nlvno0] = lvind;
	}
	if(nov > 0){
	  for(j in 1:nseen[gg,m]){
	    idx[nlvno0+j] = nlv + obsvar[gg,m,j]; //nlv + obsvar[i,(j - nlv)];
	  }
        }
	nidx = nlvno0 + nseen[gg,m];

        if(diagpsi){
          for(j in 1:nidx){
            psivecinv[gg,idx[j]] = 1/psi[idx[j],idx[j],gg];
          }
          psimatinv[gg] = diag_matrix(psivecinv[gg]);

          siginv[gg,m,1:nidx,1:nidx] = (iden[idx[1:nidx],idx[1:nidx]] - to_matrix(B[idx[1:nidx],idx[1:nidx],gg])') * psimatinv[gg,idx[1:nidx],idx[1:nidx]] * (iden[idx[1:nidx],idx[1:nidx]] - to_matrix(B[idx[1:nidx],idx[1:nidx],gg]));

	  if(fullbeta){
	    ldetcomp[gg,m] = log_determinant(iden[idx[1:nidx],idx[1:nidx]] - to_matrix(B[idx[1:nidx],idx[1:nidx],gg]));
	    ldetcomp[gg,m] = -2 * ldetcomp[gg,m] + sum(log(diagonal(to_matrix(psi[idx[1:nidx],idx[1:nidx],gg]))));
	  } else {
            ldetcomp[gg,m] = sum(log(diagonal(to_matrix(psi[idx[1:nidx],idx[1:nidx],gg]))));
  	  }
        } else {
          psimat[gg] = to_matrix(psi[,,gg]) + to_matrix(psi[,,gg])' - diag_matrix(diagonal(to_matrix(psi[,,gg])));

	  ldetcomp[gg,m] = log_determinant(psimat[gg,idx[1:nidx],idx[1:nidx]]);
	  if(fullbeta){
	    ldetcomp[gg,m] = ldetcomp[gg,m] - 2 * log_determinant(iden[idx[1:nidx],idx[1:nidx]] - to_matrix(B[idx[1:nidx],idx[1:nidx],gg]));
	  }

	  psimatinv[gg,1:nidx,1:nidx] = inverse_spd(psimat[gg,idx[1:nidx],idx[1:nidx]]);
          siginv[gg,m,1:nidx,1:nidx] = (iden[idx[1:nidx],idx[1:nidx]] - to_matrix(B[idx[1:nidx],idx[1:nidx],gg])') * psimatinv[gg,1:nidx,1:nidx] * (iden[idx[1:nidx],idx[1:nidx]] - to_matrix(B[idx[1:nidx],idx[1:nidx],gg]));
        }
      }
    }

    // now that ldetcomp and siginv computed for each pattern,
    // obtain log-likelihood
    xvectm = 0;
    ldetsum = 0;
    for(i in 1:N){
      if(nlvno0 > 0){
        idx[1:nlvno0] = lvind;
      }
      if(nov > 0){
        for(j in 1:nseen[g[i],obspatt[i]]){
	  idx[nlvno0+j] = nlv + obsvar[g[i],obspatt[i],j]; //nlv + obsvar[i,(j - nlv)];
	}
      }
      nidx = nlvno0 + nseen[g[i],obspatt[i]];

      xvec[1:nidx] = x[i,1:nidx]';
      xvectm = xvectm + (xvec[1:nidx] - evlv[g[i],idx[1:nidx]])' * siginv[g[i],obspatt[i],1:nidx,1:nidx] * (xvec[1:nidx] - evlv[g[i],idx[1:nidx]]);
      ldetsum = ldetsum + ldetcomp[g[i],obspatt[i]];
    }

    return -0.5 * (ldetsum + xvectm);
  }
