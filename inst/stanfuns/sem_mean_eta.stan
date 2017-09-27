  vector[] sem_mean_eta(real[,,] alpha, matrix eta, real[,,] B, int[] g, int[] regind, int[] exoind, int k, int N, int nexo, int Ng, int[] lvind){
    matrix[k,k] iden;
    matrix[k,k] ibinv[Ng];
    vector[k] evlv[N];
    real alphvec[k,1,1];

    iden = diag_matrix(rep_vector(1.0, k));

    alphvec = alpha;
    alphvec[lvind,1,1] = rep_array(0,size(lvind));

    for(j in 1:Ng){
      if(nexo == 0){
        ibinv[j] = inverse(iden - to_matrix(B[,,j]));
      } else {
        ibinv[j,1:size(regind),1:size(regind)] = inverse(iden[regind,regind] - to_matrix(B[regind,regind,j]));
      }
    }

    if(nexo == 0){
      for(i in 1:N){
        evlv[i] = ibinv[g[i]] * (to_vector(alphvec[,1,g[i]]) + to_matrix(B[,,g[i]]) * eta[i,]');
      }
    } else {
      for(i in 1:N){
        evlv[i,regind] = ibinv[g[i],1:size(regind),1:size(regind)] * (to_vector(alphvec[regind,1,g[i]]) + to_matrix(B[regind,exoind,g[i]]) * eta[i,]');
        evlv[i,exoind] = to_vector(alphvec[exoind,1,g[i]]);
      }
    }

    return evlv;
  }
