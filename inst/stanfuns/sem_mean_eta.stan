  vector[] sem_mean_eta(real[,,] alpha, matrix eta, real[,,] B, real[,,] gamma, int[] g, int k, int N, int Ng, int[] lvind){
    matrix[k,k] iden;
    matrix[k,k] ibinv[Ng];
    vector[k] evlv[N];
    real alphvec[k,1,Ng];

    iden = diag_matrix(rep_vector(1.0, k));

    alphvec = alpha;

    for(j in 1:Ng){
      ibinv[j] = inverse(iden - to_matrix(B[,,j]));
      alphvec[lvind,1,j] = rep_array(0,size(lvind));
    }

    for(i in 1:N){
      evlv[i] = ibinv[g[i]] * (to_vector(alphvec[,1,g[i]]) + to_matrix(B[,,g[i]]) * eta[i,]');
    }

    return evlv;
  }
