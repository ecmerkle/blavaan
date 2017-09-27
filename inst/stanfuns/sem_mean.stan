  vector[] sem_mean(real[,,] alpha, real[,,] B, int[] g, int[] regind, int[] exoind, int k, int Ng, int nexo){
    matrix[k,k] iden;
    vector[k] evlv[Ng];

    iden = diag_matrix(rep_vector(1.0, k));

    for(j in 1:Ng){
      if(nexo == 0){
        evlv[j] = inverse(iden - to_matrix(B[,,j])) * (to_vector(alpha[,1,j]) + to_matrix(B[,,j]) * to_vector(alpha[,1,j]));
      } else {
        evlv[j,regind] = inverse(iden[regind,regind] - to_matrix(B[regind,regind,j])) * (to_vector(alpha[regind,1,j]) + to_matrix(B[regind,exoind,j]) * to_vector(alpha[exoind,1,j]));
        evlv[j,exoind] = to_vector(alpha[exoind,1,j]);
      }
    }

    return evlv;
  }
