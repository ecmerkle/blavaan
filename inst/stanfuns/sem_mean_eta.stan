  vector[] sem_mean_eta(real[,,] alpha, matrix eta, real[,,] B, real[,,] gamma, int[] g, int k, int N, int Ng, int nlv, int[] lvind, int[] lv0ind){
    matrix[k,k] iden;
    matrix[k,k] ibinv[Ng];
    vector[k] evlv[N];
    real alphvec[k,1,Ng];
    int idx[(k - nlv + size(lvind))];
    int nov;
    int nlvno0;

    nov = k - nlv;
    nlvno0 = size(lvind);

    iden = diag_matrix(rep_vector(1.0, k));

    alphvec = alpha;

    if(size(lvind) > 0){
      idx[1:nlvno0] = lvind;
    }
    if(nov > 0){
      for(j in 1:nov){
        idx[nlvno0+j] = nlv + j; //nlvno0 + j?
      }
    }

    for(j in 1:Ng){
      ibinv[j,lv0ind,lv0ind] = inverse(iden[lv0ind,lv0ind] - to_matrix(B[lv0ind,lv0ind,j]));
    }

    for(i in 1:N){
      // this line took way too long to get right:
      evlv[i,lv0ind] = ibinv[g[i],lv0ind,lv0ind] * (to_vector(alphvec[lv0ind,1,g[i]]) + to_matrix(B[lv0ind,idx,g[i]]) * eta[i,idx]');
    }

    return evlv;
  }
