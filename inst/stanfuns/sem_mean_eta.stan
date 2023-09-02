  array[] vector sem_mean_eta(array[,,] real alpha, matrix eta, array[,,] real B, array[,,] real gamma, array[] int g, int k, int N, int Ng, int nlv, array[] int lvind, array[] int lv0ind){
    matrix[k,k] iden;
    array[Ng] matrix[k,k] ibinv;
    array[N] vector[k] evlv;
    array[k,1,Ng] real alphvec;
    array[(k - nlv + size(lvind))] int idx;
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
