  matrix fill_upper(matrix x){
    for(i in 1:(rows(x) - 1)){
      for(j in (i+1):rows(x)){
        x[i,j] = x[j,i];
      }
    }
    return x;
  }