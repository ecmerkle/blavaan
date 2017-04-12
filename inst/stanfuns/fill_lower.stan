  matrix fill_lower(matrix x){
    matrix[rows(x),cols(x)] newx;

    newx = x;
    for(i in 1:(rows(x) - 1)){
      for(j in (i+1):rows(x)){
        newx[j,i] = x[i,j];
      }
    }
    return newx;
  }
