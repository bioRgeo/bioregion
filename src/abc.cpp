// Inspired by bdist.c from the R package ecodist (version 2.0.7)

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix abc(NumericMatrix comat) {
  int nrow = comat.nrow();
  int ncol = comat.ncol();
  NumericMatrix res(nrow*(nrow-1)/2,5);
  int l=0;
  for(int i = 0; i < (nrow - 1); i++) {
    for(int j = (i + 1); j < (nrow); j++) {
      double minsum = 0;
      double sumi = 0;
      double sumj = 0;
      for(int   k = 0; k < ncol; k++) {
        if(comat(i,k) < comat(j,k))
          minsum += comat(i,k);
        else
          minsum += comat(j,k);
        sumi += comat(i,k);
        sumj += comat(j,k);
      }
      res(l,0)=(i+1);
      res(l,1)=(j+1);
      res(l,2)=minsum;
      res(l,3)=sumi-minsum;
      res(l,4)=sumj-minsum;
      l++;
    }
  }

  return res;
}
