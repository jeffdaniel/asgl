#include <stdio.h>
#include <math.h>
#include <R.h>

// Euclidean norm of x
double l2_norm(double *x, int p) {
  double x_norm = 0;
  for (int i = 0; i < p; i++) {
    x_norm = x_norm + pow(x[i], 2);
  }
  x_norm = sqrt(x_norm);
  return(x_norm);
}

// Cross product of the jth column of x with y
double crossprod(double *x, double *y, int n, int j) {
  double val = 0;
  int nn = n * j;
  for (int i = 0; i < n; i++) {
    val += x[nn + i] * y[i];
  }
  return(val);
}

// Soft-thresholding operator
double soft_threshold(double z, double l) {
  if (z > l)  return(z-l);
  if (z < -l) return(z+l);
  return(0);
}

// -----------------------------------------------------------------------------

void linGradCalc(int *nrow, double *eta, double *y, double *ldot) {
  for (int i = 0; i < nrow[0]; i++) {
    ldot[i] = (eta[i] - y[i])/nrow[0];
  }
}

double linNegLogLikelihoodCalc(int *nrow, double *eta, double *y) {
  double squareSum = 0;
  for (int i = 0; i < nrow[0]; i++) {
    squareSum = squareSum + pow(eta[i] - y[i], 2)/2;
  }
  return squareSum/nrow[0];
}

void linSolver(double *X, double *y, int *nrow, int *ncol, int *numGroup,
               double *beta, int *rangeGroupInd, int *groupLen,
               double *lambda1, double *lambda2, int *innerIter,
               double *thresh, double *ldot, double *gamma, double *eta,
               int* betaIsZero, int* isActive, double *step,
               double *grpWeights, double *indWeights) {

  int n = nrow[0];
  int p = ncol[0];
  int count = 0;

  double* theta = Calloc(p, double);
  double* etaNew  = Calloc(n, double);
  double* etaNull = Calloc(n, double);
  double zeroCheck = 0;
  double check = 0;
  double t = step[0];
  double diff = 1;
  double norm = 0;
  double uOp = 0;
  double Lnew, Lold = 0;
  double sqNormG = 0;
  double iProd = 0;

  // OUTER LOOP
  for (int i = 0; i < numGroup[0]; i++) {

    int g = groupLen[i];
    int groupStart = rangeGroupInd[i];
    int groupEnd = rangeGroupInd[i] + groupLen[i];
    double* grad = Calloc(g, double);

    if (isActive[i] == 0) { // Group is inactive

      // Compute gradient
      for (int j = groupStart; j < groupEnd; j++) {
        for (int k = 0; k < n; k++) {
          etaNull[k] = eta[k] - X[k + n * j] * beta[j];
        }
      }
      linGradCalc(nrow, etaNull, y, ldot); // Compute log-likelihood gradient
      for (int j = 0; j < g; j++) {
        grad[j] = soft_threshold(crossprod(X, ldot, n, groupStart + j),
                                 lambda1[0] * indWeights[groupStart + j]);
      }

      // Check if group coefficients are identically zero
      zeroCheck = l2_norm(grad, groupLen[i]);
      if (zeroCheck <= pow(lambda2[0] * grpWeights[i], 2) * g) {
        // Inequality satisfied --> all coefficients set to zero
        if (betaIsZero[i] == 0) {
          for (int k = 0; k < n; k++) {
	          eta[k] = etaNull[k];
	        }
        }
        betaIsZero[i] = 1;
        for(int j = 0; j < groupLen[i]; j++) {
          beta[j + rangeGroupInd[i]] = 0;
        }
      } else {
        // Otherwise --> proceed to inner loop
        isActive[i] = 1;
      }
    }

    if (isActive[i] == 1) { // Compute group coefficient estimates

      for (int j = 0; j < p; j++) {
        theta[j] = beta[j];
      }

  	  betaIsZero[i] = 0;
      double* z = Calloc(g, double);
      double* U = Calloc(g, double);
      double* G = Calloc(g, double);
      double* betaNew = Calloc(g, double);

  	  do { // Iterate until convergence

        // Compute current log-likelihood
        linGradCalc(nrow, eta, y ,ldot);
        for (int j = 0; j < g; j++) {
          grad[j] = crossprod(X, ldot, n, groupStart + j);
        }
        Lold = linNegLogLikelihoodCalc(nrow, eta, y);

        do { // Perform back-tracking to determine optimal descent step size

          for (int j = 0; j < g; j++) {
	          z[j] = soft_threshold(beta[groupStart + j] - t * grad[j],
                                  lambda1[0] * indWeights[groupStart + j] * t);
	        }
	        norm = l2_norm(z, groupLen[i]);
          if (norm > 0) {
            uOp = (1 - lambda2[0] * grpWeights[i] * sqrt((double)g) * t / norm);
            if (uOp < 0) {
              uOp = 0;
            }
          } else {
            uOp = 0;
          }

    		  for (int j = 0; j < g; j++) {
    		    U[j] = uOp * z[j];
    		    G[j] = 1/t * (beta[groupStart + j] - U[j]);
          }

    		  for (int k = 0; k < n; k++) {
    		    etaNew[k] = eta[k];
    			  for(int j = 0; j < g; j++) {
    			    etaNew[k] -= t*G[j] * X[k + n * (groupStart + j)];
    			  }
    		  }

		      Lnew = linNegLogLikelihoodCalc(nrow, etaNew, y);
    		  sqNormG = pow(l2_norm(G, g), 2);
    		  iProd = 0;
    		  for (int j = 0; j < g; j++) {
    		    iProd += grad[j] * G[j];
    		  }

	        diff = Lold - Lnew - t * iProd + t/2 * sqNormG;
          t *= gamma[0];

        } while (diff < 0);

        // Nesterov momentum update
        t /= gamma[0];
        count++; check = 0;
        for (int j = 0; j < g; j++) {
	        check += fabs(theta[j + rangeGroupInd[i]] - U[j]);
	        for (int k = 0; k < nrow[0]; k++) {
	          eta[k] = eta[k] -
	            X[k + nrow[0] * (j + rangeGroupInd[i])]*beta[groupStart + j];
	        }
    		  beta[groupStart + j] = U[j] +
    		    count/(count+3) * (U[j] - theta[groupStart + j]);
    		  theta[groupStart + j] = U[j];
	        for(int k = 0; k < nrow[0]; k++) {
	           eta[k] = eta[k] +
	             X[k + nrow[0] * (groupStart + j)]*beta[groupStart + j];
	        }
        }

  	  } while (count < innerIter[0] && check > thresh[0]);

  	  Free(z);
  	  Free(U);
  	  Free(G);
  	  Free(betaNew);
    }
    Free(grad);
  }
  Free(theta);
  Free(etaNew);
  Free(etaNull);
}

void fit_gaussian(double *alpha, double *allbeta, double *beta, int *rangeLambdaInd,
         double *X, double* y, int *nrow, int *ncol,
         int *numGroup, int *rangeGroupInd, int *groupLen, double *lambda,
         double *lambda1, double *lambda2, int *nlam, int *innerIter,
         int *outerIter, double *thresh, double *outerThresh,
         double *eta, double *gamma, int *betaIsZero, double *step,
         double *grpWeights, double *indWeights){

  int n = nrow[0];
  int p = ncol[0];
  int g = numGroup[0];

  int* isActive = Calloc(g, int);
  double* ldot = Calloc(n, double);
  double* oldBeta = Calloc(p, double);

  // Move down regularization path
  for (int i = 0; i < nlam[0]; i++) {

    lambda1[0] = lambda[i] * alpha[0];
    lambda2[0] = lambda[i] * (1 - alpha[0]);

    int check, count = 0;
    do { // Iterate until convergence

      for (int j = 0; j < p; j++) {
 	      oldBeta[j] = beta[j];
      }

      // Update solution
      linSolver(X, y, nrow, ncol, numGroup, beta, rangeGroupInd,
                groupLen, lambda1, lambda2, innerIter, thresh, ldot,
                gamma, eta, betaIsZero, isActive, step,
                grpWeights, indWeights);

        // Check convergence criteria
        count++; check = 0;
        for (int j = 0; j < p; j++) {
     	    check += fabs(oldBeta[j] - beta[j]);
        }
    } while (count < outerIter[0] && check > outerThresh[0]);

    int start = rangeLambdaInd[i];
    int end = rangeLambdaInd[i] + ncol[0];
    for (int j = start; j < end; j++) {
      allbeta[j] = beta[j-start];
    }
  }

  Free(ldot);
  Free(isActive);
  Free(oldBeta);
}

// -----------------------------------------------------------------------------

void pCalc(int *nrow, double *eta, double *prob) {
  for (int i = 0; i < nrow[0]; i++) {
    if (eta[i] > 10) {
      prob[i] = 1;
    } else if (eta[i] < -10) {
      prob[i] = 0;
    } else {
      prob[i] = exp(eta[i]) / (1 + exp(eta[i]));
    }
  }
}

void logitGradCalc(int *nrow, double *prob, int *y, double *ldot) {
  for (int i = 0; i < nrow[0]; i++) {
    ldot[i] = -(y[i] - prob[i]) / nrow[0];
  }
}

double logitNegLogLikelihoodCalc(int *nrow, double *prob, int *y) {
  double logLik = 0;

  for(int i = 0; i < nrow[0]; i++)
    {
logLik = logLik + y[i] * log(prob[i]) + (1 - y[i]) * log(1 - prob[i]);
  }

  return -logLik/nrow[0];
}

void betaZeroSolve(int *nrow, double *betaZero, double *eta, double *prob,
                   double *thresh, int *innerIter, int *y) {
  double diff = 10;
  double num = 0;
  double denom = 0;
  int count = 0;

  while(pow(diff,2) > pow(thresh[0],2) && count < innerIter[0])
    {
pCalc(nrow, eta, prob);
diff = 0;

for(int i = 0; i < nrow[0]; i++)
  {
    num = num + y[i] - prob[i];
    denom = denom + prob[i] * (1 - prob[i]);
  }
diff = num/denom;
betaZero[0] = betaZero[0] + diff;

for(int i = 0; i < nrow[0]; i++)
  {
    eta[i] = eta[i] + diff;
  }
    }
}

void logitSolver(double *X, int *y, int *nrow, int *ncol,
                 int *numGroup, double *beta, int *rangeGroupInd,
                 int *groupLen, double *lambda1, double *lambda2,
                 int *innerIter, double *thresh, double *ldot,
                 double *gamma, double *eta, int* betaIsZero, int* isActive,
                 double *prob, double *betaZero, double *step,
                 double *grpWeights, double *indWeights) {

  // Intermediate quantities
  int n = nrow[0];
  int p = ncol[0];
  double* theta = Calloc(p, double);
  double* etaNew = Calloc(n, double);
  double* etaNull = Calloc(n, double);


  double zeroCheck = 0;
  double check = 0;
  int count = 0;
  double t = step[0];
  double diff = 1;
  double norm = 0;
  double uOp = 0;
  double Lnew = 0;
  double Lold = 0;
  double sqNormG = 0;
  double iProd = 0;

for (int i = 0; i < numGroup[0]; i++) {

  int g = groupLen[i];
  double* grad = Calloc(g, double);

  // Setting up null gradient calc to check if group is 0
  for (int j = rangeGroupInd[i]; j < rangeGroupInd[i] + groupLen[i]; j++) {
    for (int k = 0; k < nrow[0]; k++) {
      etaNull[k] = eta[k] - X[k + nrow[0] * j] * beta[j];
    }
  }

  // Calculating Null Gradient
  pCalc(nrow, etaNull, prob);
  logitGradCalc(nrow, prob, y, ldot);

  // double *grad = NULL;
  // grad = new double[groupLen[i]];

  for (int j = 0; j < groupLen[i]; j++) {
    grad[j] = crossprod(X, ldot, nrow[0], j + rangeGroupInd[i]);
    grad[j] = soft_threshold(grad[j], lambda1[0] * indWeights[(j + rangeGroupInd[i])]);
  }

  zeroCheck = l2_norm(grad, groupLen[i]);

  if (zeroCheck <= pow(lambda2[0] * grpWeights[i],2) * groupLen[i]) {
    if (betaIsZero[i] == 0) {
      for (int k = 0; k < nrow[0]; k++) {
        for (int j = rangeGroupInd[i]; j < rangeGroupInd[i] + groupLen[i]; j++) {
          eta[k] = eta[k] - X[k + nrow[0] * j] * beta[j];
        }
      }
    }
    betaIsZero[i] = 1;
    for(int j = 0; j < groupLen[i]; j++) {
      beta[j + rangeGroupInd[i]] = 0;
    }
  } else {
    if (isActive[i] == 0) {
      //groupChange = 1;
    }
    isActive[i] = 1;

    for (int k = 0; k < ncol[0]; k++) {
      theta[k] = beta[k];
    }

	  betaIsZero[i] = 0;
    double* z = Calloc(g, double);
    double* U = Calloc(g, double);
    double* G = Calloc(g, double);
    double* betaNew = Calloc(g, double);

	  count = 0;
	  do {
      count++;
      pCalc(nrow, eta, prob);
      logitGradCalc(nrow, prob, y ,ldot);

      for (int j = 0; j < groupLen[i]; j++) {
        grad[j] = crossprod(X, ldot, nrow[0], j + rangeGroupInd[i]);
      }

      // diff = -1;
      //	      t = 0.5;
      pCalc(nrow, eta, prob);
      Lold = logitNegLogLikelihoodCalc(nrow, prob, y);

      // Back-tracking
      do {
        for(int j = 0; j < groupLen[i]; j++) {
          z[j] = beta[j + rangeGroupInd[i]] - t * grad[j];
          z[j] = soft_threshold(z[j], lambda1[0] * indWeights[(j + rangeGroupInd[i])] * t);
        }

        norm = l2_norm(z, groupLen[i]);

  		  if (norm != 0) {
  		    uOp = (1 - lambda2[0] * grpWeights[i]*sqrt( (double) groupLen[i])*t/norm);
  		  } else{
  		    uOp = 0;
  		  }
  		  if (uOp < 0) {
          uOp = 0;
        }

  		  for (int j = 0; j < groupLen[i]; j++) {
  		    U[j] = uOp*z[j];
  		    G[j] = 1/t *(beta[j + rangeGroupInd[i]] - U[j]);
        }

  		  for (int k = 0; k < nrow[0]; k++) {
  		    etaNew[k] = eta[k];
  			  for(int j = 0; j < groupLen[i]; j++) {
  			    etaNew[k] = etaNew[k] - t*G[j] * X[k + nrow[0]*(rangeGroupInd[i] + j)];
  			  }
  		  }

  		  pCalc(nrow, etaNew, prob);
  		  Lnew = logitNegLogLikelihoodCalc(nrow, prob, y);

  		  sqNormG = pow(l2_norm(G, groupLen[i]), 2);
  		  iProd = 0;
  		  for (int j = 0; j < groupLen[i]; j++) {
  		    iProd = iProd + grad[j] * G[j];
  		  }

        diff = Lold - Lnew - t * iProd + t/2 * sqNormG;
        t = t * gamma[0];
      } while (diff < 0);

      t = t / gamma[0];
      check = 0;
      for (int j = 0; j < groupLen[i]; j++) {
        check = check + fabs(theta[j + rangeGroupInd[i]] - U[j]);
        for (int k = 0; k < nrow[0]; k++) {
          eta[k] = eta[k] - X[k + nrow[0] * (j + rangeGroupInd[i])]*beta[j + rangeGroupInd[i]];
        }
  		  beta[j + rangeGroupInd[i]] = U[j] + count/(count+3) * (U[j] - theta[j + rangeGroupInd[i]]);
  		  theta[j + rangeGroupInd[i]] = U[j];

        for(int k = 0; k < nrow[0]; k++) {
          eta[k] = eta[k] + X[k + nrow[0] * (j + rangeGroupInd[i])]*beta[j + rangeGroupInd[i]];
        }
      }
	  } while (count < innerIter[0] && check > thresh[0]);

 	  Free(z);
	  Free(U);
	  Free(G);
	  Free(betaNew);
  }
  Free(grad);
}

  betaZeroSolve(nrow, betaZero, eta, prob, thresh, innerIter, y);

  Free(theta);
  Free(etaNew);
  Free(etaNull);
}

void fit_binomial(double *alpha, double *allbeta, double *beta, int *rangeLambdaInd,
         double *X, int* y, int *nrow, int *ncol,
         int *numGroup, int *rangeGroupInd, int *groupLen, double *lambda,
         double *lambda1, double *lambda2, int *nlam, int *innerIter,
         int *outerIter, double *thresh, double *outerThresh,
         double *eta, double *gamma, int *betaIsZero, double *betaZero,
         double *step, double *grpWeights, double *indWeights,
         double *intercepts) {

  int n = nrow[0];
  int p = ncol[0];
  int g = numGroup[0];

  int* isActive = Calloc(g, int);
  double oldBetaZero = 0;
  double* ldot = Calloc(n, double);
  double* prob = Calloc(n, double);
  double* oldBeta = Calloc(p, double);

  for (int i = 0; i < nlam[0]; i++) {

    lambda1[0] = lambda[i] * alpha[0];
    lambda2[0] = lambda[i] * (1 - alpha[0]);
    betaZero[0] = intercepts[i];

    int check, count = 0;
    do { // Iterate until convergence

      oldBetaZero = betaZero[0];
      for (int j = 0; j < p; j++) {
   	    oldBeta[j] = beta[j];
      }

      // Update solution
      logitSolver(X, y, nrow, ncol, numGroup, beta, rangeGroupInd,
                  groupLen, lambda1, lambda2, innerIter, thresh, ldot,
                  gamma, eta, betaIsZero, isActive, prob,
                  betaZero, step, grpWeights, indWeights);

      // Check convergence criteria
      count++; check = fabs(oldBetaZero - betaZero[0]);
      for (int j = 0; j < p; j++) {
   	    check += fabs(oldBeta[j] - beta[j]);
      }
    } while (count < outerIter[0] && check > outerThresh[0]);

    intercepts[i] = betaZero[0];
    if (i < (nlam[0] - 1)) {
      intercepts[i + 1] = intercepts[i];
    }
    int start = rangeLambdaInd[i];
    int end = rangeLambdaInd[i] + ncol[0];
    for (int j = start; j < end; j++) {
      allbeta[j] = beta[j-start];
    }
  }

  Free(ldot);
  Free(isActive);
  Free(oldBeta);
}
