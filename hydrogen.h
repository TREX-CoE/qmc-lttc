#include <math.h>
#include <stdio.h>

double potential(double *r, const int l) {
      double pot;

      pot = 0.0;
      for (int i = 0; i < l; ++i) {
	    pot += r[i] * r[i];
      }

      if (pot > 0.0) {
	   return -1.0 / sqrt(pot); 
      }
      return 1e6;

}

double psi(double a, double *r, const int l) {
      double psival, rnorm;

      psival = 0.0;
      rnorm = 0.0;
      for (int i = 0; i < l; ++i) {
	    rnorm += r[i]*r[i];
      }
      rnorm = sqrt(rnorm);

      psival = exp(-a * rnorm);

      return psival;
}

double kinetic(double a, double *r, const int l) {
      double rnorm;

      rnorm = 0.0;
      for (int i = 0; i < l; ++i) {
	    rnorm += r[i]*r[i];
      }
      rnorm = sqrt(rnorm);

      return a * (1 / rnorm - 0.5 * a);
      
}

double e_loc(double a, double *r, const int l) {
      return kinetic(a, r, l) + potential(r, l);
}

void drift(double a, double *r, double *d, const int l) {
      double rnorm, fact;

      rnorm = 0.0;
      for (int i = 0; i < l; ++i) {
	    rnorm += r[i]*r[i];
      }
      rnorm = sqrt(rnorm);

      fact = -a / rnorm;
      for (int i = 0; i < l; ++i) {
	    d[i] = r[i] * fact;
      }
}
