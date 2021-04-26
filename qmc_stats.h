#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void ave_error(double *x, const int n, double obs[2]) {
      // obs[1] --> mean and obs[2] --> err
      double var;

      if (n == 1) {
	    obs[0] = x[0];
	    obs[1] = 0.0;
      }
      else {
	    for (int i = 0; i < n; ++i) {
		  obs[0] += x[i];
	    }
            obs[0] = obs[0] / n;

            var = 0.0;
            for (int i = 0; i < n; ++i) {
            	  var += (x[i] - obs[0]) * (x[i] - obs[0]);
            }
            var = var / (n - 1);

            obs[1] = sqrt(var / n);
      }

}

void random_gauss(double *z, const int n) {
      // Box Muller method
      const double two_pi = 2.0 * acos(-1.0);
      double u[n + 1];

      srand(time(NULL));

      for (int i = 0; i < n + 1; ++i) {
          u[i] = (double) rand() / RAND_MAX;
      }


      if (!(n & 1)) {
	    // n is even
	    for (int i = 0; i < n; i+=2) {
		  z[i] = sqrt(-2.0 * log(u[i]));
		  z[i] = z[i] * cos(two_pi * u[i + 1]);
		  z[i + 1] = z[i] * sin(two_pi * u[i + 1]);
	    }
      }
      else {
	    // n is odd
	    for (int i = 0; i < n - 1; i+=2) {
		  z[i] = sqrt(-2.0 * log(u[i]));
		  z[i] = z[i] * cos(two_pi * u[i + 1]);
		  z[i + 1] = z[i] * sin(two_pi * u[i + 1]);
	    }
	    z[n] = sqrt(-2.0 * log(u[n]));
	    z[n] = z[n] * cos(two_pi * u[n + 1]);
      }
}
