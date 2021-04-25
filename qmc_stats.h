#include <math.h>
#include <stdio.h>

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
