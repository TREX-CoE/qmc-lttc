#include "hydrogen.h"
#include "qmc_stats.h"
#include <stdlib.h>
#include <time.h>

double uniform_montecarlo(double a, const int nmax) {
      double r[3], energy, norm, w, rnd;

      norm = 0.0;
      energy = 0.0;
      for (int i = 0; i < nmax; ++i) {
	    for (int j = 0; j < 3; ++j) {
		  rnd = (double) rand() / RAND_MAX;
		  r[j] = -5.0 + 10.0 * rnd;
	    }
	    w = psi(a, r, 3) * psi(a, r, 3);
	    norm += w;
	    energy += w * e_loc(a, r, 3);
      }

      return energy / norm;

}

int main() {
      const double a = 1.2;
      const long nmax = 1e5;
      int nruns = 30;

      srand(time(NULL));

      double x[nruns], obs[2];

      for (int i = 0; i < nruns; ++i) {
	    x[i] = uniform_montecarlo(a, nmax);
      }

      ave_error(x, nruns, obs);

      printf("E = %.5lf +/- %.5lf\n", obs[0], obs[1]);
      
      return 0;
}

