#include "hydrogen.h"
#include "qmc_stats.h"
#include <stdlib.h>
#include <time.h>

void metropolis_montecarlo(double a, const int nmax, double dt,
			   double *energy, double *accept) {
      int n_accept;
      double r_old[3], r_new[3], psi_old, psi_new, rnd;
      double aval, u;

      // Initial position
      for (int j = 0; j < 3; ++j) {
	    rnd = (double) rand() / RAND_MAX;
	    r_old[j] = dt * (2.0 * rnd - 1.0);
      }
      psi_old = psi(a, r_old, 3) * psi(a, r_old, 3);

      *energy = 0.0;
      *accept = 0.0;
      n_accept = 0;
      for (int i = 0; i < nmax; ++i) {
	    // Compute and accumulate the local energy
	    *energy += e_loc(a, r_old, 3);

	    // Compute new position
	    for (int j = 0; j < 3; ++j) {
		 rnd = (double) rand() / RAND_MAX;
		 r_new[j] = r_old[j] + dt * (2.0 * rnd - 1.0);
	    }

	    // New WF and acceptance probability
	    psi_new = psi(a, r_new, 3) * psi(a, r_new, 3);
	    aval = psi_new / psi_old;

	    u = (double) rand() / RAND_MAX;

	    if (u <= aval) {
		  for (int j = 0; j < 3; ++j) {
			r_old[j] = r_new[j];
		  }
		  psi_old = psi_new;
		  n_accept += 1;

	    }

      }
      *energy /= nmax;
      *accept = (double) n_accept / nmax;
}

int main() {
      const double a = 1.2;
      const double dt = 1.0;
      const long nmax = 1e5;
      int nruns = 30;

      srand(time(NULL));

      double x[nruns], y[nruns], obs[2];

      for (int i = 0; i < nruns; ++i) {
	    metropolis_montecarlo(a, nmax, dt, &x[i], &y[i]);
      }

      ave_error(x, nruns, obs);
      printf("E = %.5lf +/- %.5lf\n", obs[0], obs[1]);

      ave_error(y, nruns, obs);
      printf("A = %.5lf +/- %.5lf\n", obs[0], obs[1]);
      
      return 0;
}
