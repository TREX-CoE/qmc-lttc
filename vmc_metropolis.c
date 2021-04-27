#include "hydrogen.h"
#include "qmc_stats.h"

void variational_montecarlo(double a, const int nmax, double dt,
			   double *energy, double *accept) {
      int n_accept;
      double r_old[3], r_new[3], psi_old, psi_new;
      double d_old[3], d_new[3], d2_old, d2_new;
      double rnd[3], fact_a, fact_b, sq_dt, q, u;

      sq_dt = sqrt(dt);

      // Initial position
      random_gauss(r_old, 3);

      drift(a, r_old, d_old, 3);
      d2_old = 0.0;
      for (int j = 0; j < 3; ++j) {
	    d2_old += d_old[j] * d_old[j];
      }

      psi_old = psi(a, r_old, 3);

      *energy = 0.0;
      *accept = 0.0;
      n_accept = 0;
      for (int i = 0; i < nmax; ++i) {
	    // Compute and accumulate the local energy
	    *energy += e_loc(a, r_old, 3);

	    // Compute new position (correct variance of sampled Gaussian)
	    random_gauss(rnd, 3);
	    for (int j = 0; j < 3; ++j) {
		  r_new[j] = r_old[j] + dt * d_old[j] + rnd[j] * sq_dt;
	    }

	    // New WF and acceptance probability
            drift(a, r_new, d_new, 3);
	    d2_new = 0.0;
            for (int j = 0; j < 3; ++j) {
		  d2_new += d_new[j] * d_new[j];
            }

	    psi_new = psi(a, r_new, 3);

	    // Compute the ratio of probabilities q
	    fact_b = 0.0;
	    for (int j = 0; j < 3; ++j) {
		  fact_b += (d_new[j] + d_old[j]) * (r_new[j] - r_old[j]);
	    }
	    fact_a = 0.5 * dt * (d2_new - d2_old) + fact_b;

	    q = psi_new / psi_old;
	    q = exp(-fact_a) * q * q;

	    u = (double) rand() / RAND_MAX;

	    if (u <= q) {
		  n_accept += 1;
		  for (int j = 0; j < 3; ++j) {
			r_old[j] = r_new[j];
			d_old[j] = d_new[j];
		  }
		  d2_old = d2_new;
		  psi_old = psi_new;
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

      double ene[nruns], acc[nruns], obs[2];

      for (int i = 0; i < nruns; ++i) {
	    variational_montecarlo(a, nmax, dt, &ene[i], &acc[i]);
      }

      ave_error(ene, nruns, obs);
      printf("E = %.5lf +/- %.5lf\n", obs[0], obs[1]);

      ave_error(acc, nruns, obs);
      printf("A = %.5lf +/- %.5lf\n", obs[0], obs[1]);
      
      return 0;
}
