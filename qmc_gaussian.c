#include "qmc_stats.h"
#include "hydrogen.h"

double gaussian(double *r) {
      const double norm_gauss = 1.0 / pow(2.0 * acos(-1.0), 1.5);
      return norm_gauss * exp(-0.5 * (r[0]*r[0] + r[1]*r[1] + r[2]*r[2]));
}

double gaussian_montecarlo(double a, const int nmax) {
      double r[3], energy, norm, w;

      norm = 0.0;
      energy = 0.0;
      for (int i = 0; i < nmax; ++i) {
	    random_gauss(r, 3);
	    w = psi(a, r, 3) * psi(a, r, 3) / gaussian(r);
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
	    x[i] = gaussian_montecarlo(a, nmax);
      }

      ave_error(x, nruns, obs);

      printf("E = %.5lf +/- %.5lf\n", obs[0], obs[1]);
      
      return 0;
}

