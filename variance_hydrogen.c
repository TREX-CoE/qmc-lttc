#include "hydrogen.h"

int main() {
      const double a[6] = {0.1, 0.2, 0.5, 1., 1.5, 2.};
      const int sizex = 50;
      double energy, r[3], sigma;
      double x[sizex], dx, w, wsum, dV, e_tmp;

      dx = 10.0 / (sizex - 1);
      for (int i = 0; i < sizex; ++i) {
	    x[i] = -5.0 + i * dx;
      }
      dV = dx * dx *dx;

      for (int l = 0; l < 6; ++l) {
	    energy = 0.0;
	    sigma = 0.0;
	    wsum = 0.0;
	    for (int i = 0; i < sizex; ++i) {
		  r[0] = x[i];
		  for (int j = 0; j < sizex; ++j) {
			r[1] = x[j];
			for (int k = 0; k < sizex; ++k) {
			      r[2] = x[k];
			      w = psi(a[l], r, 3) * psi(a[l], r, 3) * dV;
			      wsum += w;
			      e_tmp = e_loc(a[l], r, 3);
			      energy += w * e_tmp;
			      sigma += w * e_tmp * e_tmp;
			}
		  }
	    }
	    energy = energy / wsum;
	    sigma = sigma / wsum;
	    sigma -= energy * energy;
	    printf("a = %.2lf  E = %.5lf  Sigma^2 = %.5lf\n", a[l], energy, sigma); 
      }

      
      return 0;
}
