#include "hydrogen.h"

int main() {
      const double a[6] = {0.1, 0.2, 0.5, 1., 1.5, 2.};
      const int sizex = 50;
      double x[sizex], dx;

      dx = 10.0 / (sizex - 1);
      for (int i = 0; i < sizex; ++i) {
	    x[i] = -5.0 + (i - 1) * dx;
      }

      FILE * fil;

      fil = fopen("./data", "w+");

      for (int i; i < 6; ++i) {
	    for (int j = 0; j < sizex; ++j) {
		  fprintf(fil, "%lf %lf\n", x[j], e_loc(a[i], &x[j], 1)); 
	    }
	    fprintf(fil, "\n\n");
      }

      fclose(fil);

      return 0;
}
