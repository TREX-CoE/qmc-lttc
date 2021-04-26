#include "qmc_stats.h"

int main() {
      double rnd[3];
      random_gauss(rnd, 3);
      for (int i = 0; i < 3; ++i) {
	    printf("val %lf\n", rnd[i]);
      }

      return 0;
}
