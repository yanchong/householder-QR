#include <cmath>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "linearalgebra.h"

/* ----------------------- householder ----------------------- */
/*  Given a matrix A of dimension m by n (with n <= m) and 
    arrays v_i of dimension m-i, for i = 1, ..., m - 1, 
    respectively, this algorithm computes n reflection vectors
    and the factor R of a full QR decomposition of A, where R 
    is a m by n upper triangular matrix. The n reflection 
    vectors are stored in the arrays v_1, ..., v_n and the 
    columns of A are overwritten by the columns of R.

    Input variables:
a: pointer to array of arrays, the ith array of
which should correspond to the ith column of the 
matrix A. During the algorithm, the columns of R 
will overwrite the columns of A.
v: pointer to array of arrays in which the ith 
reflection vector of dimension m - i will be 
stored.
m: number of rows of A.
n: number of columns of A.

Features: The number of flops for this implementation is
~ 2 * m * n^2 - (2/3) * n^3 and requires O(1) additional 
memory.                                                    */

void printqt(std::vector<std::vector<double>>& qt, int m) {
  for(int x = 0; x < m; x++) {
    for(int y = 0; y < m; y++) {
      if (fabs(qt[x][y]) < 1e-10) {
        printf("%9.6g ", 0.0);
      } else {
        printf("%9.6g ", qt[x][y]);
      }
    }
    printf("\n");
  }
  printf("\n");
}
void householder (double ** a, double ** v, double ** q, int m, int n) {
  int i, j;
  double vnorm, vTa, vpartdot;

  for(i = 0; i < n; i++) {
    /* set v[i] equal to subvector a[i][i : m] */
    partialvec_copy(a[i], v[i], m - i, i);

    /* vpartdot = ||v[i]||^2 - v[i][0] * v[i][0]; since vpartdot 
       is unaffected by the change in v[i][0], storing this value 
       prevents the need to recalculate the entire norm of v[i] 
       after updating v[i][0] in the following step              */
    vpartdot = partialdot_product(v[i], v[i], m - i, 1);

    /* set v[i][0] = v[i][0] + sign(v[i][0]) * ||v[i]|| */
    if(v[i][0] < 0) {
      v[i][0] -= sqrt(v[i][0] * v[i][0] + vpartdot);
    }
    else {
      v[i][0] += sqrt(v[i][0] * v[i][0] + vpartdot);
    }

    /* normalize v[i] */
    vnorm = sqrt(v[i][0] * v[i][0] + vpartdot);
    scalar_div(v[i], vnorm, m - i, v[i]);
    double scale = v[i][0];
    scale = 1;
    scalar_div(v[i], scale, m - i, v[i]);

    for(j = i; j < n; j++) {
      /* set a[j][i:m] = a[j][i:m] - 2 * (v[i]^T a[j][i:m]) * v[i] */
      vTa = subdot_product(a[j], v[i], m - i, i);
      vTa *= (2 * scale * scale);
      partialscalar_sub(v[i], vTa, m - i, i, a[j]);
    }
    std::vector<std::vector<double>> qt(m);
    for (int x = 0; x < m; x++) {
      qt[x].resize(m);
      for (int y = 0; y < m; y++) {
        if (x == y) {
          qt[x][y] = 1;
        } else {
          qt[x][y] = 0;
        }
      }
    }
    for (int x = i; x < m; x++) {
      for (int y = i; y < m; y++) {
        qt[x][y] -= 2 * v[i][x-i] * v[i][y-i] * scale * scale;
      }
    }
    std::vector<std::vector<double>> qnew(m);
    for (int x = 0; x < m; x++) {
      qnew[x].resize(m);
    }
    for (int x = 0; x < m; x++) {
      for (int y = 0; y < m; y++) {
        double sum = 0.0;
        for (int z = 0; z < m; z++) {
          sum += q[x][z] * qt[z][y];
        }
        qnew[x][y] = sum;
      }
    }
    for (int x = 0; x < m; x++) {
      for (int y = 0; y < m; y++) {
        q[x][y] = qnew[x][y];
      }
    }
  }
}


int main () {
  int i, j, m, n;
  double x;

  /* let user set the dimension of matrix A */
  //printf("Enter the dimension m (where A is a m by n matrix): ");
  m = 4;
  //scanf("%i", &m);
  //printf("Enter the dimension n (where A is a m by n matrix): ");
  //scanf("%i", &n);
  n = 4;

  /* check if m < n */
  if(m < n) {
    printf("For a successful factorization, this implementation "
        "requires n <= m.\nTerminating program.\n");
    return 0;
  }

  /* allocate memory for A and vectors v */
  double ** a = new double*[n];
  double ** v = new double*[n];
  for(i = 0; i < n; i++) {
    a[i] = new double[m];
    v[i] = new double[m - i];
  }
  /* initialize the values in matrix A */
  for(i = 0; i < n; i++) {
    for(j = 0; j < m; j++) {
      if(j < i) {
        a[i][j] = 0;
      }
      else {
        a[i][j] = j - i + 1; // this choice of values was arbitrary
      }
    }
  }
  // q is column first 2d-array
  double ** q = new double*[m];
  for (i = 0; i < m; i++) {
    q[i] = new double[m];
  }
  for (i = 0; i < m; i++) {
    for (j = 0; j < m; j++) {
      if (i == j) {
        q[i][j] = 1;
      } else {
        q[i][j] = 0;
      }
    }
  }


  /* print the matrix A before calling houheholder */
  printf("A = \n");
  for(i = 0; i < m; i++) {
    for(j = 0; j < n; j++) {
      if (fabs(a[j][i]) < 1e-10) {
        printf("%9.6g ", 0.0);
      } else {
        printf("%9.6g ", a[j][i]);
      }
    }
    printf("\n");
  }
  printf("\n");

  householder(a, v, q, m, n);

  /* print the matrix R (stored in A) after calling houheholder */
  printf("Q = \n");
  for(i = 0; i < m; i++) {
    for(j = 0; j < m; j++) {
      if (fabs(q[i][j]) < 1e-10) {
        printf("%9.6g ", 0.0);
      } else {
        printf("%9.6g ", q[i][j]);
      }
    }
    printf("\n");
  }
  printf("\n");

  /* print the matrix R (stored in A) after calling houheholder */
  printf("R = \n");
  for(i = 0; i < m; i++) {
    for(j = 0; j < n; j++) {
      if (fabs(a[j][i]) < 1e-10) {
        printf("%9.6g ", 0.0);
      } else {
        printf("%9.6g ", a[j][i]);
      }
    }
    printf("\n");
  }
  printf("\n");

  printf("QR = \n");
  for (i = 0; i < m; i++) {
    for (j = 0; j < n; j++) {
      double sum = 0;
      for (int k = 0; k < m; k++) {
        sum += q[i][k] * a[j][k];
      }
      if (fabs(sum) < 1e-10) {
        printf("%9.6g ", 0.0);
      } else {
        printf("%9.6g ", sum);
      }
    }
    printf("\n");
  }
  printf("\n");

  /* print the vectors v after calling householder */
  for(i = 0; i < n; i++) {
    printf("v[%i] = ", i);
    for(j = 0; j < m - i; j++) {
      if (fabs(v[i][j]) < 1e-10) {
        printf("%9.6g ", 0.0);
      } else {
        printf("%9.6g ", v[i][j]);
      }
    }
    printf("\n");
  }
  printf("\n");

  /* print numerical evidence that v's are normalized */
  printf("Numerical verification that v_1, ..., v_%i are "
      "normalized:\n", n);
  for(i = 1; i < n; i++) {
    x = dot_product(v[i - 1], v[i - 1], m - i + 1);
    printf("||v[%i]|| = %lg, ", i, x);
    if(i % 5 == 0) {
      printf("\n");
    }
  }
  x = dot_product(v[n - 1], v[n - 1], m - n + 1);
  printf("||v[%i]|| = %lg.", n, x);
  if(n % 5 != 0) printf("\n");
  printf("\n");

  /* free memory */
  for(i = 0; i < n; i++) {
    delete[] a[i];
    delete[] v[i];
  }
  delete[] a;
  delete[] v;
  for (i = 0; i < m; i++) {
    delete [] q[i];
  }
  delete [] q;
  return 0;
}
