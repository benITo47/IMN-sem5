
#include "gsl/gsl_blas.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_math.h"
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

class Crank_Nicolson {
private:
  static const int NX = 40;
  static const int NY = 40;
  static const int N = (NX + 1) * (NY + 1);
  static const int DELTA = 1;
  static const int DELTA_T = 1;
  static const int TA = 40;
  static const int TB = 0;
  static const int TC = 30;
  static const int TD = 0;
  static const int IT_MAX = 2000;

  gsl_matrix *A;
  gsl_matrix *B;
  gsl_vector *c;
  gsl_vector *T;
  gsl_permutation *p;
  gsl_vector *d;

  double kB = 0.1;
  double kD = 0.6;

  int l(int i, int j) { return i + j * (NX + 1); }

public:
  Crank_Nicolson() {
    A = gsl_matrix_calloc(N, N);
    B = gsl_matrix_calloc(N, N);
    T = gsl_vector_calloc(N);
    c = gsl_vector_calloc(N);
    p = gsl_permutation_calloc(N);
    d = gsl_vector_calloc(N);

    intialize_arrays();
  }

  ~Crank_Nicolson() {
    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_vector_free(c);
    gsl_vector_free(T);
    gsl_vector_free(d);
    gsl_permutation_free(p);
  }

  void intialize_arrays() {

    // Internal area
    for (int i = 1; i < NX; i++) {
      for (int j = 1; j < NY; j++) {
        double a = DELTA_T / (2 * pow(DELTA, 2));

        gsl_matrix_set(A, l(i, j), l(i, j) - NX - 1, a);
        gsl_matrix_set(A, l(i, j), l(i, j) - 1, a);
        gsl_matrix_set(A, l(i, j), l(i, j) + 1, a);
        gsl_matrix_set(A, l(i, j), l(i, j) + NX + 1, a);
        gsl_matrix_set(A, l(i, j), l(i, j), -4 * a - 1);

        gsl_matrix_set(B, l(i, j), l(i, j) - NX - 1, -a);
        gsl_matrix_set(B, l(i, j), l(i, j) - 1, -a);
        gsl_matrix_set(B, l(i, j), l(i, j) + 1, -a);
        gsl_matrix_set(B, l(i, j), l(i, j) + NX + 1, -a);
        gsl_matrix_set(B, l(i, j), l(i, j), 4 * a - 1);
      }
    }

    initalizeEdges();

    // Initial conditions
    for (int j = 0; j <= NY; j++) {
      gsl_vector_set(T, l(0, j), TA);
      gsl_vector_set(T, l(NX, j), TC);
      for (int i = 1; i < NX; i++) {
        gsl_vector_set(T, l(i, j), 0);
      }
    }
  }

  void initalizeEdges() {

    // Dirichlet  conditions (left and right edges)
    for (int i = 0; i <= NX; i += NX) {
      for (int j = 0; j <= NY; j++) {
        gsl_matrix_set(A, l(i, j), l(i, j), 1);
        gsl_matrix_set(B, l(i, j), l(i, j), 1);
        gsl_vector_set(c, l(i, j), 0);
      }
    }

    // Neumann  conditions (upper edge for n+1)
    for (int i = 1; i < NX; i++) {
      gsl_matrix_set(A, l(i, NY), l(i, NY) - NX - 1, -1 / (kB * DELTA));
      gsl_matrix_set(A, l(i, NY), l(i, NY), 1 + 1 / (kB * DELTA));
      gsl_vector_set(c, l(i, NY), TB);
      for (int j = 0; j <= NY; j++) {
        gsl_matrix_set(B, l(i, NY), j, 0);
      }
    }

    // Neumann  conditions (lower edge for n+1 )
    for (int i = 1; i < NX; i++) {
      gsl_matrix_set(A, l(i, 0), l(i, 0), 1 + 1 / (kD * DELTA));
      gsl_matrix_set(A, l(i, 0), l(i, 0) + NX + 1, -1 / (kD * DELTA));
      gsl_vector_set(c, l(i, 0), TD);
      for (int j = 0; j <= NY; j++) {
        gsl_matrix_set(B, l(i, 0), j, 0);
      }
    }
  }
  void solve() {

    // Crank-Nicolson algorithm
    int signum;
    gsl_linalg_LU_decomp(A, p, &signum);

    for (int it = 0; it <= IT_MAX; it++) {
      gsl_blas_dgemv(CblasNoTrans, 1, B, T, 0, d);
      gsl_blas_daxpy(1, c, d);
      gsl_linalg_LU_solve(A, p, d, T);

      if (it == 100 || it == 200 || it == 500 || it == 1000 || it == 2000) {
        save_to_file(it);
      }
    }
  }

  void save_to_file(int it) {
    string name = "output/" + to_string(it) + ".dat";
    fstream p1;
    p1.open(name, ios::out);
    if (p1.good()) {
      for (int i = 0; i <= NX; i++) {
        for (int j = 0; j <= NY; j++) {
          p1 << i << "\t" << j << "\t" << gsl_vector_get(T, l(i, j)) << "\n";
        }
        p1 << "\n";
      }
    }
    p1.close();

    string name2 = "output/laplasjan_" + to_string(it) + ".dat";
    fstream p2;
    p2.open(name2, ios::out);
    if (p2.good()) {
      for (int i = 1; i <= NX - 1; i++) {
        for (int j = 1; j <= NY - 1; j++) {
          p2 << i << "\t" << j << "\t"
             << ((gsl_vector_get(T, l(i, j) + 1) -
                  2 * gsl_vector_get(T, l(i, j)) +
                  gsl_vector_get(T, l(i, j) - 1)) /
                 pow(DELTA, 2)) +
                    ((gsl_vector_get(T, l(i, j) + NX + 1) -
                      2 * gsl_vector_get(T, l(i, j)) +
                      gsl_vector_get(T, l(i, j) - NX - 1)) /
                     pow(DELTA, 2))
             << "\n";
        }
        p2 << "\n";
      }
    }
    p2.close();
  }
};

int main() {
  Crank_Nicolson simulation;
  simulation.solve();
  return 0;
}
