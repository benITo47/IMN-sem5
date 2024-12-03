extern "C" {
#include "mgmres.h"
}
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

const double DELTA = 0.1;

double rho(int i, int j, double xmax, double ymax, double sigma) {
  return exp(-pow((DELTA * i - 0.25 * xmax) / sigma, 2) -
             pow((DELTA * j - 0.5 * ymax) / sigma, 2)) -
         exp(-pow((DELTA * i - 0.75 * xmax) / sigma, 2) -
             pow((DELTA * j - 0.5 * ymax) / sigma, 2));
}

int calculateJ(int l, int nx) { return l / (nx + 1); }

int calculateI(int l, int nx) { return l - calculateJ(l, nx) * (nx + 1); }

double chooseEpsilon(int l, int nx, double epsilon1, double epsilon2) {
  if (calculateI(l, nx) <= nx / 2.0) {
    return epsilon1;
  } else {
    return epsilon2;
  }
}

void poissonEquation(int nx, int ny, double epsilon1, double epsilon2,
                     double v1, double v2, double v3, double v4,
                     bool rhoPresent, bool printSpace, std::string filename) {
  const int N = (nx + 1) * (ny + 1);
  double xMax = DELTA * nx;
  double yMax = DELTA * ny;

  double sigma = xMax / 10.0;

  std::vector<double> a =
      std::vector(5 * N, 0.0); // Niezerowe element macierzowe
  std::vector<double> b = std::vector(N, 0.0);
  std::vector<double> V = std::vector(N, 0.0);
  std::vector<int> ja = std::vector(N * 5, 0); // Informacje o numerach kolumn
  std::vector<int> ia = std::vector(
      N + 1, -1); // Wskazniki do elementow rozpoczynajacych dany wiersz

  int k = -1;
  for (int l = 0; l < N; l++) {
    int brzeg = 0;
    double vb = 0;
    int i = calculateI(l, nx);
    int j = calculateJ(l, nx);
    // lewy brzeg
    if (i == 0) {
      brzeg = 1;
      vb = v1;
    }
    // gorny brzeg
    if (j == ny) {
      brzeg = 1;
      vb = v2;
    }
    // prawy brzeg
    if (i == nx) {
      brzeg = 1;
      vb = v3;
    }
    // dolny brzeg
    if (j == 0) {
      brzeg = 1;
      vb = v4;
    }

    if (brzeg == 1) {
      b[l] = vb;
    } else if (rhoPresent) {
      b[l] = -rho(i, j, xMax, yMax, sigma);
    }

    ia[l] = -1;
    // Lewa skrajan przekatna
    if (l - nx - 1 >= 0 && brzeg == 0) {
      k++;
      if (ia[l] < 0)
        ia[l] = k;
      a[k] = chooseEpsilon(l, nx, epsilon1, epsilon2) / pow(DELTA, 2);
      ja[k] = l - nx - 1;
    }
    // Poddiagonala
    if (l - 1 >= 0 && brzeg == 0) {
      k++;
      if (ia[l] < 0)
        ia[l] = k;
      a[k] = chooseEpsilon(l, nx, epsilon1, epsilon2) / pow(DELTA, 2);
      ja[k] = l - 1;
    }
    // diagonala
    k++;
    if (ia[l] < 0) {
      ia[l] = k;
    }
    if (brzeg == 0) {
      a[k] = -(2 * chooseEpsilon(l, nx, epsilon1, epsilon2) +
               chooseEpsilon(l + 1, nx, epsilon1, epsilon2) +
               chooseEpsilon(l + nx + 1, nx, epsilon1, epsilon2)) /
             pow(DELTA, 2);
    } else {
      a[k] = 1;
    }
    ja[k] = l;

    // Naddiagonala
    if (l < N && brzeg == 0) {
      k++;
      a[k] = chooseEpsilon(l + 1, nx, epsilon1, epsilon2) / pow(DELTA, 2);
      ja[k] = l + 1;
    }

    // prawa skarjan przekatna
    if (l < N - nx - 1 && brzeg == 0) {
      k++;
      a[k] = chooseEpsilon(l + nx + 1, nx, epsilon1, epsilon2) / pow(DELTA, 2);
      ja[k] = l + nx + 1;
    }
  }

  int nz_num = k + 1;
  ia[N] = nz_num;

  pmgmres_ilu_cr(N, nz_num, ia.data(), ja.data(), a.data(), V.data(), b.data(),
                 500, 500, 1e-8, 1e-8);

  if (printSpace) {
    std::ofstream fA("../data/matrix_A.dat");
    for (int l = 0; l < 5 * N; l++) {
      fA << l << " " << std::fixed << std::setprecision(6) << a[l] << "\n";
    }

    std::ofstream fb("../data/vector_b.dat");
    for (int l = 0; l < N; l++) {
      fb << l << " " << calculateI(l, nx) << " " << calculateJ(l, nx) << " "
         << std::fixed << std::setprecision(6) << b[l] << "\n";
    }
  } else {

    std::ofstream f("../data/" + filename);
    double limit = 0.0;
    for (int l = 0; l < N; l++) {
      if (DELTA * calculateJ(l, nx) > limit) {
        f << "\n";
      }
      f << std::fixed << std::setprecision(6) << DELTA * calculateJ(l, nx)
        << " " << DELTA * calculateI(l, nx) << " " << V[l] << "\n";
      limit = DELTA * calculateJ(l, nx);
    }
  }
}
int main() {
  poissonEquation(4, 4, 1, 1, 10, -10, 10, -10, false, true, "");

  // Dla ε1 = ε2 = 1, V1 = V3 = −V2 = −V4 = 10, ρ1=p2 = 0
  //
  // dla nx = ny = 50
  poissonEquation(50, 50, 1, 1, 10, -10, 10, -10, false, false, "n_50.dat");

  // dla nx = ny = 100
  poissonEquation(100, 100, 1, 1, 10, -10, 10, -10, false, false, "n_100.dat");

  // dla nx = ny = 200
  poissonEquation(200, 200, 1, 1, 10, -10, 10, -10, false, false, "n_200.dat");

  // Dla nx = ny = 100, V1 = V3 = V2 = V4 = 0,
  //
  // dla ε1 = 1 i ε2 = 1
  poissonEquation(100, 100, 1, 1, 0, 0, 0, 0, true, false, "e_1_1.dat");

  // dla ε1 = 1 i ε2 = 2
  poissonEquation(100, 100, 1, 2, 0, 0, 0, 0, true, false, "e_1_2.dat");

  // dla ε1 = 1 i ε2 = 10
  poissonEquation(100, 100, 1, 10, 0, 0, 0, 0, true, false, "e_1_10.dat");

  return 0;
}
