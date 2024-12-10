#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

const int nx = 400;
const int ny = 90;
const int i1 = 200;
const int i2 = 210;
const int j_1 = 50;
const double delta = 0.01;
const double sigma = 10 * delta;
const double xA = 0.45;
const double yA = 0.45;
const int IT_MAX = 10000;

void calculate_v_field(double psi[nx + 1][ny + 1], double vx[nx + 1][ny + 1],
                       double vy[nx + 1][ny + 1]) {
  for (int i = 1; i < nx; i++) {
    for (int j = 1; j < ny; j++) {
      vx[i][j] = (psi[i][j + 1] - psi[i][j - 1]) / (2 * delta);
      vy[i][j] = -(psi[i + 1][j] - psi[i - 1][j]) / (2 * delta);
    }
  }
  for (int i = i1; i <= i2; i++) {
    for (int j = 0; j <= j_1; j++) {
      vx[i][j] = vy[i][j] = 0;
    }
  }
  for (int i = 1; i < nx; i++) {
    vx[i][0] = vy[i][ny] = 0;
  }
  for (int j = 0; j <= ny; j++) {
    vx[0][j] = vx[1][j];
    vx[nx][j] = vx[nx - 1][j];
  }
  fstream p1;
  p1.open("plots/V.txt", ios::out);
  if (p1.good()) {
    for (int i = 0; i <= nx; i++) {
      for (int j = 0; j <= ny; j++) {
        p1 << i << "\t" << j << "\t" << vx[i][j] << "\t" << vy[i][j] << "\n";
      }
      p1 << "\n";
    }
  }
  p1.close();
}

double v_max(double vx[nx + 1][ny + 1], double vy[nx + 1][ny + 1]) {
  double vmax = sqrt(pow(vx[0][0], 2) + pow(vy[0][0], 2));
  for (int i = 0; i < nx + 1; i++) {
    for (int j = 0; j < ny + 1; j++) {
      double v = sqrt(pow(vx[i][j], 2) + pow(vy[i][j], 2));
      if (v > vmax) {
        vmax = v;
      }
    }
  }
  return vmax;
}

void solve_ad(double D, double psi[nx + 1][ny + 1], double vx[nx + 1][ny + 1],
              double vy[nx + 1][ny + 1], double u0[nx + 1][ny + 1],
              double u1[nx + 1][ny + 1], double delta_T, fstream &p1) {
  char p = '0';
  for (int i = 0; i <= nx; i++) {
    for (int j = 0; j <= ny; j++) {
      u0[i][j] = exp(-(pow(i * delta - xA, 2) + pow(j * delta - yA, 2)) /
                     (2 * pow(sigma, 2))) /
                 (2 * M_PI * pow(sigma, 2));
    }
  }
  for (int it = 0; it <= IT_MAX; it++) {
    for (int i = 0; i <= nx; i++) {
      for (int j = 0; j <= ny; j++) {
        u1[i][j] = u0[i][j];
      }
    }

    for (int k = 1; k <= 20; k++) {
      for (int i = 0; i <= nx; i++) {
        for (int j = 1; j <= ny; j++) {
          if (i >= i1 && i <= i2 && j <= j_1) {
            continue;
          } else if (i == 0) {
            u1[i][j] =
                (1 / (1 + ((2 * D * delta_T) / pow(delta, 2)))) *
                (u0[i][j] -
                 (delta_T / 2) * vx[i][j] *
                     (((u0[i + 1][j] - u0[nx][j]) / (2 * delta)) +
                      (u1[i + 1][j] - u1[nx][j]) / (2 * delta)) -
                 (delta_T / 2) * vy[i][j] *
                     ((u0[i][j + 1] - u0[i][j - 1]) / (2 * delta) +
                      (u1[i][j + 1] - u1[i][j - 1]) / (2 * delta)) +
                 (delta_T / 2) * D *
                     ((u0[i + 1][j] + u0[nx][j] + u0[i][j + 1] + u0[i][j - 1] -
                       4 * u0[i][j]) /
                          pow(delta, 2) +
                      (u1[i + 1][j] + u1[nx][j] + u1[i][j + 1] + u1[i][j - 1]) /
                          pow(delta, 2)));
          } else if (i == nx) {
            u1[i][j] =
                (1 / (1 + ((2 * D * delta_T) / pow(delta, 2)))) *
                (u0[i][j] -
                 (delta_T / 2) * vx[i][j] *
                     (((u0[0][j] - u0[i - 1][j]) / (2 * delta)) +
                      (u1[0][j] - u1[i - 1][j]) / (2 * delta)) -
                 (delta_T / 2) * vy[i][j] *
                     ((u0[i][j + 1] - u0[i][j - 1]) / (2 * delta) +
                      (u1[i][j + 1] - u1[i][j - 1]) / (2 * delta)) +
                 (delta_T / 2) * D *
                     ((u0[0][j] + u0[i - 1][j] + u0[i][j + 1] + u0[i][j - 1] -
                       4 * u0[i][j]) /
                          pow(delta, 2) +
                      (u1[0][j] + u1[i - 1][j] + u1[i][j + 1] + u1[i][j - 1]) /
                          pow(delta, 2)));
          } else {
            u1[i][j] = (1 / (1 + ((2 * D * delta_T) / pow(delta, 2)))) *
                       (u0[i][j] -
                        (delta_T / 2) * vx[i][j] *
                            (((u0[i + 1][j] - u0[i - 1][j]) / (2 * delta)) +
                             (u1[i + 1][j] - u1[i - 1][j]) / (2 * delta)) -
                        (delta_T / 2) * vy[i][j] *
                            ((u0[i][j + 1] - u0[i][j - 1]) / (2 * delta) +
                             (u1[i][j + 1] - u1[i][j - 1]) / (2 * delta)) +
                        (delta_T / 2) * D *
                            ((u0[i + 1][j] + u0[i - 1][j] + u0[i][j + 1] +
                              u0[i][j - 1] - 4 * u0[i][j]) /
                                 pow(delta, 2) +
                             (u1[i + 1][j] + u1[i - 1][j] + u1[i][j + 1] +
                              u1[i][j - 1]) /
                                 pow(delta, 2)));
          }
        }
      }
    }
    for (int i = 0; i <= nx; i++) {
      for (int j = 0; j <= ny; j++) {
        u0[i][j] = u1[i][j];
      }
    }
    double c = 0.0;
    double x_sr = 0.0;
    for (int i = 0; i <= nx; i++) {
      for (int j = 0; j <= ny; j++) {
        c += u0[i][j];
        x_sr += (i * delta) * u0[i][j];
      }
    }
    c *= pow(delta, 2);
    x_sr *= pow(delta, 2);
    p1 << delta_T * it << "\t" << c << "\t" << x_sr << "\n";

    fstream p2;
    if (it == 2000 || it == 4000 || it == 6000 || it == 8000 || it == 10000) {
      p++;
      string name = "";
      if (D == 0.0) {
        name += "plots/graph";
        name += p;
        name += "_00.txt";
      } else {
        name += "plots/graph";
        name += p;
        name += "_01.txt";
      }
      p2.open(name, ios::out);
      if (p2.good()) {
        for (int i = 0; i <= nx; i++) {
          for (int j = 0; j <= ny; j++) {
            p2 << i << "\t" << j << "\t" << u1[i][j] << "\n";
          }
          p2 << "\n";
        }
      }
      p2.close();
    }
  }
}

int main() {
  double psi[nx + 1][ny + 1];
  double vx[nx + 1][ny + 1];
  double vy[nx + 1][ny + 1];
  fstream p, p1, p2;
  int i, j;
  double psi2;
  p.open("psi.dat", ios::in);
  if (p.good()) {
    while (!p.eof() && p >> i >> j >> psi2) {
      psi[i][j] = psi2;
    }
  }
  calculate_v_field(psi, vx, vy);

  double vmax = v_max(vx, vy);
  double delta_T = delta / (4 * vmax);
  double u0[nx + 1][ny + 1];
  double u1[nx + 1][ny + 1];

  p1.open("plots/ct_00.txt", ios::out);
  if (p1.good()) {
    solve_ad(0.0, psi, vx, vy, u0, u1, delta_T, p1);
  }

  p2.open("plots/ct_01.txt", ios::out);
  if (p2.good()) {
    solve_ad(0.1, psi, vx, vy, u0, u1, delta_T, p2);
  }

  p.close();
  return 0;
}
