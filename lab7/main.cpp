

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>

class NavierStokes {
private:
  double DELTA; // Rozmiar kroku siatki
  double RHO;   // Gęstość cieczy
  double MU;    // Lepkość cieczy
  double Q;     // Strumień objętości
  int NX;       // Liczba punktów w kierunku x
  int NY;       // Liczba punktów w kierunku y
  int I1;       // Indeks przeszkody w x
  int J1;       // Indeks przeszkody w y
  int J2;
  int IT_MAX; // Maksymalna liczba iteracji

  std::vector<std::vector<double>> PSI;  // Funkcja strumienia
  std::vector<std::vector<double>> ZETA; // Funkcja wirowości
  std::vector<std::vector<double>> U;    // Funkcja wirowości
  std::vector<std::vector<double>> V;    // Funkcja wirowości

  double q_out() {
    return Q *
           (pow(DELTA * NY, 3) - pow(DELTA * J1, 3) -
            3 * DELTA * J1 * pow(DELTA * NY, 2) +
            3 * pow(DELTA * J1, 2) * DELTA * NY) /
           pow(DELTA * NY, 3);
  }

  bool is_obstacle(int i, int j) {
    if ((i == 0 && j >= J1 && j <= NY) || // Ściana A
        (j == NY) ||                      // Ściana B
        (i == NX) ||                      // Ściana C
        (i >= I1 && j == 0) ||            // Ściana D
        (i == I1 && j >= 0 && j <= J1) || // Ściana E
        (i >= 0 && i <= I1 && j == J1))   // Ściana F
    {
      return true;
    }
    return false;
  }
  // Warunki brzegowe dla PSI
  void psiBoundaries() {
    // Ściana A (wejście)
    for (int j = J1; j <= NY; j++) {

      PSI[0][j] = Q / (2.0 * MU) *
                  (std::pow(j * DELTA, 3) / 3.0 -
                   std::pow(j * DELTA, 2) / 2.0 * (J1 * DELTA + NY * DELTA) +
                   j * DELTA * J1 * DELTA * NY * DELTA);
    }

    // Ściana C (wyjście)
    for (int j = 0; j <= NY; j++) {
      PSI[NX][j] = q_out() / (2.0 * MU) *
                       (std::pow(j * DELTA, 3) / 3.0 -
                        std::pow(j * DELTA, 2) / 2.0 * NY * DELTA) +
                   Q / (12.0 * MU) * std::pow(J1 * DELTA, 2) *
                       (-J1 * DELTA + 3.0 * NY * DELTA);
    }

    // Ściana B
    for (int i = 1; i < NX; i++) {
      PSI[i][NY] = PSI[0][NY];
    }

    // Ściana D
    for (int i = I1; i < NX; i++) {
      PSI[i][0] = PSI[0][J1];
    }

    // Ściana E
    for (int j = 1; j <= J1; j++) {
      PSI[I1][j] = PSI[0][J1];
    }

    // Ściana F
    for (int i = 1; i <= I1; i++) {
      PSI[i][J1] = PSI[0][J1];
    }
  }

  // Warunki brzegowe dla ZETA
  void zetaBoundaries() {
    // Ściana A (wejście)
    for (int j = J1; j <= NY; j++) {
      ZETA[0][j] = Q / (2.0 * MU) * (2 * j * DELTA - J1 * DELTA - NY * DELTA);
    }

    // Ściana C (wyjście)
    for (int j = 0; j <= NY; j++) {
      ZETA[NX][j] = q_out() / (2.0 * MU) * (2 * j * DELTA - NY * DELTA);
    }

    // Ściana B
    for (int i = 1; i < NX; i++) {
      ZETA[i][NY] = 2.0 / std::pow(DELTA, 2) * (PSI[i][NY - 1] - PSI[i][NY]);
    }

    // Ściana D
    for (int i = I1 + 1; i < NX; i++) {
      ZETA[i][0] = 2.0 / std::pow(DELTA, 2) * (PSI[i][1] - PSI[i][0]);
    }

    // Ściana E
    for (int j = 1; j < J1; j++) {
      ZETA[I1][j] = 2.0 / std::pow(DELTA, 2) * (PSI[I1 + 1][j] - PSI[I1][j]);
    }

    // Ściana F
    for (int i = 1; i <= I1; i++) {
      ZETA[i][J1] = 2.0 / std::pow(DELTA, 2) * (PSI[i][J1 + 1] - PSI[i][J1]);
    }

    // Wierzchołek E/F
    ZETA[I1][J1] = 0.5 * (ZETA[I1 - 1][J1] + ZETA[I1][J1 - 1]);
  }

public:
  NavierStokes(double delta, double rho, double mu, double q, int nx, int ny,
               int i1, int j1, int iter_max)
      : DELTA(delta), RHO(rho), MU(mu), Q(q), NX(nx), NY(ny), I1(i1), J1(j1),
        J2(j1 + 2), IT_MAX(iter_max) {
    PSI = std::vector<std::vector<double>>(NX + 1,
                                           std::vector<double>(NY + 1, 0));
    ZETA = std::vector<std::vector<double>>(NX + 1,
                                            std::vector<double>(NY + 1, 0));
    U = std::vector<std::vector<double>>(NX + 1,
                                         std::vector<double>(NY + 1, 0));
    V = std::vector<std::vector<double>>(NX + 1,
                                         std::vector<double>(NY + 1, 0));
  }

  void solve() {
    psiBoundaries();

    double omega;
    double zeta_tmp;
    double gamma;

    for (int it = 1; it <= IT_MAX; it++) {
      omega = it < 2000 ? 0 : 1;
      for (int i = 1; i < NX; i++) {
        for (int j = 1; j < NY; j++) {
          if (!is_obstacle(i, j)) {
            PSI[i][j] = 0.25 * (PSI[i + 1][j] + PSI[i - 1][j] + PSI[i][j + 1] +
                                PSI[i][j - 1] - DELTA * DELTA * ZETA[i][j]);

            double zeta_tmp = 0.25 * (ZETA[i + 1][j] + ZETA[i - 1][j] +
                                      ZETA[i][j + 1] + ZETA[i][j - 1]);

            ZETA[i][j] =
                omega ? (zeta_tmp - RHO / (16.0 * MU) *
                                        ((PSI[i][j + 1] - PSI[i][j - 1]) *
                                             (ZETA[i + 1][j] - ZETA[i - 1][j]) -
                                         (PSI[i + 1][j] - PSI[i - 1][j]) *
                                             (ZETA[i][j + 1] - ZETA[i][j - 1])))
                      : zeta_tmp;

            U[i][j] = (PSI[i][j + 1] - PSI[i][j - 1]) / (2.0 * DELTA);
            V[i][j] = -(PSI[i + 1][j] - PSI[i - 1][j]) / (2.0 * DELTA);
          }
        }
      }
      zetaBoundaries();
      gamma = 0.0;
      for (int i = 1; i < NX; ++i) {
        gamma += PSI[i + 1][J2] + PSI[i - 1][J2] + PSI[i][J2 + 1] +
                 PSI[i][J2 - 1] - 4.0 * PSI[i][J2] -
                 DELTA * DELTA * ZETA[i][J2];
      }
      if (it % 100 == 0) {
        std::cout << "Iteracja " << it << ", Gamma: " << gamma << std::endl;
      }
    }
  }

  void saveDataToFile(const std::string &filename) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
      std::cerr << "Error opening file: " << filename << std::endl;
      return;
    }

    // #   plik z danymi: x, y, psi(x,y), ksi(x,y), u(x,y), v(x,y)
    for (int i = 1; i < NX; ++i) {
      for (int j = 1; j < NY; ++j) {

        outFile << i * DELTA << " " << j * DELTA << " " << PSI[i][j] << " "
                << ZETA[i][j] << " " << U[i][j] << " " << V[i][j] << "\n";
      }
      outFile << "\n";
    }
    outFile.close();
  }
};

int main() {
  NavierStokes zad1(0.01, 1, 1, -1000, 200, 90, 50, 55, 20000);

  NavierStokes zad2(0.01, 1, 1, -4000, 200, 90, 50, 55, 20000);
  NavierStokes zad3(0.01, 1, 1, 4000, 200, 90, 50, 55, 20000);

  zad1.solve();
  zad1.saveDataToFile("data/wyniki_Qm1000.dat");

  zad2.solve();
  zad2.saveDataToFile("data/wyniki_Qm4000.dat");
  zad3.solve();
  zad3.saveDataToFile("data/wyniki_Q4000.dat");

  return 0;
}
