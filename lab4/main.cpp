
#include <chrono>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <ostream>
#include <vector>

const double EPSILON = 1;
const double DELTA = 0.1;
const double N_X = 150;
const double N_Y = 100;
const double V_1 = 10;
const double V_2 = 0;
const double X_MAX = N_X * DELTA;
const double Y_MAX = N_Y * DELTA;
const double SIGMA_X = 0.1 * X_MAX;
const double SIGMA_Y = 0.1 * Y_MAX;
const double TOL = 1e-8;

// Funkcja do obliczania gęstości ładunku ρ(x, y)
double gestosc_RO(double x, double y) {
  double ro1 = exp(-(pow((x - 0.35 * X_MAX), 2) / pow(SIGMA_X, 2)) -
                   (pow(y - 0.5 * Y_MAX, 2)) / (pow(SIGMA_Y, 2)));
  double ro2 = exp(-(pow((x - 0.65 * X_MAX), 2) / pow(SIGMA_X, 2)) -
                   (pow(y - 0.5 * Y_MAX, 2)) / (pow(SIGMA_Y, 2)));
  return ro1 + ro2;
}

void relaksacja_globalna() {
  std::vector<double> omega_gs = {0.6, 1.0};

  for (double omega_g : omega_gs) {
    // Inicjalizacja tablic potencjałów
    std::vector<std::vector<double>> v_stare(N_X + 1,
                                             std::vector<double>(N_Y + 1, 0.0));
    std::vector<std::vector<double>> v_nowe(N_X + 1,
                                            std::vector<double>(N_Y + 1, 0.0));

    // Ustawienie warunków brzegowych Dirichleta na dolnym brzegu (V1) i górnym
    // (V2)
    for (int i = 0; i <= N_X; ++i) {
      v_stare[i][0] = V_1;
      v_stare[i][N_Y] = V_2;
      v_nowe[i][0] = V_1;
      v_nowe[i][N_Y] = V_2;
    }

    // Inicjalizacja zmiennych do obliczeń
    std::vector<std::vector<double>> err(N_X, std::vector<double>(N_Y, 0.0));
    double s_prev;
    double s_next = 0.0;
    int it = 0;

    // Otwarcie plików do zapisów
    std::ofstream f_s("dane/globalna_s_" +
                      std::to_string(omega_g).substr(0, 3) + ".txt");
    std::ofstream f_v("dane/globalna_v_" +
                      std::to_string(omega_g).substr(0, 3) + ".txt");
    std::ofstream f_err("dane/globalna_err_" +
                        std::to_string(omega_g).substr(0, 3) + ".txt");

    std::cout << "Calculating with global relaxation for omega_g = " << omega_g
              << "\n";

    while (true) {
      ++it;

      // Obliczanie nowych wartości potencjału V w każdym węźle (oprócz brzegów)
      for (int i = 1; i < N_X; ++i) {
        for (int j = 1; j < N_Y; ++j) {
          v_nowe[i][j] =
              0.25 *
              (v_stare[i + 1][j] + v_stare[i - 1][j] + v_stare[i][j + 1] +
               v_stare[i][j - 1] +
               DELTA * DELTA * gestosc_RO(i * DELTA, j * DELTA) / EPSILON);
        }
      }

      // Nałożenie warunków brzegowych Neumanna (po bokach)
      for (int j = 1; j < N_Y; ++j) {
        v_nowe[0][j] = v_nowe[1][j];         // Lewy brzeg
        v_nowe[N_X][j] = v_nowe[N_X - 1][j]; // Prawy brzeg
      }

      // Aktualizacja potencjału V przy użyciu relaksacji globalnej
      for (int i = 0; i <= N_X; ++i) {
        for (int j = 0; j <= N_Y; ++j) {
          v_stare[i][j] =
              (1.0 - omega_g) * v_stare[i][j] + omega_g * v_nowe[i][j];
        }
      }

      // Obliczanie całki funkcjonalnej S dla warunku stopu
      s_prev = s_next;
      s_next = 0.0;
      for (int i = 0; i < N_X; ++i) {
        for (int j = 0; j < N_Y; ++j) {
          double dx = (v_stare[i + 1][j] - v_stare[i][j]) / DELTA;
          double dy = (v_stare[i][j + 1] - v_stare[i][j]) / DELTA;
          s_next += DELTA * DELTA *
                    (0.5 * dx * dx + 0.5 * dy * dy -
                     gestosc_RO(i * DELTA, j * DELTA) * v_stare[i][j]);
        }
      }

      f_s << it << " " << s_next << "\n";
      std::cout << "\r[" << it << "] " << s_next << std::flush;

      // Sprawdzenie warunku stopu na podstawie zmiany całki S
      if (std::abs((s_next - s_prev) / s_prev) < TOL) {
        break;
      }
    }

    std::cout << std::endl;

    // Obliczanie błędu numerycznego w każdym węźle i zapisywanie wyników
    for (int i = 1; i < N_X; ++i) {
      for (int j = 1; j < N_Y; ++j) {
        err[i][j] = (v_nowe[i + 1][j] - 2.0 * v_nowe[i][j] + v_nowe[i - 1][j]) /
                        (DELTA * DELTA) +
                    (v_nowe[i][j + 1] - 2.0 * v_nowe[i][j] + v_nowe[i][j - 1]) /
                        (DELTA * DELTA) +
                    gestosc_RO(i * DELTA, j * DELTA) / EPSILON;

        // Zapis błędu oraz potencjału do plików wynikowych
        f_err << i * DELTA << " " << j * DELTA << " " << err[i][j] << "\n";
        f_v << i * DELTA << " " << j * DELTA << " " << v_nowe[i][j] << "\n";
      }
      f_err << "\n";
      f_v << "\n";
    }

    f_s.close();
    f_v.close();
    f_err.close();
  }
}

void relaksacja_lokalna() {
  std::vector<double> omega_ls = {1.0, 1.4, 1.8, 1.9};

  // Pętla dla różnych wartości omega_L
  for (double omega_l : omega_ls) {
    std::vector<std::vector<double>> V(N_X + 1,
                                       std::vector<double>(N_Y + 1, 0.0));

    // Ustawienie warunków brzegowych Dirichleta na dolnym brzegu (V1) i górnym
    // brzegu (V2)
    for (int i = 0; i <= N_X; ++i) {
      V[i][0] = V_1;
      V[i][N_Y] = V_2;
    }

    std::vector<std::vector<double>> err(N_X, std::vector<double>(N_Y, 0.0));
    double s_prev;
    double s_next = 0.0;
    int it = 0;

    std::ofstream f_s("dane/lokalna_s_" + std::to_string(omega_l).substr(0, 3) +
                      ".txt");

    while (true) {
      ++it;

      // Aktualizacja wartości potencjału w każdym węźle (oprócz brzegów)
      for (int i = 1; i < N_X; ++i) {
        for (int j = 1; j < N_Y; ++j) {
          V[i][j] =
              (1.0 - omega_l) * V[i][j] +
              omega_l * 0.25 *
                  (V[i + 1][j] + V[i - 1][j] + V[i][j + 1] + V[i][j - 1] +
                   DELTA * DELTA * gestosc_RO(i * DELTA, j * DELTA) / EPSILON);
        }
      }

      // Nałożenie warunków brzegowych von Neumanna (po bokach)
      for (int j = 1; j < N_Y; ++j) {
        V[0][j] = V[1][j];
        V[N_X][j] = V[N_X - 1][j];
      }

      // Obliczanie całki funkcjonalnej S dla warunku stopu
      s_prev = s_next;
      s_next = 0.0;
      for (int i = 0; i < N_X; ++i) {
        for (int j = 0; j < N_Y; ++j) {
          double dx = (V[i + 1][j] - V[i][j]) / DELTA;
          double dy = (V[i][j + 1] - V[i][j]) / DELTA;
          s_next += DELTA * DELTA *
                    (0.5 * dx * dx + 0.5 * dy * dy -
                     gestosc_RO(i * DELTA, j * DELTA) * V[i][j]);
        }
      }

      // Zapis wyników do pliku oraz wyświetlenie iteracji
      f_s << it << " " << s_next << "\n";
      std::cout << "\r[" << it << "] " << s_next << std::flush;

      // Sprawdzenie warunku stopu na podstawie zmiany całki S
      if (std::abs((s_next - s_prev) / s_prev) < TOL) {
        break;
      }
    }

    std::cout << std::endl;

    f_s.close();
  }
}
int main() {
  auto start = std::chrono::high_resolution_clock::now();
  relaksacja_lokalna();
  relaksacja_globalna();
  auto end = std::chrono::high_resolution_clock::now();

  std::cout << "\n Executed in: "
            << std::chrono::duration_cast<std::chrono::seconds>(end - start)
            << std::endl;
  return 0;
}
