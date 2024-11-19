#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

const int NX = 128, NY = 128;                // Rozmiary siatki
const double DX = 0.2, DY = 0.2;             // Krok siatki
const double XMAX = NX * DX, YMAX = NY * DY; // Rozmiary obszaru
const double TOL = 1e-8;                     // Warunek stopu
std::vector<int> KS = {16, 8, 4, 2, 1};      // Kroki siatki

// Funkcje warunków brzegowych
double VB1(double y) { return sin(M_PI * y / YMAX); }
double VB2(double x) { return -sin(2 * M_PI * x / XMAX); }
double VB3(double y) { return sin(M_PI * y / YMAX); }
double VB4(double x) { return sin(2 * M_PI * x / XMAX); }

// Relaksacja wielosiatkowa
void multigridRelaxation() {

  std::vector<std::vector<double>> V(NX + 1, std::vector<double>(NY + 1, 0.0));

  // Ustawiamy warunki brzegowe
  for (int i = 0; i <= NX; ++i) {
    double x = i * DX;
    V[i][0] = VB4(x);  // Dolny brzeg
    V[i][NY] = VB2(x); // Górny brzeg
  }
  for (int j = 0; j <= NY; ++j) {
    double y = j * DY;
    V[0][j] = VB1(y);  // Lewy brzeg
    V[NX][j] = VB3(y); // Prawy brzeg
  }

  int iteration = 0;
  for (int k : KS) {
    int step = k;

    double S_prev = 0.0;
    double S_next;

    // Pliki do zapisu danych funkcji całkowej i mapy potencjału
    std::ofstream file("data/s_" + std::to_string(k) + ".txt");
    std::ofstream potential_file("data/v_" + std::to_string(k) + ".txt");

    while (true) {
      iteration++;

      // Iteracja relaksacji
      std::vector<std::vector<double>> V_old = V;
      for (int i = step; i <= NX - step; i += step) {
        for (int j = step; j <= NY - step; j += step) {
          V[i][j] = 0.25 * (V[i + step][j] + V[i - step][j] + V[i][j + step] +
                            V[i][j - step]);
        }
      }

      // Obliczanie funkcji całkowej

      S_next = 0.0;
      for (int i = 0; i <= NX - step; i += step) {
        for (int j = 0; j <= NY - step; j += step) {
          double term1_x = (V[i + step][j] - V[i][j] + V[i + step][j + step] -
                            V[i][j + step]) /
                           (2.0 * DX * step);
          double term1_y = (V[i][j + step] - V[i][j] + V[i + step][j + step] -
                            V[i + step][j]) /
                           (2.0 * DY * step);
          S_next += 0.5 * step * step * DX * DY *
                    (term1_x * term1_x + term1_y * term1_y);
        }
      }

      // Zapisujemy numer iteracji i wartość całki
      file << iteration << "\t" << S_next << "\n";

      // Sprawdzenie warunku stopu
      if (iteration > 0 && std::abs(S_next - S_prev) / S_prev < TOL)
        break;

      S_prev = S_next;
    }

    // Zapis mapy potencjału
    for (int i = 0; i <= NX; ++i) {
      for (int j = 0; j <= NY; ++j) {
        potential_file << i * DX << "\t" << j * DY << "\t" << V[i][j] << "\n";
      }
      potential_file << "\n"; // Dodatkowy odstęp dla czytelności
    }

    // Zagęszczanie siatki
    if (k > 1) {
      for (int i = 0; i < NX; i += step) {
        for (int j = 0; j < NY; j += step) {
          V[i + step / 2][j] = 0.5 * (V[i][j] + V[i + step][j]); // Pionowe
          V[i][j + step / 2] = 0.5 * (V[i][j] + V[i][j + step]); // Poziome
          V[i + step / 2][j + step / 2] =
              0.25 * (V[i][j] + V[i + step][j] + V[i][j + step] +
                      V[i + step][j + step]); // Nowe węzły
        }
      }
    }
    for (int i = 0; i <= NX; ++i) {
      double x = i * DX;
      V[i][0] = VB4(x);  // Dolny brzeg
      V[i][NY] = VB2(x); // Górny brzeg
    }
    for (int j = 0; j <= NY; ++j) {
      double y = j * DY;
      V[0][j] = VB1(y);  // Lewy brzeg
      V[NX][j] = VB3(y); // Prawy brzeg
    }
  }
}

int main() {
  multigridRelaxation();
  return 0;
}
