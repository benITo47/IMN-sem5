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
void relaksacja() {

  std::vector<std::vector<double>> V(NX + 1, std::vector<double>(NY + 1, 0.0));

  for (int j = 0; j <= NY; ++j) {
    double y = j * DY;
    V[0][j] = VB1(y); // Lewy brzeg
  }
  for (int i = 0; i <= NX; ++i) {
    double x = i * DX;
    V[i][NY] = VB2(x); // Górny brzeg
  }
  for (int j = 0; j <= NY; ++j) {
    double y = j * DY;
    V[NX][j] = VB3(y); // Prawy brzeg
  }
  for (int i = 0; i <= NX; ++i) {
    double x = i * DX;
    V[i][0] = VB4(x); // Dolny brzeg
  }
  int iteration = 0;
  for (int k : KS) {
    double S_prev = 0.0;
    double S_next;

    // Pliki do zapisu danych funkcji całkowej i mapy potencjału
    std::ofstream file("data/s_" + std::to_string(k) + ".txt");
    std::ofstream potential_file("data/v_" + std::to_string(k) + ".txt");

    while (true) {
      iteration++;

      // Iteracja relaksacji
      for (int i = k; i < NX; i += k) {
        for (int j = k; j < NY; j += k) {
          V[i][j] =
              0.25 * (V[i + k][j] + V[i - k][j] + V[i][j + k] + V[i][j - k]);
        }
      }
      // Zapis do pliku
      for (int i = 0; i <= NX; i += k) {
        for (int j = 0; j <= NY; j += k) {
          potential_file << i * DX << "\t" << j * DY << "\t" << V[i][j] << "\n";
        }
        potential_file << "\n"; // Odstep dla GNUPLOTa
      }

      // Obliczanie funkcji całkowej

      S_next = 0.0;
      for (int i = 0; i <= NX - k; i += k) {
        for (int j = 0; j <= NY - k; j += k) {
          double term1_x =
              (V[i + k][j] - V[i][j] + V[i + k][j + k] - V[i][j + k]) /
              (2.0 * DX * k);
          double term1_y =
              (V[i][j + k] - V[i][j] + V[i + k][j + k] - V[i + k][j]) /
              (2.0 * DY * k);
          S_next +=
              0.5 * k * k * 0.2 * 0.2 * (term1_x * term1_x + term1_y * term1_y);
        }
      }

      // Zapisujemy numer iteracji i wartość całki
      file << iteration << "\t" << S_next << "\n";

      // Sprawdzenie warunku stopu
      if (iteration > 0 && std::abs(S_next - S_prev) / S_prev < TOL)
        break;

      S_prev = S_next;
    }

    // Zagęszczanie siatki
    if (k > 1) {
      int k2 = k / 2;
      for (int i = 0; i <= NX - k; i += k) {
        for (int j = 0; j <= NY - k; j += k) {
          V[i + k2][j + k2] =
              0.25 * (V[i][j] + V[i + k][j] + V[i][j + k] + V[i + k][j + k]);
          if (i != NX - k)
            V[i + k][j + k2] = 0.5 * (V[i + k][j] + V[i + k][j + k]);

          if (j != NY - k)
            V[i + k2][j + k] = 0.5 * (V[i][j + k] + V[i + k][j + k]);
          if (j != 0)
            V[i + k2][j] = 0.5 * (V[i][j] + V[i + k][j]);
          if (i != 0)
            V[i][j + k2] = 0.5 * (V[i][j] + V[i][j + k]);
        }
      }
    }
  }
}

int main() {
  relaksacja();
  return 0;
}
