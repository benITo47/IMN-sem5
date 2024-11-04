#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>

// Parametry początkowe
double X0 = 0.01;     // Początkowa pozycja x,
double V0 = 0;        // Początkowa prędkość,
double delta_T0 = 1;  // Początkowy Krok czasowy
double S = 0.75;      // Współczynnik bezpieczeństwa,
int P = 2;            // Rząd dokładności metody,
double T_MAX = 40;    // Jak długo będziemy symulować - maksymalny czas
double ALFA = 5;      // Parametr alfa - mówi jak "sztywny" jest problem
double TOL1 = 1e-2;   // Tolerancja pierwsza
double TOL2 = 1e-5;   // Tolerancja druga
double DELTA = 1e-10; // Minimalna dokładność

// Funkcja, która zwraca pochodną x - czyli po prostu v w tym przypadku
double dx_deriv(double x, double v) { return v; }

// Funkcja, która zwraca pochodną v
double dv_deriv(double x, double v) { return ALFA * (1 - x * x) * v - x; }

// Funkcja do kontroli kroku czasowego
double kontrola_kroku(double Ex, double Ev, double TOL, double dt) {
  return pow((S * TOL) / (std::max(fabs(Ex), fabs(Ev))), 1.0 / (P + 1)) * dt;
}

// Metoda trapezów - liczy x i v w kolejnym kroku, ale wymaga iteracji
std::pair<double, double> metodaTrapezow(double x, double v, double dt,
                                         double alfa) {
  double x0 = x, v0 = v, xn1 = x0, vn1 = v0, dx, dv;

  // Pętla iteracyjna Newtona - poprawia wyniki aż będą wystarczająco dobre
  do {
    double fx = xn1 - x0 - dt / 2.0 * (dx_deriv(x0, v0) + dx_deriv(xn1, vn1));
    double gx = vn1 - v0 - dt / 2.0 * (dv_deriv(x0, v0) + dv_deriv(xn1, vn1));

    double a11 = 1;
    double a12 = -dt / 2.0;
    double a21 = -dt / 2.0 * (-2 * ALFA * xn1 * vn1 - 1);
    double a22 = 1 - dt / 2.0 * ALFA * (1 - xn1 * xn1);

    double wyznacznik = a11 * a22 - a12 * a21;

    // Wyliczamy poprawki dx i dv
    dx = (-fx * a22 + gx * a12) / wyznacznik;
    dv = (-gx * a11 + fx * a21) / wyznacznik;

    xn1 += dx;
    vn1 += dv;
  } while (fabs(dx) > DELTA || fabs(dv) > DELTA);

  return {xn1, vn1}; // Zwracamy nowe wartości x i v
}

// Metoda RK2 - liczy x i v w kolejnym kroku,
std::pair<double, double> metodaRK2(double x, double v, double dt,
                                    double alfa) {
  // Liczymy współczynniki k1 i k2
  double k1x = dx_deriv(x, v);
  double k1v = dv_deriv(x, v);

  double k2x = dx_deriv(x + dt * k1x, v + dt * k1v);
  double k2v = dv_deriv(x + dt * k1x, v + dt * k1v);

  //  na podstawie k1 i k2 dostajemy nowe x i v
  double xn1 = x + dt / 2.0 * (k1x + k2x);
  double vn1 = v + dt / 2.0 * (k1v + k2v);

  return {xn1, vn1};
}

void rozwiaz(
    std::function<std::pair<double, double>(double, double, double, double)>
        funkcja,
    double tol, std::string filename) {

  // Inicjalizacja  parametrów
  std::vector<double> t = {0};
  std::vector<double> dt = {delta_T0};
  std::vector<double> x = {X0};
  std::vector<double> v = {V0};

  // Otwarcie pliku do zapisu
  std::ofstream file(filename);
  if (!file) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }

  // Pętla czasu - symulujemy aż do T_MAX
  while (t.back() < T_MAX) {
    // Robimy dwa kroki po dt
    std::pair<double, double> para1 =
        funkcja(x.back(), v.back(), dt.back(), ALFA);
    double x2_n1 = para1.first, v2_n1 = para1.second;

    std::pair<double, double> para2 = funkcja(x2_n1, v2_n1, dt.back(), ALFA);
    double x2_n2 = para2.first, v2_n2 = para2.second;

    // Robimy jeden krok 2 * dt, żeby porównać
    std::pair<double, double> para3 =
        funkcja(x.back(), v.back(), 2 * dt.back(), ALFA);
    double x1_n2 = para3.first, v1_n2 = para3.second;

    // Liczymy błędy Ex i Ev
    double Ex = (x2_n2 - x1_n2) / (pow(2, P) - 1);
    double Ev = (v2_n2 - v1_n2) / (pow(2, P) - 1);

    // Jak błędy są małe, to akceptujemy krok i zapisujemy wyniki
    if (std::max(fabs(Ex), fabs(Ev)) < tol) {
      t.push_back(t.back() + 2 * dt.back());
      x.push_back(x2_n2);
      v.push_back(v2_n2);
      file << std::fixed << std::setprecision(6) << t.back() << "\t"
           << dt.back() << "\t" << x.back() << "\t" << v.back() << "\n";
    }

    // Aktualizujemy krok czasowy zgodnie z kontrolą kroku
    dt.back() = kontrola_kroku(Ex, Ev, tol, dt.back());
    if (dt.back() < DELTA) { // Jeśli krok za mały, to kończymy
      break;
    }
  }

  // Zamykamy plik
  file.close();
}

int main() {
  rozwiaz(metodaTrapezow, TOL1, "trapezy_tol1.txt");
  rozwiaz(metodaTrapezow, TOL2, "trapezy_tol2.txt");

  rozwiaz(metodaRK2, TOL1, "RK2_tol1.txt");
  rozwiaz(metodaRK2, TOL2, "RK2_tol2.txt");

  return 0;
}
