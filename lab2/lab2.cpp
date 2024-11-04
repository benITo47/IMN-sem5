

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

const double beta = 0.001;
const int N = 500;
const double GAMMA = 0.1;
const double t_max = 100.0;
const double dt = 0.1;
const double u0 = 1.0;
const double TOL = 1e-6;
const int max_iter = 20;
const double alpha = (beta * N - GAMMA);

double f(double u) { return (alpha)*u - beta * u * u; }

// Metoda trapezów z iteracją Picarda
std::vector<double> metoda_trapezow_picarda() {
  std::vector<double> u_values; // Wektor wynikow
  double u_n = u0;              // Inicjalizacja warunkiem początkowym
  u_values.push_back(u_n);

  std::ofstream file("picard_iter.dat"); // otwiera plik
  // Glowna pętla co dt
  for (double t = dt; t <= t_max; t += dt) {
    double u_next =
        u_n; // Zaczynamy od ustalenia u_n+1 na u_n z poprzedniej iteracji

    int current_iteracje = 0;
    for (int iter = 0; iter < max_iter; ++iter) {
      double u_temp =
          u_n + (dt / 2.0) * (f(u_n) +
                            f(u_next)); // zgodnie ze wzorem 10 obliczamy u_n+1;
                            
      current_iteracje = iter;
      if (std::fabs(u_temp - u_next) <=
          TOL) // Obliczamy roznice, jesli mniejsza niz zalozona to przerywamy
               // bo dokladnosc jest zadowalajaca
        break;
      u_next = u_temp; // Przypisujemy wyliczone U_n+1
    }
    file << t << ", " << current_iteracje << "\n";
    u_values.push_back(u_next); // dodajemy obliczone u do wynikow
    u_n = u_next;               // ustawiamy U_n do kolejnej iteracji
  }

  file.close();
  return u_values;
}

// Metoda trapezow Picarda
std::vector<double> metoda_trapezow_newtona() {
  std::vector<double> u_values; // Wektor wynikow
  double u_n = u0;              // Inicjalizacja warunkiem poczatkowym
  u_values.push_back(u_n);

  std::ofstream file("newton_iter.dat"); // otwiera plik

  for (double t = dt; t <= t_max; t += dt) {
    double u_next = u_n; // ustalamy, u_n+1 na u_n z poprzedniej iteracji

    int current_iteracje = 0;
    for (int iter = 0; iter < max_iter; ++iter) {
      // definiujemy F i dF dla metody newtona
      double F = u_next - u_n - (dt / 2) * (f(u_n) + f(u_next));
      double dF = 1 - (dt / 2) * ((alpha)-2 * beta * u_next);
      double u_temp = u_next - F / dF;
      
      current_iteracje = iter;
      if (std::fabs(u_temp - u_next) <=
          TOL) // Obliczamy roznice miedzy iteracjami, jesli mniejsza niz
               // tolerancja to przerywamy bo dokladnosc jest zadowalajaca
        break;
      u_next = u_temp; // Aktualizujemy u_next dla kolejnej iteracji w petli
                       // wewnetrznej
    }
    file << t << ", " << current_iteracje << "\n";
    u_values.push_back(u_next); // dodajemy obliczone u_n+1 do wynikow
    u_n = u_next;               // ustalamy U_n dla kolejnej iteracji
  }

  file.close();
  return u_values;
}

std::vector<double> metoda_niejawna_RK() {
  std::vector<double> u_values; // Wektor wynikow
  double u_n = u0;              // Inicjalizacja warunkiem początkowym
  u_values.push_back(u_n);

  std::ofstream file("rk_iter.dat"); // otwiera plik

  // Wspolczynniki tablicy Butchera

  const double a11 = 0.25;
  const double a12 = 0.25 - sqrt(3) / 6.0;
  const double a21 = 0.25 + sqrt(3) / 6.0;
  const double a22 = 0.25;
  const double b1 = 0.5;
  const double b2 = 0.5;
  const double c1 = 0.5 - sqrt(3) / 6.0;
  const double c2 = 0.5 - sqrt(3) / 6.0;

  for (double t = dt; t <= t_max; t += dt) {
    double U1 = u_n, U2 = u_n; // ustalamy U1 i U2 na u_n z poprzedniej iteracji

    int current_iteracje = 0;
    for (int i = 0; i < max_iter; i++) {
      // Wyznaczamy F1 i F2
      double F1 = U1 - u_n - dt * (a11 * f(U1) + a12 * f(U2));
      double F2 = U1 - u_n - dt * (a21 * f(U1) + a22 * f(U2));
      // Obliczamy elementy Jakobianu
      double m11 = 1 - dt * a11 * (alpha - 2 * beta * U1);
      double m12 = -dt * a12 * (alpha - 2 * beta * U2);
      double m21 = -dt * a21 * (alpha - 2 * beta * U1);
      double m22 = 1 - dt * a22 * (alpha - 2 * beta * U2);

      // Obliczamy wspolczynnik macierzy
      double determinant = m11 * m22 - m12 * m21;

      // Jezeli mniejszy od toleracnji to przerywamy liczenie
      if (fabs(determinant) <= TOL)
        break;

      // Wyznaczamy dU1 i dU2 wykorzystujac obliczone wartosci;
      double dU1 = (m12 * F2 - m22 * F1) / determinant;
      double dU2 = (m21 * F1 - m11 * F2) / determinant;

      // Aktualizujemy U1 i U2
      U1 += dU1;
      U2 += dU2;

      // Przerywamy iteracje jezeli delty są mniejsze od tolerancji
      if (std::fabs(dU1) <= TOL && std::fabs(dU2) <= TOL)
        break;

      current_iteracje = i;
    }

    // Ustawiamy u_n+1 dla kolejnej iteracji
    u_n = u_n + dt * (b1 * f(U1) + b2 * f(U2));
    u_values.push_back(u_n);

    file << t << ", " << current_iteracje << "\n";
  }

  file.close();
  return u_values;
}
// Funkcja zapisująca wyniki do pliku
void zapisz_do_pliku(const std::string &filename,
                     const std::vector<double> &u_values) {
  std::ofstream file(filename); // otwiera plik
  double t = 0.0;
  for (double u : u_values) {
    double z = N - u;
    file << t << ", " << u << ", " << z << "\n"; // zapisujemy kolejno t,u,z
    t += dt;
  }
  file.close();
}

int main() {
  // Obliczanie wyników za pomocą  metod
  auto u_picarda = metoda_trapezow_picarda();
  auto u_newtona = metoda_trapezow_newtona();
  auto u_rk = metoda_niejawna_RK();
  // Zapis wyników do plików dat
  zapisz_do_pliku("wyniki_picarda.dat", u_picarda);
  zapisz_do_pliku("wyniki_newtona.dat", u_newtona);
  zapisz_do_pliku("wyniki_rk.dat", u_rk);
  return 0;
}
