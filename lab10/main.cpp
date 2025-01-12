#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

class Varlete_Equation {
private:
  const int NX = 150;
  const int NT = 1000;
  const double DELTA = 0.1;
  const double DELTA_T = 0.05;
  const double X_A = 7.5;
  const double SIGMA = 0.5;
  const double ALFA;
  const double BETA;

  std::vector<double> U;
  std::vector<double> U_prev;
  std::vector<double> V;
  std::vector<double> Vp;
  std::vector<double> a;

  std::ofstream energy_file;
  std::ofstream displacement_file;

public:
  Varlete_Equation(double alfa, double beta) : ALFA(alfa), BETA(beta) {
    U = std::vector<double>(NX, 0.0);
    U_prev = std::vector<double>(NX, 0.0);
    V = std::vector<double>(NX, 0.0);
    Vp = std::vector<double>(NX, 0.0);
    a = std::vector<double>(NX, 0.0);

    std::filesystem::path data_dir("data");
    if (!std::filesystem::exists(data_dir)) {
      std::filesystem::create_directory(data_dir);
    }

    std::ostringstream alpha_stream, beta_stream;
    alpha_stream << std::fixed << std::setprecision(1) << ALFA;
    beta_stream << std::fixed << std::setprecision(1) << BETA;

    std::string energy_filename = "data/E_alpha_" + alpha_stream.str() +
                                  "_beta_" + beta_stream.str() + ".txt";
    std::string displacement_filename = "data/u_alpha_" + alpha_stream.str() +
                                        "_beta_" + beta_stream.str() + ".txt";

    energy_file.open(energy_filename);
    displacement_file.open(displacement_filename);

    if (!energy_file.is_open() || !displacement_file.is_open()) {
      throw std::runtime_error("Unable to open output files.");
    }
  }

  ~Varlete_Equation() {
    if (energy_file.is_open()) {
      energy_file.close();
    }
    if (displacement_file.is_open()) {
      displacement_file.close();
    }
  }

  double calculate_Af(int i) {
    const double xF = 2.5 / DELTA;
    const double x = DELTA * i;
    return (std::cos(50.0 * DELTA_T * i / (DELTA_T * NT)) * (x == xF));
  }

  void initialize_edge_condition() {
    U[0] = 0.0;
    V[0] = 0.0;
    U[NX - 1] = 0.0;
    V[NX - 1] = 0.0;
  }

  void initialize_start_condition() {
    for (int i = 0; i < NX; i++) {
      double x = i * DELTA;
      U[i] = std::exp(-(std::pow((x - X_A), 2) / (2 * std::pow(SIGMA, 2))));
      V[i] = 0.0;
    }
  }

  void calculate_a() {
    for (int i = 1; i < NX - 1; i++) {
      a[i] = (U[i + 1] - 2 * U[i] + U[i - 1]) / (DELTA * DELTA) -
             BETA * (U[i] - U_prev[i]) / DELTA_T + ALFA * calculate_Af(i);
    }
  }

  double calculate_E() {
    const double E_n = DELTA / 4.0 *
                       (std::pow((U[1] - U[0]) / DELTA, 2) +
                        std::pow((U[NX - 1] - U[NX - 2]) / DELTA, 2));
    double E_sum = 0.0;

    for (int i = 1; i < NX - 1; i++) {
      E_sum +=
          std::pow(V[i], 2) + std::pow((U[i + 1] - U[i - 1]) / (2 * DELTA), 2);
    }

    return E_n + (DELTA / 2.0 * E_sum);
  }

  void solve() {
    initialize_edge_condition();
    if (ALFA != 1.0) {
      initialize_start_condition();
    }
    U_prev = U;

    calculate_a();

    for (int t = 1; t <= NT; t++) {
      for (int i = 1; i < NX - 1; i++) {
        Vp[i] = V[i] + DELTA_T / 2.0 * a[i];
      }

      U_prev = U;

      for (int i = 1; i < NX - 1; i++) {
        U[i] = U[i] + DELTA_T * Vp[i];
      }
      calculate_a();

      for (int i = 1; i < NX - 1; i++) {
        V[i] = Vp[i] + DELTA_T / 2.0 * a[i];
      }

      double E = calculate_E();
      energy_file << t * DELTA_T << "\t" << E << "\n";

      for (int i = 0; i < U.size(); i++) {
        displacement_file << t * DELTA_T << "\t" << i * DELTA << "\t" << U[i]
                          << "\n";
      }

      displacement_file << "\n";
    }
  }
};

int main() {
  try {
    Varlete_Equation ve1(0.0, 0.0);
    ve1.solve();

    Varlete_Equation ve2(0.0, 0.1);
    ve2.solve();

    Varlete_Equation ve3(0.0, 1.0);
    ve3.solve();

    Varlete_Equation ve4(1.0, 1.0);
    ve4.solve();

    return 0;
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
}
