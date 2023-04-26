#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

int main(void) {
  std::size_t n_x;
  double time, l, area, c, rho, T_0, T_e, alpha, lambda, q_L;
  std::ifstream input;
  input.exceptions(std::ifstream::badbit | std::ifstream::failbit);
  try {
    input.open("input.dat");
  } catch (const std::ifstream::failure &ex) {
    std::cerr << ex.what() << std::endl;
    std::cerr << ex.code() << std::endl;
  }
  input >> n_x;
  input >> time;
  input >> l >> area;
  input >> c >> rho;
  input >> T_0 >> T_e;
  input >> alpha >> lambda;
  input >> q_L;
  input.close();
  std::cout << "N = " << n_x << std::endl;
  std::cout << "Time = " << time << std::endl;
  std::cout << "L = " << l << ", A = " << area << std::endl;
  std::cout << "c = " << c << ", rho = " << rho << std::endl;
  std::cout << "T_0 = " << T_0 << ", T_e = " << T_e << std::endl;
  std::cout << "alpha = " << alpha << ", lambda = " << lambda << std::endl;
  std::cout << "q_L = " << q_L << std::endl;

  const double a = lambda / (rho * c);
  const double delta_x = l / (n_x - 1);
  const double delta_t = std::pow(delta_x, 2.0) / (6.0 * a);
  const std::size_t n_t = time / delta_t + 1;

  std::vector<double> T(n_x);
  T.front() = T.at(1) + q_L * delta_x / lambda;
  T.back() = (alpha * T_e / lambda + T.at(n_x - 2) / delta_x) /
             (alpha / lambda + 1.0 / delta_x);

  std::vector<double> T_next(n_x);
  for (std::size_t n = 1; n < n_t; ++n) {
    for (std::size_t i = 1; i < n_x - 1; ++i)
      T_next.at(i) =
          T.at(i) +
          a * delta_t * (T.at(i - 1) - 2.0 * T.at(i) + T.at(i + 1)) /
              (2.0 * std::pow(delta_x, 2.0)) +
          2.0 * alpha * delta_t * std::sqrt(M_PI / area) * (T_e - T.at(i));
    T_next.front() = T_next.at(1) + q_L * delta_x / lambda;
    T_next.back() = (alpha * T_e / lambda + T_next.at(n_x - 2) / delta_x) /
                    (alpha / lambda + 1.0 / delta_x);
    T = T_next;
  }

  std::vector<double> x;
  for (std::size_t i = 0; i < n_x; ++i)
    x.push_back(delta_x * i);
  std::ofstream output;
  output.exceptions(std::ofstream::badbit | std::ofstream::failbit);
  try {
    output.open("result.dat");
  } catch (const std::ofstream::failure &ex) {
    std::cerr << ex.what() << std::endl;
    std::cerr << ex.code() << std::endl;
  }
  for (std::size_t i = 0; i < n_x; ++i)
    output << x.at(i) << '\t' << T.at(i) << std::endl;
  output.close();

  std::system("gnuplot -p plot.plt");

  return EXIT_SUCCESS;
}
