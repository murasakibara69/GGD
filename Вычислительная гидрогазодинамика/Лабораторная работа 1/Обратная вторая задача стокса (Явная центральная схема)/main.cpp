#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

int main(void) {
  std::size_t n_y;
  double time, nu, a, omega;
  std::ifstream input;
  input.exceptions(std::ifstream::badbit | std::ifstream::failbit);
  try {
    input.open("input.dat");
  } catch (const std::ifstream::failure &ex) {
    std::cerr << ex.what() << std::endl;
    std::cerr << ex.code() << std::endl;
  }
  input >> n_y;
  input >> time;
  input >> nu;
  input >> a >> omega;
  input.close();
  std::cout << "N = " << n_y << std::endl;
  std::cout << "Time = " << time << std::endl;
  std::cout << "nu = " << nu << std::endl;
  std::cout << "A = " << a << ", omega = " << omega << std::endl;

  const double delta = std::sqrt(2.0 * nu / omega);
  const double delta_y = 2.0 * delta / (n_y - 1);
  const double delta_t = std::pow(delta_y, 2.0) / (6.0 * nu);
  const std::size_t n_t = time / delta_t + 1;

  const size_t j = n_y / 2 + 1;
  std::vector<double> u;
  for (std::size_t i = 0; i < n_y; ++i)
    u.push_back(0.0);

  double t = delta_t;
  std::vector<double> u_next(n_y);
  for (std::size_t n = 1; n < n_t; ++n) {
    for (std::size_t i = 1; i < n_y - 1; ++i)
      u_next.at(i) = u.at(i) +
                     nu * delta_t / (2.0 * std::pow(delta_y, 2.0)) *
                         (u.at(i - 1) - 2.0 * u.at(i) + u.at(i + 1)) -
                     a * std::cos(omega * t) * delta_t;
    u_next.front() = 0.0;
    for (std::size_t i = j + 1; i < n_y; ++i)
      u_next.at(i) = u_next.at(j);
    u = u_next;
    t += delta_t;
  }

  std::vector<double> y;
  for (std::size_t i = 0; i < n_y; ++i)
    y.push_back(delta_y * i);
  std::ofstream output;
  output.exceptions(std::ofstream::badbit | std::ofstream::failbit);
  try {
    output.open("result.dat");
  } catch (const std::ofstream::failure &ex) {
    std::cerr << ex.what() << std::endl;
    std::cerr << ex.code() << std::endl;
  }
  for (std::size_t i = 0; i < n_y; ++i)
    output << u.at(i) << '\t' << y.at(i) << std::endl;
  output.close();

  std::system("gnuplot -p plot.plt");

  return EXIT_SUCCESS;
}
