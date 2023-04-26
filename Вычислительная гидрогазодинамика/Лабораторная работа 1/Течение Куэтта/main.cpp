#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

void sweep_method(const std::vector<double> &a, const std::vector<double> &b,
                  const std::vector<double> &c, const std::vector<double> &d,
                  std::vector<double> &x);

int main(void) {
  std::size_t n_y;
  double h, time, nu, a, u_0;
  std::ifstream input;
  input.exceptions(std::ifstream::badbit | std::ifstream::failbit);
  try {
    input.open("input.dat");
  } catch (const std::ifstream::failure &ex) {
    std::cerr << ex.what() << std::endl;
    std::cerr << ex.code() << std::endl;
  }
  input >> n_y;
  input >> h;
  input >> time;
  input >> nu;
  input >> a;
  input >> u_0;
  input.close();
  std::cout << "N = " << n_y << std::endl;
  std::cout << "H = " << h << std::endl;
  std::cout << "Time = " << time << std::endl;
  std::cout << "nu = " << nu << std::endl;
  std::cout << "A = " << a << std::endl;
  std::cout << "U_0 = " << u_0 << std::endl;

  const double delta_y = h / (n_y - 1);
  const double delta_t = std::pow(delta_y, 2.0) / 2.0;
  const std::size_t n_t = time / delta_t + 1;

  std::vector<double> u;
  for (std::size_t i = 0; i < n_y; ++i)
    u.push_back(0.0);
  u.back() = u_0;

  const double C_A = nu / std::pow(delta_y, 2.0);
  const double C_D = -2.0 / delta_t;
  const double C_B = C_D - 2.0 * C_A;

  std::vector<double> A(n_y - 2, C_A), B(n_y - 2, C_B), C = A, D(n_y - 2);

  std::vector<double> u_next(n_y);
  for (std::size_t n = 1; n < n_t; ++n) {
    D.front() = C_D * u.at(1) - A.front() * u.front() - a;
    for (std::size_t i = 1; i < D.size() - 1; ++i)
      D.at(i) = C_D * u.at(i + 1) - a;
    D.back() = C_D * u.at(D.size()) - C.back() * u.back() - a;
    std::vector<double> x(u.begin() + 1, u.end() - 1);
    sweep_method(A, B, C, D, x);
    x.push_back(u.back());
    x.insert(x.begin(), u.front());
    u = x;
    for (std::size_t i = 1; i < n_y - 1; ++i)
      u_next.at(i) = u.at(i) +
                     nu * delta_t / (2.0 * std::pow(delta_y, 2.0)) *
                         (u.at(i - 1) - 2.0 * u.at(i) + u.at(i + 1)) +
                     2.0 * a / delta_t;
    u_next.front() = 0.0;
    u_next.back() = u_0;
    u = u_next;
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

void sweep_method(const std::vector<double> &a, const std::vector<double> &b,
                  const std::vector<double> &c, const std::vector<double> &d,
                  std::vector<double> &x) {
  std::vector<double> alpha, beta, gamma;
  gamma.push_back(b.front());
  alpha.push_back(-c.front() / gamma.front());
  beta.push_back(d.front() / gamma.front());
  for (std::size_t i = 1; i < a.size() - 1; ++i) {
    gamma.push_back(b.at(i) + a.at(i) * alpha.back());
    alpha.push_back(-c.at(i) / gamma.back());
    beta.push_back((d.at(i) - a.at(i) * beta.back()) / gamma.back());
  }
  gamma.push_back(b.back() + a.back() * alpha.back());
  beta.push_back((d.back() - a.back() * beta.back()) / gamma.back());
  x.back() = beta.back();
  for (std::size_t i = x.size() - 2; i < SIZE_MAX; --i) {
    x.at(i) = alpha.at(i) * x.at(i + 1) + beta.at(i);
  }
}
