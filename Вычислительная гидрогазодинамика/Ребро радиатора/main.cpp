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

  const double C_A = a / std::pow(delta_x, 2.0);
  const double C_D = -2.0 / delta_t;
  const double C_B = C_D - 2.0 * C_A;

  std::vector<double> A(n_x - 2, C_A), B(n_x - 2, C_B), C(A.begin(), A.end()),
      D(n_x - 2);

  std::vector<double> T_next(n_x);
  for (std::size_t n = 1; n < n_t; ++n) {
    D.front() = -A.front() * T.front() + C_D * T.at(1) -
                2.0 * std::sqrt(M_PI / area) * alpha * (T_e - T.at(1));
    for (std::size_t i = 1; i < D.size() - 1; ++i)
      D.at(i) = C_D * T.at(i + 1) -
                2.0 * std::sqrt(M_PI / area) * alpha * (T_e - T.at(i + 1));
    D.back() = -C.back() * T.back() + C_D * T.at(D.size()) -
               2.0 * std::sqrt(M_PI / area) * alpha * (T_e - T.at(D.size()));
    std::vector<double> y(T.begin() + 1, T.end() - 1);
    sweep_method(A, B, C, D, y);
    y.push_back(T.back());
    y.insert(y.cbegin(), T.front());
    T = y;
    for (std::size_t i = 1; i < n_x - 1; ++i)
      T_next.at(i) =
          T.at(i) +
          a * delta_t * (T.at(i - 1) - 2.0 * T.at(i) + T.at(i + 1)) /
              (2.0 * std::pow(delta_x, 2.0)) +
          4.0 * alpha * delta_t * std::sqrt(M_PI / area) * (T_e - T.at(i));
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
