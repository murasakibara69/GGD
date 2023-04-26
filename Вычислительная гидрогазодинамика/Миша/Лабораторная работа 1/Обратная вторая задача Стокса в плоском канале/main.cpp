#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <vector>

int main(void) {
  std::size_t n_y;
  double h, time, nu, a, omega;
  std::ifstream in;
  in.open("input.dat");
  if (in.is_open()) {
    in >> n_y;
    in >> h;
    in >> time;
    in >> nu;
    in >> a >> omega;
  }
  in.close();

  double delta_y = h / (n_y - 1);
  double delta_t = std::pow(delta_y, 2.0) / (6.0 * nu);
  std::size_t n_t = time / delta_t + 1;

  std::vector<double> u(n_y);
  std::vector<double> t;
  for (std::size_t n = 0; n < n_t; n++)
    t.push_back(n * delta_t);

  std::vector<double> u1(n_y), u2(n_y);
  for (std::size_t n = 1; n < n_t; n++) {
    u1.front() = 0.0;
    for (std::size_t i = 1; i < n_y - 1; i++)
      u1[i] = (-a * std::cos(omega * t[n]) +
               nu * (u1[i - 1] - u[i] + u[i + 1]) / std::pow(delta_y, 2.0) +
               u[i] / delta_t) /
              (1.0 / delta_t + nu / std::pow(delta_y, 2.0));
    u1.back() = 0.0;
    u2.front() = 0.0;
    for (std::size_t i = 0; i < n_y - 2; i++)
      u2[i + 1] = (-a * std::cos(omega * t[n]) - (u2[i] - u[i]) / delta_t) *
                      std::pow(delta_y, 2.0) / nu +
                  (u[i - 1] - u[i] - u2[i]);
    u2.back() = 0.0;
    for (std::size_t i = 0; i < n_y; i++)
      u[i] = (u1[i] + u[i]) / 2.0;
  }

  std::vector<double> u_next(n_y);
  for (std::size_t n = 1; n < n_t; n++) {
    u_next.front() = 0.0;
    for (std::size_t i = 1; i < n_y - 1; i++)
      u_next[i] = u[i] - a * std::cos(omega * t[n]) * delta_t +
                  nu * delta_t * (u[i - 1] - 2.0 * u[i] + u[i + 1]) /
                      std::pow(delta_y, 2.0);
    u_next.back() = 0.0;
    u = u_next;
  }
  u_next.clear();

  std::vector<double> y;
  for (std::size_t i = 0; i < n_y; i++)
    y.push_back(i * delta_y);
  std::ofstream out;
  out.open("result.dat");
  if (out.is_open()) {
    for (std::size_t i = 0; i < n_y; i++)
      out << u[i] << '\t' << y[i] << std::endl;
  }
  out.close();

  return 0;
}
