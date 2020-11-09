#include <branch_and_bound_method.hpp>
#include <iostream>
#include <vector>

int main() {
  std::vector<std::vector<double>> A = {{2, 1, 2}, {1, 2, 0}, {0, 0.5, 1}};
  std::vector<double> b = {6, 6, 2}, c = {2, 5, 3};
  Simplix_table table{A, b, c};
  auto answer = branch_and_bound_method(table);
  std::cout << "F=" << answer[0] << std::endl;
  for (size_t i = 1; i < answer.size(); i++) {
    std::cout << "x" << i << "=" << answer[i] << " ";
  }
  std::cout<<std::endl;
  brute_force_method(A, b, c);
}
