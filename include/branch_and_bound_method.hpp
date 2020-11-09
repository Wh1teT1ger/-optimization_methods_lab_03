
#ifndef SIMPLEX_METHOD_SIMPLEX_METHOD_HPP
#define SIMPLEX_METHOD_SIMPLEX_METHOD_HPP

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>

class Simplix_table {
  std::vector<size_t> _basic_variables;
  std::vector<size_t> _free_variables;
  std::vector<std::vector<double>> _table;

public:
  Simplix_table(std::vector<std::vector<double>> const &A,
                std::vector<double> const &b, std::vector<double> const &c);

  std::vector<size_t> basic_variables();

  void empty();

  bool is_empty();

  std::vector<size_t> free_variables();

  std::vector<std::vector<double>> table();

  void transformation(size_t string, size_t colunm);

  void print_table();

  size_t permissive_string(size_t permissive_column);

  std::vector<double> answer();

  void simplix_method();

  void add_constraint(size_t number_x, char sign, double number);


};

bool is_integer(double value);

std::vector<double> branch_and_bound_method(Simplix_table simplix_table);

void brute_force_method(const std::vector<std::vector<double>>& A,
                        const std::vector<double>& b, const std::vector<double> &c);

bool check_constraints(const int &i, const int &j, const int &k, const std::vector<std::vector<double>> &A,
                       const std::vector<double> &b);

#endif // SIMPLEX_METHOD_SIMPLEX_METHOD_HPP
