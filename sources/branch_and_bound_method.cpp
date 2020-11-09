
#include <branch_and_bound_method.hpp>

Simplix_table::Simplix_table(std::vector<std::vector<double>> const &A,
                             std::vector<double> const &b,
                             std::vector<double> const &c) {
  for (size_t i = 1; i <= c.size(); i++) {
    _free_variables.push_back(i);
    _basic_variables.push_back(i + c.size());
  }
  std::vector<std::vector<double>> table(c.size() + 1,
                                         std::vector<double>(c.size() + 1));
  for (size_t i = 0; i < b.size(); i++) {
    table[i][0] = b[i];
  }
  table[c.size()][0] = 0;
  for (size_t i = 0; i < c.size(); i++) {
    for (size_t j = 0; j < c.size(); j++) {
      table[i][j + 1] = A[i][j];
    }
  }
  for (size_t i = 0; i < b.size(); i++) {
    table[b.size()][i + 1] = -c[i];
  }
  _table = table;
}

std::vector<size_t> Simplix_table::basic_variables() {
  return _basic_variables;
}

std::vector<size_t> Simplix_table::free_variables() { return _free_variables; }

std::vector<std::vector<double>> Simplix_table::table() { return _table; }

void Simplix_table::transformation(size_t string, size_t colunm) {
  size_t variable = _basic_variables[string];
  _basic_variables[string] = _free_variables[colunm - 1];
  _free_variables[colunm - 1] = variable;
  for (size_t i = 0; i < _table.size(); i++) {
    for (size_t j = 0; j < _table.at(0).size(); j++) {
      if (i != string && j != colunm) {
        _table.at(i).at(j) = _table.at(i).at(j) -
                             _table.at(i).at(colunm) * _table.at(string).at(j) /
                                 _table.at(string).at(colunm);
      }
    }
  }
  for (size_t i = 0; i < _table.size(); i++) {
    if (i != string) {
      _table.at(i).at(colunm) =
          -_table.at(i).at(colunm) / _table.at(string).at(colunm);
    }
  }
  for (size_t i = 0; i < _table.at(0).size(); i++) {
    if (i != colunm) {
      _table.at(string).at(i) =
          _table.at(string).at(i) / _table.at(string).at(colunm);
    }
  }
  _table[string][colunm] = 1 / _table[string][colunm];
}

size_t Simplix_table::permissive_string(size_t permissive_column) {
  size_t permissive_string;
  double min;
  bool check = true;
  double simplex_relation;
  for (size_t i = 0; i < table().size() - 1; i++) {
    if (table()[i][permissive_column] != 0) {
      simplex_relation = table()[i][0] / table()[i][permissive_column];
      if (simplex_relation >= 0 && (simplex_relation < min || check)) {
        min = simplex_relation;
        check = false;
        permissive_string = i;
      }
    }
  }
  return permissive_string;
}

void Simplix_table::print_table() {
  std::cout << std::setw(5) << "     " << std::setw(8) << 'x' << 0;
  for (size_t i = 0; i < _free_variables.size(); i++) {
    std::cout << std::setw(7) << 'x' << _free_variables.at(i);
  }
  std::cout << std::endl;
  for (size_t i = 0; i < table().size(); i++) {
    if (i != _table.size() - 1)
      std::cout << std::setw(5) << 'x' << _basic_variables.at(i);
    else
      std::cout << "     F";
    for (size_t j = 0; j < _table.at(i).size(); j++) {
      std::cout << std::setw(8) << _table.at(i).at(j);
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

std::vector<double> Simplix_table::answer() {
  std::vector<double> answer(_free_variables.size() + _basic_variables.size() +
                             1);
  answer[0] = table().at(table().size() - 1).at(0);
  for (size_t i = 0; i < table().size() - 1; i++) {
    answer.at(_basic_variables.at(i)) = table().at(i).at(0);
  }
  for (size_t i = 0; i < table().at(0).size() - 1; i++) {
    answer.at(_free_variables.at(i)) = 0;
  }
  return answer;
}

void Simplix_table::empty() {
  _basic_variables.clear();
  _free_variables.clear();
  _table.clear();
}

bool Simplix_table::is_empty() {
  return _table.empty() && _free_variables.empty() && _basic_variables.empty();
}

void Simplix_table::simplix_method() {
  bool cheak;

  do {
    size_t column = 0;
    cheak = false;
    for (size_t i = 0; i < _basic_variables.size(); i++) {
      if (table().at(i).at(0) < 0) {
        for (size_t j = 1; j <= _free_variables.size(); j++) {
          if (table().at(i).at(j) < 0) {
            column = j;
            size_t string = permissive_string(column);
            transformation(string, column);
            cheak = true;
            break;
          }
          if (j == _free_variables.size()) {
            print_table();
            empty();
            return;
          }
        }
      }
    }
  } while (cheak);
  double min = 0;
  do {
    min = 0;
    size_t permissive_column = 0;
    for (size_t i = 1; i <= _free_variables.size(); i++) {
      if (table().at(_basic_variables.size()).at(i) < 0 &&
          table().at(_basic_variables.size()).at(i) < min) {
        min = table().at(_basic_variables.size()).at(i);
        permissive_column = i;
      }
    }
    if (min == 0)
      break;
    size_t perm_string = permissive_string(permissive_column);
    transformation(perm_string, permissive_column);

  } while (min < 0);

  print_table();
}

void Simplix_table::add_constraint(size_t number_x, char sign, double number) {
  _basic_variables.push_back(_basic_variables.size() + _free_variables.size() +
                             1);
  if (sign == '>') {
    _table.insert(_table.end() - 1, _table.at(number_x));
    _table.at(_table.size() - 2).at(0) -= number;
  } else if (sign == '<') {
    _table.insert(_table.end() - 1, _table.at(number_x));
    for (auto &value : _table.at(_table.size() - 2)) {
      value *= -1;
    }
    _table.at(_table.size() - 2).at(0) += number;
  }
}

bool is_integer(double value) {
  auto int_value = static_cast<int>(value);
  return value == int_value;
}

std::vector<double> branch_and_bound_method(Simplix_table simplix_table) {
  std::cout << "initial simplix table" << std::endl;
  simplix_table.print_table();
  simplix_table.simplix_method();
  if (simplix_table.is_empty()) {
    return std::vector<double>{0};
  }
  auto answer = simplix_table.answer();
  for (size_t i = 0; i < simplix_table.basic_variables().size(); i++) {
    if (!is_integer(simplix_table.table().at(i).at(0)) &&
        simplix_table.basic_variables().at(i) <=
            simplix_table.free_variables().size()) {
      Simplix_table simplix_table1 = simplix_table;
      Simplix_table simplix_table2 = simplix_table;
      simplix_table1.add_constraint(i, '<',
                                    static_cast<int>(static_cast<int>(
                                        simplix_table.table().at(i).at(0))));
      simplix_table2.add_constraint(i, '>',
                                    static_cast<int>(static_cast<int>(
                                        simplix_table.table().at(i).at(0))) +
                                        1);
      std::cout << "X <= "
                << static_cast<int>(
                       static_cast<int>(simplix_table.table().at(i).at(0)))
                << " for X" << simplix_table.basic_variables().at(i)
                << std::endl;
      auto answer1 = branch_and_bound_method(simplix_table1);
      std::cout << "X>= "
                << static_cast<int>(
                       static_cast<int>(simplix_table.table().at(i).at(0))) +
                       1
                << " for X" << simplix_table.basic_variables().at(i)
                << std::endl;
      auto answer2 = branch_and_bound_method(simplix_table2);
      if (answer1.size() == 1 && answer2.size() == 1) {
        return std::vector<double>{0};
      } else if (answer1.size() == 1) {
        return answer2;
      } else if (answer2.size() == 1) {
        return answer1;
      } else {
        if (answer1.at(0) > answer2.at(0)) {
          return answer1;
        } else {
          return answer2;
        }
      }
    }
  }
  return answer;
}

bool check_constraints(const int &i, const int &j, const int &k, const std::vector<std::vector<double>> &A,
                       const std::vector<double> &b) {
  return A.at(0).at(0) * i + A.at(0).at(1) * j + A.at(0).at(2) * k <= b.at(0) &&
         A.at(1).at(0) * i + A.at(1).at(1) * j + A.at(1).at(2) * k <= b.at(1) &&
         A.at(2).at(0) * i + A.at(2).at(1) * j + A.at(2).at(2) * k <= b.at(2);
}

void brute_force_method(const std::vector<std::vector<double>>& A,
                        const std::vector<double>& b, const std::vector<double> &c) {
  for (auto i = 0, j = 0, k = 0;;) {
    if (check_constraints(i, j, k, A, b)) {
      std::cout << "F(" << i << ", " << j << ", " << k
                << ")=" << c.at(0) * i + c.at(1) * j + c.at(2) * k << std::endl;
      k++;
    } else {
      if (k > 0) {
        k = 0;
        j++;
      } else if (k == 0 && j > 0) {
        j = 0;
        i++;
      } else {
        break;
      }
    }
  }
}