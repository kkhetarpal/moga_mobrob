#ifndef AWI_H
#define AWI_H

#include <sstream>
#include <vector>
#include <iostream>

#define SHOW(x) std::cout << #x << " :: " << x << std::endl

template <typename t>
std::string num2str(t const num){
  std::stringstream ss;
  ss << num;
  return ss.str();
}

template<typename T>
void show(std::vector<std::vector<T> > a) {
  for(size_t i = 0; i < a.size(); ++i) {
    std::cout << i << "\t";
    for(size_t j = 0; j < a[i].size(); ++j) {
      std::cout << a[i][j] << "\t";
    }
    std::cout << "\n";
  }
}

template<typename T>
void show(T * base, size_t const n) {
  for(size_t i = 0; i < n; ++i) {
    std::cout << base[i] << "\t";
  }
  std::cout << std::endl;
}

template<typename T>
void show(std::vector<T> const & v) {
  for(size_t i = 0; i < v.size(); ++i) {
    std::cout << v[i] << "\t";
  }
  std::cout << std::endl;
}

template<typename T>
void show(T * base, size_t const nr, size_t const nc) {
  for(size_t i = 0; i < nr; ++i) {
    for(size_t j = 0; j < nc; ++j) {
      std::cout << base[i*nc+j] << "\t";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

#endif
