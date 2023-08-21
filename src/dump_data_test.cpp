//
// Created by Hao-Xin on 2023/6/16.
//

#include <iostream>
#include <fstream>
#include <cstdlib>

int main() {
  std::string file_name = "test";
  std::ofstream ofs(file_name, std::ofstream::binary);
  double a = 1.02;
  ofs << a;
//  ofs.write((const char *) renyi2_a.data(), renyi2_a.size() * sizeof(double));
//  ofs.write((const char *) renyi2_b.data(), renyi2_b.size() * sizeof(double));
  ofs << std::endl;
  ofs.close();
}