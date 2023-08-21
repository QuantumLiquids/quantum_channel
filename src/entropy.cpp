/*
 * File Name: entropy.cpp
 * Description: entanglement entropy, mutual information and so on for ground state of XXZ model with quantum channel.
 * Created by Hao-Xin on 2023/05/08.
 */
#define GQTEN_COUNT_FLOPS 1

#include "hei_type.h"
#include "params_case.h"
#include "myutil.h"
#include "calc_entropy.h"

int main(int argc, char *argv[]) {
  double lambda(atof(argv[1]));
  std::cout << "lambda=" << lambda << std::endl;
  std::string channel_type = argv[2];
  std::cout << "channel_type=" << channel_type << std::endl;
  std::size_t thread_num = atoi(argv[3]);
  std::cout << "thread_num = " << thread_num << std::endl;

  Tensor sz({pb_in, pb_out});
  Tensor sp({pb_in, pb_out});
  Tensor sm({pb_in, pb_out});
  sz(0, 0) = 0.5;
  sz(1, 1) = -0.5;
  sp(0, 1) = 1;
  sm(1, 0) = 1;
  size_t L = GetNumofMps();
  SiteVecT sites(L, pb_out);

  if (L % 2 == 1) {
    std::cout << "we now assume L is even!" << std::endl;
    return 0;
  }

  FiniteMPST mps(sites);
  mps.Load();

  gqten::hp_numeric::SetTensorTransposeNumThreads(thread_num);
  gqten::hp_numeric::SetTensorManipulationThreads(thread_num);

  double thermal_renyi2_entropy;
  std::vector<double> renyi2_a, renyi2_b;
  CalRenyi2(mps, lambda, thermal_renyi2_entropy, renyi2_a, renyi2_b, channel_type);

  std::string
      file_name = "renyi2_entropyL" + std::to_string(L) + "channel" + channel_type + "lambda" + std::to_string(lambda);
  std::ofstream ofs(file_name, std::ofstream::binary);
  ofs.write((const char *) (&thermal_renyi2_entropy), sizeof(double));
  ofs.write((const char *) renyi2_a.data(), renyi2_a.size() * sizeof(double));
  ofs.write((const char *) renyi2_b.data(), renyi2_b.size() * sizeof(double));
  ofs << std::endl;
  ofs.close();

//  std::cout << "thermal renyi2 entropy = " << thermal_renyi2_entropy << std::endl;


//  std::cout << "entanglement renyi2 entropy (subsystem A): " << std::endl;
//  std::copy(renyi2_a.begin(), renyi2_a.end(), std::ostream_iterator<double>(std::cout, " "));
//  std::cout << std::endl;

//  std::cout << "entanglement renyi2 entropy (subsystem B): " << std::endl;
//  std::copy(renyi2_b.begin(), renyi2_b.end(), std::ostream_iterator<double>(std::cout, " "));
//  std::cout << std::endl;


//  double thermal_renyi3_entropy, negativity3;
//  CalRenyi3(mps, lambda, thermal_renyi3_entropy, negativity3);
//  std::cout << "thermal renyi3 entropy = " << thermal_renyi3_entropy << std::endl;
//  std::cout << "negativity3 entropy = " << negativity3 << std::endl;
  return 0;
}