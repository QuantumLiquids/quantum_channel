/*
 * File Name: quench.cpp
 * Description: find the entanglement entropy of state exp(-beta H) |GS> where |GS> is ground state of XXZ model and H = \sum_i (-1)^i \sigma_z^i.
 * Created by Hao-Xin on 2023/05/19.
 */
#define GQTEN_COUNT_FLOPS 1

#include "hei_type.h"
#include "params_case.h"
#include "myutil.h"

void ActProjectors(
    FiniteMPST &mps,
    Tensor &p_a,
    Tensor &p_b
) {
  size_t L = mps.size();
  for (size_t i = 0; i < L; i = i + 2) {
    Tensor temp;
    Contract(mps(i), &p_a, {{1}, {0}}, &temp);
    temp.Transpose({0, 2, 1});
    mps[i] = std::move(temp);
  }
  for (size_t i = 1; i < L; i = i + 2) {
    Tensor temp;
    Contract(mps(i), &p_b, {{1}, {0}}, &temp);
    temp.Transpose({0, 2, 1});
    mps[i] = std::move(temp);
  }
}

int main(int argc, char *argv[]) {
  double beta(atof(argv[1]));
  std::cout << "beta=" << beta << std::endl;
  Tensor sz({pb_in, pb_out});
  Tensor sp({pb_in, pb_out});
  Tensor sm({pb_in, pb_out});
  Tensor projector_a({pb_in, pb_out}); //even site;
  Tensor projector_b({pb_in, pb_out}); //odd_site;
  sz(0, 0) = 0.5;
  sz(1, 1) = -0.5;
  sp(0, 1) = 1;
  sm(1, 0) = 1;
  projector_a(0, 0) = std::exp(-beta);
  projector_a(1, 1) = std::exp(beta);
  projector_b(0, 0) = std::exp(beta);
  projector_b(1, 1) = std::exp(-beta);
  size_t L = GetNumofMps();
  SiteVecT sites(L, pb_out);

  if (L % 2 == 1) {
    std::cout << "we now assume L is even!" << std::endl;
    return 0;
  }

  FiniteMPST mps(sites);
  mps.Load();

  ActProjectors(mps, projector_a, projector_b);
  std::vector<double> ee_list = mps.GetEntanglementEntropy(1);

  std::copy(ee_list.begin(), ee_list.end(), std::ostream_iterator<double>(std::cout, " "));

  std::cout << std::endl;
  return 0;
}