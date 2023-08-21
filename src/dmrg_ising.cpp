/*
 * File Name: dmrg_ising.cpp
 * Description: ground state of 1D Ising model
 * Created by Hao-Xin on 2023/05/08.
 *
 */

#include "ising_type.h"
#include "params_case.h"
#include "myutil.h"

int main(int argc, char *argv[]) {
  CaseParams params(argv[1]);

  Tensor sigma_z({pb_in, pb_out});
  Tensor sigma_x({pb_in, pb_out});
  Tensor sp({pb_in, pb_out});
  Tensor sm({pb_in, pb_out});
  sigma_z(0, 0) = 1;
  sigma_z(1, 1) = -1;
  sp(0, 1) = 1;
  sm(1, 0) = 1;
  sigma_x(0, 1) = 1;
  sigma_x(1, 0) = 1;
  size_t L = params.L;
  SiteVecT sites(L, pb_out);
  gqmps2::MPOGenerator<TenElemT, QNT> mpo_gen(sites, qn0);
  for (size_t i = 0; i < L - 1; ++i) {
//    mpo_gen.AddTerm(1.0, sp, i, sm, i + 1);
//    mpo_gen.AddTerm(1.0, sm, i, sp, i + 1);
    mpo_gen.AddTerm(1.0, sigma_x, i, sigma_x, i + 1);
  }
  for (size_t i = 0; i < L; ++i) {
    mpo_gen.AddTerm(1.0, sigma_z, i);
  }
  auto mpo = mpo_gen.Gen();

  gqmps2::SweepParams sweep_params(
      params.Sweeps,
      params.Dmin, params.Dmax, params.CutOff,
      gqmps2::LanczosParams(params.LanczErr, params.MaxLanczIter)
  );

  FiniteMPST mps(sites);
  std::vector<size_t> stat_labs;
  for (size_t i = 0; i < L; ++i) { stat_labs.push_back(i % 2); }

  if (gqmps2::IsPathExist(gqmps2::kMpsPath)) {
    if (L == GetNumofMps()) {
      std::cout << "The number of mps files is consistent with mps size." << std::endl;
      std::cout << "Directly use mps from files." << std::endl;
    } else {
      gqmps2::DirectStateInitMps(mps, stat_labs);
      std::cout << "MPS files number incorrect. Initial mps as direct product state." << std::endl;
      mps.Dump(sweep_params.mps_path, true);
    }
  } else {
    gqmps2::DirectStateInitMps(mps, stat_labs);
    std::cout << "No MPS File. Initial mps as direct product state." << std::endl;
    mps.Dump(sweep_params.mps_path, true);
  }

  gqten::hp_numeric::SetTensorTransposeNumThreads(params.Threads);
  gqten::hp_numeric::SetTensorManipulationThreads(params.Threads);

  auto e0 = gqmps2::TwoSiteFiniteVMPS(mps, mpo, sweep_params);
  std::cout << "E0/site: " << e0 / L << std::endl;    // E0/site: -0.43412

  mps.Load(sweep_params.mps_path);
  auto ee2_list = mps.GetEntanglementEntropy(2);

  std::copy(ee2_list.begin(), ee2_list.end(), std::ostream_iterator<double>(std::cout, " "));

  return 0;
}
