/*
 * File Name: dmrg.cpp
 * Description: ground state of 1D XXZ model
 * Created by Hao-Xin on 2023/05/08.
 *
 */

#include "hei_type.h"
#include "params_case.h"
#include "myutil.h"

int main(int argc, char *argv[]) {
  CaseParams params(argv[1]);

  Tensor sz({pb_in, pb_out});
  Tensor sp({pb_in, pb_out});
  Tensor sm({pb_in, pb_out});
  sz(0, 0) = 0.5;
  sz(1, 1) = -0.5;
  sp(0, 1) = 1;
  sm(1, 0) = 1;
  size_t L = 16;
  size_t Lx = 4;
  size_t Ly = 4;
  SiteVecT sites(L, pb_out);
  gqmps2::MPOGenerator<TenElemT, QNT> mpo_gen(sites, qn0);
  for (size_t i = 0; i < L - 1; ++i) {
    size_t x = i / Ly;
    size_t y = i % Ly;
    if (y != Ly - 1) {
      mpo_gen.AddTerm(1.0 * params.Delta, sz, i, sz, i + 1);
      mpo_gen.AddTerm(0.5, sp, i, sm, i + 1);
      mpo_gen.AddTerm(0.5, sm, i, sp, i + 1);
    }

    if (x != Lx - 1) {
      mpo_gen.AddTerm(1.0 * params.Delta, sz, i, sz, i + 4);
      mpo_gen.AddTerm(0.5, sp, i, sm, i + 4);
      mpo_gen.AddTerm(0.5, sm, i, sp, i + 4);
    }
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
