/*
 * File Name: dmrg.cpp
 * Description: ground state of effective model of brownian.
 * Created by Hao-Xin on 2023/07/12.
 *
 */

#include "brownian_type.h"
#include "params_case.h"
#include "myutil.h"
#include "brownian_type.h"
#include "params_case.h"
#include "myutil.h"

int main(int argc, char* argv[]) {
  // Check the command-line argument
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <filename>" << std::endl;
    return 1;
  }

  // Read parameters from the file
  CaseParams params(argv[1]);
  InitialLocalHilbertSpace();

  // Initialize operators
  std::vector<Tensor> sz_list(4, Tensor({pb_in, pb_out}));
  std::vector<Tensor> sp_list = sz_list;
  std::vector<Tensor> sm_list = sz_list;
  Tensor id = Tensor({pb_in, pb_out});

  for (size_t a = 0; a < 4; a++) {
    Tensor& sz = sz_list[a];
    for (size_t dim = 0; dim < local_hilbert_dim; dim++) {
      auto binary_form = BinaryTransform(dim, 4);
      int up = binary_form[a];
      sz({dim, dim}) = up ? 1.0 : -1.0;
    }
    assert(Div(sz) == qn0);
  }

  for (size_t i = 0; i < local_hilbert_dim; i++) {
    id({i, i}) = 1.0;
  }

  for (size_t i = 0; i < 4; i++) {
    Tensor& sp = sp_list[i];
    Tensor& sm = sm_list[i];
    for (size_t dim1 = 0; dim1 < local_hilbert_dim; dim1++) {
      auto binary_form = BinaryTransform(dim1, 4);
      if (binary_form[i] == 0) {
        size_t dim2 = dim1 + (1 << i);
        sm({dim1, dim2}) = 2.0;
        sp({dim2, dim1}) = 2.0;
      }
    }
  }

  // Calculate S^2 operator
  Tensor S2op;
  Tensor tmp1, tmp2, tmp3, tmp4;
  Contract(&sz_list[0], &sz_list[1], {{1}, {0}}, &tmp1);
  Contract(&sz_list[2], &sz_list[3], {{1}, {0}}, &tmp2);
  Contract(&sz_list[0], &sz_list[3], {{1}, {0}}, &tmp3);
  Contract(&sz_list[1], &sz_list[2], {{1}, {0}}, &tmp4);
  tmp1 += tmp2;
  tmp1 += (-tmp3);
  tmp1 += (-tmp4);
  S2op = tmp1;

  Tensor Sp_prod = Tensor({pb_in, pb_out});
  Tensor Sm_prod = Tensor({pb_in, pb_out});
  Sp_prod({0b1111, 0b0000}) = 16.0;
  Sm_prod({0b0000, 0b1111}) = 16.0;


  // Obtain the system size
  size_t L = params.L;
  SiteVecT sites(L, pb_out);

  // Create and initialize the MPS
  FiniteMPST mps(sites);

  // Check if the MPS file exists
  if (gqmps2::IsPathExist(gqmps2::kMpsPath) && L == GetNumofMps()) {
    // Do nothing, use the existing MPS file
  } else {
    std::cout << "No correct MPS File." << std::endl;
    std::cout << "Exit" << std::endl;
    return 1;
  }

  // Set the number of threads for tensor transpose and manipulation
  gqten::hp_numeric::SetTensorTransposeNumThreads(params.Threads);
  gqten::hp_numeric::SetTensorManipulationThreads(params.Threads);

  std::vector<std::vector<size_t>> sites_set;
  for (size_t i = L / 4 + 1; i < std::min(3 * L / 4 + 2, L); i++) {
    sites_set.push_back({L / 4, i});
  }

  mps.Load();
  // Measure the S^2 operator
  MeasureTwoSiteOp(mps, {S2op, S2op}, id, sites_set, "S2sym");
  MeasureTwoSiteOp(mps, {Sp_prod, Sm_prod}, id, sites_set, "PM");
  MeasureTwoSiteOp(mps, {Sm_prod, Sp_prod}, id, sites_set, "MP");

  // Measure the spin operators
//  for (size_t i = 0; i < 4; i++) {
//    MeasureTwoSiteOp(mps, {sp_list[i], sm_list[i]}, id, sites_set, "pm" + std::to_string(i));
//    MeasureTwoSiteOp(mps, {sm_list[i], sp_list[i]}, id, sites_set, "mp" + std::to_string(i));
//  }

  return 0;
}
