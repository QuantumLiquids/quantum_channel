/*
 * File Name: dmrg.cpp
 * Description: ground state of effective model of brownian.
 * Created by Hao-Xin on 2023/06/25.
 *
 */

//#define PLAIN_TRANSPOSE 1

#include "brownian_type.h"
#include "params_case.h"
#include "myutil.h"

int main(int argc, char *argv[]) {
  namespace mpi = boost::mpi;
  mpi::environment env(mpi::threading::multiple);
  if (env.thread_level() < mpi::threading::multiple) {
    std::cout << "thread level of env is not right." << std::endl;
    env.abort(-1);
  }
  mpi::communicator world;

  CaseParams params(argv[1]);
  InitialLocalHilbertSpace();

  std::vector<Tensor> sz_list(4, Tensor({pb_in, pb_out}));
  auto sp_list = sz_list, sm_list = sz_list;
  Tensor id = Tensor({pb_in, pb_out});
  for (size_t a = 0; a < 4; a++) {
    Tensor &sz = sz_list[a];
    for (size_t dim = 0; dim < local_hilbert_dim; dim++) {
      auto binary_form = BinaryTransform(dim, 4);
      int up = binary_form[a];
      if (up) {
        sz({dim, dim}) = 1.0;
      } else {
        sz({dim, dim}) = -1.0;
      }
    }
    assert(Div(sz) == qn0);
  }

  for (size_t i = 0; i < local_hilbert_dim; i++) {
    id({i, i}) = 1.0;
  }

  for (size_t i = 0; i < 4; i++) {
    Tensor &sp = sp_list[i];
    Tensor &sm = sm_list[i];
    for (size_t dim1 = 0; dim1 < local_hilbert_dim; dim1++) {
      auto binary_form = BinaryTransform(dim1, 4);
      if (binary_form[i] == 0) {
        size_t dim2 = dim1 + (1 << i);
        sm({dim1, dim2}) = 2.0;
        sp({dim2, dim1}) =
            2.0; // follow the matrix, don't care the whether operator act on state or conjugate of state...
      }
    }
  }

  size_t L = params.L;
  SiteVecT sites(L, pb_out);
  gqmps2::MPOGenerator<TenElemT, QNT> mpo_gen(sites, qn0);

  Tensor sz_alternating_sum = sz_list[0] + (-sz_list[1]) + sz_list[2] + (-sz_list[3]);
  Tensor sz_sum = sz_list[0] + sz_list[1] + sz_list[2] + sz_list[3];
  Tensor sz_alternating_sum_square;
  Contract(&sz_alternating_sum, &sz_alternating_sum, {{1}, {0}}, &sz_alternating_sum_square);

  for (size_t i = 0; i < L; ++i) {
    int sign_factor = ((int) i % 2) * 2 - 1;
    mpo_gen.AddTerm(params.gamma * double(sign_factor), sz_sum, i);
    mpo_gen.AddTerm(params.h / 2, sz_alternating_sum_square, i);
  }

  std::vector<std::pair<Tensor, Tensor>> spin_coupling_operator_pairs(12);
  std::vector<double> spin_coupling_coefficients(12);
  std::vector<std::pair<Tensor, Tensor>> spin_coupling_square_operator_pairs(144);
  std::vector<double> spin_coupling_square_coefficients(144);
  for (size_t a = 0; a < 4; a++) {
    spin_coupling_operator_pairs[3 * a] = std::make_pair(sp_list[a], sm_list[a]);
    spin_coupling_operator_pairs[3 * a + 1] = std::make_pair(sm_list[a], sp_list[a]);
    spin_coupling_operator_pairs[3 * a + 2] = std::make_pair(sz_list[a], sz_list[a]);

    spin_coupling_coefficients[3 * a] = -0.5 * (2 * ((int) a % 2) - 1);
    spin_coupling_coefficients[3 * a + 1] = -0.5 * (2 * ((int) a % 2) - 1);
    spin_coupling_coefficients[3 * a + 2] = -params.Delta * (2 * ((int) a % 2) - 1);
  }
  for (size_t i = 0; i < 12; i++) {
    for (size_t j = 0; j < 12; j++) {
      size_t k = 12 * i + j;
      std::pair<Tensor, Tensor> &op_pair1 = spin_coupling_operator_pairs[i];
      std::pair<Tensor, Tensor> &op_pair2 = spin_coupling_operator_pairs[j];
      Tensor product_op_a, product_op_b;
      if (i == j) {
        if (i % 3 == 2) { //sz * sz
          product_op_a = id;
          product_op_b = id;
        } else {
          product_op_a = Tensor({pb_in, pb_out});
          product_op_b = Tensor({pb_in, pb_out});
        }
      } else {
        Contract(&op_pair1.first, &op_pair2.first, {{1}, {0}}, &product_op_a);
        Contract(&op_pair1.second, &op_pair2.second, {{1}, {0}}, &product_op_b);
      }

      spin_coupling_square_operator_pairs[k] = std::make_pair(product_op_a, product_op_b);
      spin_coupling_square_coefficients[k] = spin_coupling_coefficients[i] * spin_coupling_coefficients[j];
    }
  }

  for (size_t i = 0; i < L - 1; ++i) {
    for (size_t term = 0; term < 144; term++) {
      if (spin_coupling_square_operator_pairs[term].first == Tensor({pb_in, pb_out})) {
        continue;
      }
      mpo_gen.AddTerm(0.5 * spin_coupling_square_coefficients[term],
                      spin_coupling_square_operator_pairs[term].first,
                      i,
                      spin_coupling_square_operator_pairs[term].second,
                      i + 1);//J = 1.0
    }
  }
  auto mpo = mpo_gen.Gen();

  gqmps2::TwoSiteVMPSSweepParams sweep_params(
      params.Sweeps,
      params.Dmin, params.Dmax, params.CutOff,
      gqmps2::LanczosParams(params.LanczErr, params.MaxLanczIter),
      params.noise
  );

  FiniteMPST mps(sites);
  std::vector<size_t> stat_labs;
  // suppose L % 16 == 0
  if (L % 16 == 0) {
    for (size_t i = 0; i < L; ++i) { stat_labs.push_back(i % 16); }
  } else if (L % 2 == 0) {
    for (size_t i = 0; i < L; ++i) {
      stat_labs.push_back((i % 2) * 15);
    }
  }

  if (world.rank() == 0) {
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
  }

  gqten::hp_numeric::SetTensorTransposeNumThreads(params.Threads);
  gqten::hp_numeric::SetTensorManipulationThreads(params.Threads);

  double e0;
  if (world.size() == 1) {
    e0 = gqmps2::TwoSiteFiniteVMPS(mps, mpo, sweep_params);
  } else {
    e0 = gqmps2::TwoSiteFiniteVMPS(mps, mpo, sweep_params, world);
  }

  if (world.rank() == 0) {
    std::cout << "E0/site: " << e0 / L << std::endl;    // E0/site: -0.43412
  }
  if (world.rank() == 0) {
    mps.Load(sweep_params.mps_path);
    auto ee_list = mps.GetEntanglementEntropy(1);
    std::copy(ee_list.begin(), ee_list.end(), std::ostream_iterator<double>(std::cout, " "));

    std::string file_name = "ee";
    std::ofstream ofs(file_name, std::ofstream::binary);
    ofs.write((const char *) ee_list.data(), ee_list.size() * sizeof(double));
    ofs << std::endl;
    ofs.close();

    std::vector<std::vector<size_t>> sites_set;
    for (size_t i = L / 4 + 1; i < std::min(3 * L / 4 + 2, L); i++) {
      sites_set.push_back({L / 4, i});
    }

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
    MeasureTwoSiteOp(mps, {S2op, S2op}, id, sites_set, "S2sym");

    Tensor Sp_prod = Tensor({pb_in, pb_out});
    Sp_prod({0b1111, 0b0000}) = 16.0;
    Tensor Sm_prod = Tensor({pb_in, pb_out});
    Sm_prod({0b0000, 0b1111}) = 16.0;

    MeasureTwoSiteOp(mps, {Sp_prod, Sm_prod}, id, sites_set, "PM");
    MeasureTwoSiteOp(mps, {Sm_prod, Sp_prod}, id, sites_set, "MP");
  }

  return 0;
}
