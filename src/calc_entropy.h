/*
 * File Name: calc_entropy.h
 * Description: some functionality entanglement entropy, mutual information and so on for ground state of XXZ model with quantum channel.
 * Created by Hao-Xin on 2023/05/08.
 */
#ifndef CHANNELING_QUANTUM_CRITICALITY_SRC_CALC_ENTROPY_H
#define CHANNELING_QUANTUM_CRITICALITY_SRC_CALC_ENTROPY_H
#define GQTEN_COUNT_FLOPS 1
#include "hei_type.h"

/**
 *
 * @param lambda
 * @return 1, 3 indexes is the physical bond, 2 index is the virtual bond
 */
Tensor GenerateChannelZMPOTensor(const double lambda) {
  IndexT virtual_bond = IndexT({QNSctT(qn0, 2)}, gqten::GQTenIndexDirType::OUT);
  Tensor mpo_ten = Tensor({pb_in, virtual_bond, pb_out});
  //Id \times Id
  mpo_ten({0, 0, 0}) = sqrt(1 - lambda);
  mpo_ten({1, 0, 1}) = sqrt(1 - lambda);

  //sigma_z \times sigma_z
  mpo_ten({0, 1, 0}) = sqrt(lambda);
  mpo_ten({1, 1, 1}) = -sqrt(lambda);

  return mpo_ten;
}

Tensor GenerateChannelXYZMPOTensor(const double lambda) {
  auto qn2up = QNT({gqten::QNCard("Sz", gqten::U1QNVal(2))});
  auto qn2dn = QNT({gqten::QNCard("Sz", gqten::U1QNVal(-2))});
  IndexT virtual_bond = IndexT({QNSctT(qn0, 2), QNSctT(qn2up, 1), QNSctT(qn2dn, 1)}, gqten::GQTenIndexDirType::OUT);
  Tensor mpo_ten = Tensor({pb_in, virtual_bond, pb_out});
  //Id
  mpo_ten({0, 0, 0}) = sqrt(1 - lambda);
  mpo_ten({1, 0, 1}) = sqrt(1 - lambda);

  //sigma_z
  mpo_ten({0, 1, 0}) = sqrt(lambda / 3);
  mpo_ten({1, 1, 1}) = -sqrt(lambda / 3);

  // sigma_plus
  mpo_ten({0, 2, 1}) = sqrt(2 * lambda / 3);

  // sigma_minus
  mpo_ten({1, 3, 0}) = sqrt(2 * lambda / 3);

  assert(mpo_ten.Div() == qn0);
  return mpo_ten;
}

Tensor GenerateChannelMPOTensor(const double lambda, const std::string type) {
  if (type == "z") {
    return GenerateChannelZMPOTensor(lambda);
  } else if (type == "xyz") {
    return GenerateChannelXYZMPOTensor(lambda);
  } else {
    std::cout << "do not support this channel" << std::endl;
    exit(1);
  }

}

/** for zz channel
 *
 * @param mps_ten
 * @param channel_mpo_tensor
 * @return
 *  two tensors, each has order like
 *
 *     2
 *     |
 * 0--MPO---3
 *     |
 *     1
 * 1, 2 is the virtual bond of MPS. 0, 3 will be contracted if one want the explict form of transfer matrix.
 */
std::vector<Tensor> GenerateDensityMatrixTransferMatrix(
    const Tensor &mps_ten,
    const Tensor &channel_mpo_tensor //get by GenerateChannelZMPOTensor(v), GenerateChannelXYZMPOTensor(v)
) {
  std::vector<Tensor> transfer_matrix(2);
  Contract(&mps_ten, &channel_mpo_tensor, {{1}, {0}}, &transfer_matrix[1]);
  transfer_matrix[1].Transpose({2, 0, 1,
                                3});// then 0 is the virtual bond of channel MPO, 1, 2 is the virtual bond of mps
  // 3 is the physics bond at start.
  transfer_matrix[0] = Dag(transfer_matrix[1]);
  transfer_matrix[0].Transpose({3, 1, 2, 0});
  return transfer_matrix;
}

///< Rho measure density matrix after Depolarizing channel
Tensor UpdateBoundaryTensorInRhoSquare(
    const Tensor &mps_ten,
    const Tensor &channel_mpo_tensor,
    const Tensor &boundary_tensor   //rank = 4, each dimension chi is the bond dimension of MPS
) {
  Tensor boundary_tensor_next;
  std::vector<Tensor> transfer_matrix = GenerateDensityMatrixTransferMatrix(mps_ten, channel_mpo_tensor);

  Tensor temp1, temp2, temp3;
  Contract(&boundary_tensor, &transfer_matrix[0], {{0}, {1}}, &temp1);
  Contract(&temp1, &transfer_matrix[1], {{5, 0}, {0, 1}}, &temp2);
  Contract(&temp2, &transfer_matrix[0], {{5, 0}, {0, 1}}, &temp3);
  Contract(&temp3, &transfer_matrix[1], {{5, 0, 1}, {0, 1, 3}}, &boundary_tensor_next);
  return boundary_tensor_next;
}

///< Rho measure density matrix after Depolarizing channel
Tensor UpdateBoundaryTensorInRhoCube(
    const Tensor &mps_ten,
    const Tensor &channel_mpo_tensor,
    const Tensor &boundary_tensor   //rank = 6, each dimension chi is the bond dimension of MPS
) {
  Tensor boundary_tensor_next;
  std::vector<Tensor> transfer_matrix = GenerateDensityMatrixTransferMatrix(mps_ten, channel_mpo_tensor);

  Tensor temp1, temp2, temp3, temp4, temp5;
  Contract(&boundary_tensor, &transfer_matrix[0], {{0}, {1}}, &temp1);
  Contract(&temp1, &transfer_matrix[1], {{7, 0}, {0, 1}}, &temp2);
  Contract(&temp2, &transfer_matrix[0], {{7, 0}, {0, 1}}, &temp3);
  Contract(&temp3, &transfer_matrix[1], {{7, 0}, {0, 1}}, &temp4);
  Contract(&temp4, &transfer_matrix[0], {{7, 0}, {0, 1}}, &temp5);
  Contract(&temp5, &transfer_matrix[1], {{7, 0, 1}, {0, 1, 3}}, &boundary_tensor_next);
  return boundary_tensor_next;
}

void CalRenyi3(
    const FiniteMPST &mps,
    const double v,
    double &thermal_renyi3, //output
    double &negativity_3, //output, half of the system
    std::string type = "z"
) {
  using namespace gqmps2;
  size_t L = mps.size();
  assert(L % 2 == 0);
  auto left_virtual_index = mps[0].GetIndex(0);
  auto left_virtual_index_inv = InverseIndex(left_virtual_index);
  Tensor left_boundary_ten =
      Tensor({left_virtual_index, left_virtual_index_inv, left_virtual_index, left_virtual_index_inv,
              left_virtual_index, left_virtual_index_inv});
  left_boundary_ten({0, 0, 0, 0, 0, 0}) = 1.0;
  Tensor channel_mpo_tensor = GenerateChannelMPOTensor(v, type);

  for (size_t i = 0; i < L / 2; i++) {
    Timer update_boundary_tensor("update_boundary_tensor");
    size_t flop_start = flop;
    left_boundary_ten = (UpdateBoundaryTensorInRhoCube(mps[i], channel_mpo_tensor, left_boundary_ten));
    double update_time = update_boundary_tensor.Elapsed();
    double Gflop = double(flop - flop_start) * 1e-9 / update_time;
    std::cout << "i = " << i << ", " << "time = " << update_time << ", Gflops = " << Gflop << std::endl;
  }

  FiniteMPST mps_transpose = mps;
  mps_transpose.Reverse();
  left_virtual_index = mps_transpose[0].GetIndex(0);
  left_virtual_index_inv = InverseIndex(left_virtual_index);
  Tensor right_boundary_ten =
      Tensor({left_virtual_index, left_virtual_index_inv, left_virtual_index, left_virtual_index_inv,
              left_virtual_index, left_virtual_index_inv});
  right_boundary_ten({0, 0, 0, 0, 0, 0}) = 1.0;
  for (size_t i = 0; i < L / 2; i++) {
    Timer update_boundary_tensor("update_boundary_tensor");
    right_boundary_ten = (UpdateBoundaryTensorInRhoCube(mps_transpose[i], channel_mpo_tensor, right_boundary_ten));
    double update_time = update_boundary_tensor.Elapsed();
    std::cout << "i = " << i << ", " << "time = " << update_time << std::endl;
  }
  Tensor overlap_ten, overlap_ten2;
  Contract(&left_boundary_ten, &right_boundary_ten, {{0, 1, 2, 3, 4, 5}, {0, 1, 2, 3, 4, 5}}, &overlap_ten);
  thermal_renyi3 = -std::log(overlap_ten()) / 2;
  right_boundary_ten.Dag();
  Contract(&left_boundary_ten, &right_boundary_ten, {{0, 1, 2, 3, 4, 5}, {1, 0, 3, 2, 5, 4}}, &overlap_ten2);
  negativity_3 = -std::log(overlap_ten2() / overlap_ten()) / 2;
}

void CalRenyi2(
    FiniteMPST &mps,
    const double v,
    double &thermal_renyi2, //output
    std::vector<double> &renyi2_a, //output, renyi entropy in A subsystem
    std::vector<double> &renyi2_b, //output
    std::string type = "z"
) {
  using namespace gqmps2;
  size_t L = mps.size();
  mps.Centralize(0);
  renyi2_a = std::vector<double>(L - 1);
  renyi2_b = std::vector<double>(L - 1);
  assert(L % 2 == 0);
  {
    auto left_virtual_index = mps[0].GetIndex(0);
    auto left_virtual_index_inv = InverseIndex(left_virtual_index);
    Tensor left_boundary_ten =
        Tensor({left_virtual_index, left_virtual_index_inv, left_virtual_index, left_virtual_index_inv});
    left_boundary_ten({0, 0, 0, 0}) = 1.0;
    Tensor channel_mpo_tensor = GenerateChannelMPOTensor(v, type);

    for (size_t i = 0; i < L; i++) {
      Timer update_boundary_tensor("update_boundary_tensor");
      size_t flop_start = flop;
      left_boundary_ten = (UpdateBoundaryTensorInRhoSquare(mps[i], channel_mpo_tensor, left_boundary_ten));
      double update_time = update_boundary_tensor.Elapsed();
      double Gflop = double(flop - flop_start) * 1.e-9 / update_time;
      std::cout << "i = " << i << ", " << "time = " << update_time << ", Gflops = " << Gflop << std::endl;
      if (i < L - 1) {
        double sum = 0.0;
        for (size_t k = 0; k < left_boundary_ten.GetShape()[0]; k++) {
          for (size_t j = 0; j < left_boundary_ten.GetShape()[2]; j++) {
            sum += left_boundary_ten({k, k, j, j});
          }
        }
        renyi2_a[i] = -std::log(sum);
      }
    }
    thermal_renyi2 = -std::log(left_boundary_ten({0, 0, 0, 0}));
  }

  mps.Reverse();
  {
    mps.Centralize(0);
    auto left_virtual_index = mps[0].GetIndex(0);
    auto left_virtual_index_inv = InverseIndex(left_virtual_index);
    Tensor left_boundary_ten =
        Tensor({left_virtual_index, left_virtual_index_inv, left_virtual_index, left_virtual_index_inv});
    left_boundary_ten({0, 0, 0, 0}) = 1.0;
    Tensor channel_mpo_tensor = GenerateChannelMPOTensor(v, type);

    for (size_t i = 0; i < L; i++) {
      Timer update_boundary_tensor("update_boundary_tensor");
      size_t flop_start = flop;
      left_boundary_ten = (UpdateBoundaryTensorInRhoSquare(mps[i], channel_mpo_tensor, left_boundary_ten));
      double update_time = update_boundary_tensor.Elapsed();
      double Gflop = double(flop - flop_start) * 1.e-9 / update_time;
      std::cout << "i = " << i << ", " << "time = " << update_time << ", Gflops = " << Gflop << std::endl;
      if (i < L - 1) {
        double sum = 0.0;
        for (size_t k = 0; k < left_boundary_ten.GetShape()[0]; k++) {
          for (size_t j = 0; j < left_boundary_ten.GetShape()[2]; j++) {
            sum += left_boundary_ten({k, k, j, j});
          }
        }
        renyi2_b[i] = -std::log(sum);
      }
    }
  }
  mps.Reverse();
}

#endif //CHANNELING_QUANTUM_CRITICALITY_SRC_CALC_ENTROPY_H