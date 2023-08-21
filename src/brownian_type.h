//
// Created by Hao-Xin on 2023/6/25.
//

#ifndef CHANNELING_QUANTUM_CRITICALITY_SRC_BROWNIAN_TYPE_H
#define CHANNELING_QUANTUM_CRITICALITY_SRC_BROWNIAN_TYPE_H

#include "gqten/gqten.h"
#include "gqmps2/gqmps2.h"
#include "myutil.h"
#include "u1u1u1u1qn.h"
//#include <numeric>

using TenElemT = gqten::GQTEN_Double;

using gqten::GQTensor;

const std::string kMpoPath = "mpo";
const std::string kMpoTenBaseName = "mpo_ten";

using TenElemT = gqten::GQTEN_Double;
//using QNT = gqten::QN<gqten::U1QNVal, gqten::U1QNVal, gqten::U1QNVal, gqten::U1QNVal>;
using QNT = gqten::special_qn::U1U1U1U1QN;
using QNSctT = gqten::QNSector<QNT>;
using IndexT = gqten::Index<QNT>;
using Tensor = gqten::GQTensor<TenElemT, QNT>;
using SiteVecT = gqmps2::SiteVec<TenElemT, QNT>;
using FiniteMPST = gqmps2::FiniteMPS<TenElemT, QNT>;

const size_t num_U1 = 4;
const size_t num_qnscts = 16; //2^num_U1
std::vector<gqten::QNCard> qn_cards(num_U1, gqten::QNCard(" ", gqten::U1QNVal(0)));
auto qn0 = QNT(qn_cards);

const size_t local_hilbert_dim = 16;
QNT qn_physical[num_qnscts];
IndexT pb_out, pb_in;

inline void InitialLocalHilbertSpace() {
  for (size_t i = 0; i < num_qnscts; i++) {
    auto binary_form = BinaryTransform(i, num_U1);
    std::vector<gqten::QNCard> qn_cards(num_U1, gqten::QNCard(" ", gqten::U1QNVal(0)));
    for (size_t j = 0; j < num_U1; j++) {
      int sz_val = (2 * int(binary_form[j])) - 1;
      qn_cards[j] = gqten::QNCard(" ", gqten::U1QNVal(sz_val));
    }
    qn_physical[i] = QNT(qn_cards);
  }
  std::vector<QNSctT> qn_scts(num_qnscts);
  for (size_t i = 0; i < num_qnscts; i++) {
    qn_scts[i] = QNSctT(qn_physical[i], local_hilbert_dim/num_qnscts);
  }
  pb_out = IndexT(
      {qn_scts},
      gqten::GQTenIndexDirType::OUT
  );
  pb_in = gqten::InverseIndex(pb_out);
//  pb_out.Show();
}

#endif //CHANNELING_QUANTUM_CRITICALITY_SRC_BROWNIAN_TYPE_H
