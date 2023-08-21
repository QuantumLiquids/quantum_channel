/*
 * File Name: ising_type.h
 * Description: Define the types and constant in Ising model.
 * Created by Hao-Xin on 2023/04/04.
 *
 *
 */

#ifndef ISING_TYPE_H
#define ISING_TYPE_H

#include "gqten/gqten.h"
#include "gqmps2/gqmps2.h"

using TenElemT = gqten::GQTEN_Double;

using gqten::GQTensor;

const std::string kMpoPath = "mpo";
const std::string kMpoTenBaseName = "mpo_ten";

using TenElemT = gqten::GQTEN_Double;
using QNT = gqten::special_qn::U1U1ZnQN<2>;// Z2
using QNSctT = gqten::QNSector<QNT>;
using IndexT = gqten::Index<QNT>;
using Tensor = gqten::GQTensor<TenElemT, QNT>;
using SiteVecT = gqmps2::SiteVec<TenElemT, QNT>;
using FiniteMPST = gqmps2::FiniteMPS<TenElemT, QNT>;

auto qn0 = QNT(0, 0, 0);//only last number is valid;
auto qn1 = QNT(0, 0, 1);
auto pb_out = IndexT(
    {QNSctT(qn0, 1), QNSctT(qn1, 1)},
    gqten::GQTenIndexDirType::OUT
);
auto pb_in = gqten::InverseIndex(pb_out);

#endif //ISING_TYPE_H
