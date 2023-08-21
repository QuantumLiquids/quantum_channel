/*
 * File Name: gqdouble.h
 * Description: Define the types and constant in Hubbard DMRG project.
 * Created by Hao-Xin on 2023/04/04.
 *
 *
 */

#ifndef HUBBARD_SRC_GQDOUBLE_H_
#define HUBBARD_SRC_GQDOUBLE_H_

#include "gqten/gqten.h"
#include "gqmps2/gqmps2.h"


using TenElemT = gqten::GQTEN_Double;



using gqten::GQTensor;

const std::string kMpoPath = "mpo";
const std::string kMpoTenBaseName = "mpo_ten";



using TenElemT = gqten::GQTEN_Double;
using QNT = gqten::special_qn::U1QN;
using QNSctT = gqten::QNSector<QNT>;
using IndexT = gqten::Index<QNT>;
using Tensor = gqten::GQTensor<TenElemT, QNT>;
using SiteVecT = gqmps2::SiteVec<TenElemT, QNT>;
using FiniteMPST = gqmps2::FiniteMPS<TenElemT, QNT>;


auto qn0 = QNT({gqten::QNCard("Sz", gqten::U1QNVal(0))});
auto qnup = QNT({gqten::QNCard("Sz", gqten::U1QNVal(1))});
auto qndn = QNT({gqten::QNCard("Sz", gqten::U1QNVal(-1))});
auto pb_out = IndexT(
    {QNSctT(qnup, 1), QNSctT(qndn, 1)},
    gqten::GQTenIndexDirType::OUT
);
auto pb_in = gqten::InverseIndex(pb_out);

#endif //HUBBARD_SRC_GQDOUBLE_H_
