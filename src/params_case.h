/*
 * File Name: params_case.h
 * Description: Declare CaseParams class used set parameters by users
 * Created by Hao-Xin on 2021/10/15.
 *
 */


#ifndef HEISENBERG_SRC_PARAMS_CASE_H
#define HEISENBERG_SRC_PARAMS_CASE_H

#include "gqmps2/case_params_parser.h"
using gqmps2::CaseParamsParserBasic;

struct CaseParams : public CaseParamsParserBasic {
  CaseParams(const char *pf) : CaseParamsParserBasic(pf) {
    Model = ParseStr("Model");
    L= ParseInt("L");
    Delta = ParseDouble("Delta");
    h = ParseDouble("h");
    gamma = ParseDouble("gamma");
    Sweeps = ParseInt("Sweeps");
    Dmin = ParseInt("Dmin");
    Dmax = ParseInt("Dmax");
    CutOff = ParseDouble("CutOff");
    LanczErr = ParseDouble("LanczErr");
    MaxLanczIter = ParseInt("MaxLanczIter");
    Threads = ParseInt("Threads");
    noise = ParseDoubleVec("noise");
  }

  std::string Model; //Ising, and XXZ
  size_t L;
  double Delta;
  double h;
  double gamma;
  size_t Sweeps;
  size_t Dmin;
  size_t Dmax;
  double CutOff;
  double LanczErr;
  size_t MaxLanczIter;
  size_t Threads;
  std::vector<double> noise;
};

#endif //HEISENBERG_SRC_PARAMS_CASE_H
