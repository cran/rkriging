#ifndef UTILS_H
#define UTILS_H

#include <Rcpp.h>
#include <nloptrAPI.h>
// [[Rcpp::depends(nloptr)]]

nlopt_opt nlopt_init(const std::string& algorithm, const::std::size_t& dim);

#endif

