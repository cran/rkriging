#include <Rcpp.h>
#include <nloptrAPI.h>
// [[Rcpp::depends(nloptr)]]
#include <unordered_map>
#include "utils.h"

// Mapping from algorithm name strings to nlopt_algorithm enum values
static const std::unordered_map<std::string, nlopt_algorithm> ALGORITHM_MAP = {
    // Global derivative-free
    {"NLOPT_GN_DIRECT", NLOPT_GN_DIRECT},
    {"NLOPT_GN_DIRECT_L", NLOPT_GN_DIRECT_L},
    {"NLOPT_GN_DIRECT_L_RAND", NLOPT_GN_DIRECT_L_RAND},
    {"NLOPT_GN_DIRECT_NOSCAL", NLOPT_GN_DIRECT_NOSCAL},
    {"NLOPT_GN_DIRECT_L_NOSCAL", NLOPT_GN_DIRECT_L_NOSCAL},
    {"NLOPT_GN_DIRECT_L_RAND_NOSCAL", NLOPT_GN_DIRECT_L_RAND_NOSCAL},
    {"NLOPT_GN_ORIG_DIRECT", NLOPT_GN_ORIG_DIRECT},
    {"NLOPT_GN_ORIG_DIRECT_L", NLOPT_GN_ORIG_DIRECT_L},
    
    // Global derivative-based
    {"NLOPT_GD_STOGO", NLOPT_GD_STOGO},
    {"NLOPT_GD_STOGO_RAND", NLOPT_GD_STOGO_RAND},
    
    // Local derivative-based
    {"NLOPT_LD_SLSQP", NLOPT_LD_SLSQP},
    {"NLOPT_LD_LBFGS", NLOPT_LD_LBFGS},
    {"NLOPT_LD_VAR1", NLOPT_LD_VAR1},
    {"NLOPT_LD_VAR2", NLOPT_LD_VAR2},
    {"NLOPT_LD_TNEWTON", NLOPT_LD_TNEWTON},
    {"NLOPT_LD_TNEWTON_RESTART", NLOPT_LD_TNEWTON_RESTART},
    {"NLOPT_LD_TNEWTON_PRECOND", NLOPT_LD_TNEWTON_PRECOND},
    {"NLOPT_LD_TNEWTON_PRECOND_RESTART", NLOPT_LD_TNEWTON_PRECOND_RESTART},
    {"NLOPT_LD_MMA", NLOPT_LD_MMA},
    {"NLOPT_LD_CCSAQ", NLOPT_LD_CCSAQ},
    {"NLOPT_LD_AUGLAG", NLOPT_LD_AUGLAG},
    {"NLOPT_LD_AUGLAG_EQ", NLOPT_LD_AUGLAG_EQ},
    
    // Local derivative-free
    {"NLOPT_LN_PRAXIS", NLOPT_LN_PRAXIS},
    {"NLOPT_LN_COBYLA", NLOPT_LN_COBYLA},
    {"NLOPT_LN_NEWUOA", NLOPT_LN_NEWUOA},
    {"NLOPT_LN_NEWUOA_BOUND", NLOPT_LN_NEWUOA_BOUND},
    {"NLOPT_LN_NELDERMEAD", NLOPT_LN_NELDERMEAD},
    {"NLOPT_LN_SBPLX", NLOPT_LN_SBPLX},
    {"NLOPT_LN_BOBYQA", NLOPT_LN_BOBYQA},
    
    // Global derivative-free with local
    {"NLOPT_GN_CRS2_LM", NLOPT_GN_CRS2_LM},
    {"NLOPT_GN_MLSL", NLOPT_GN_MLSL},
    {"NLOPT_GN_MLSL_LDS", NLOPT_GN_MLSL_LDS},
    {"NLOPT_GN_ISRES", NLOPT_GN_ISRES},
    
    // Global derivative-based with local
    {"NLOPT_GD_MLSL", NLOPT_GD_MLSL},
    {"NLOPT_GD_MLSL_LDS", NLOPT_GD_MLSL_LDS},
    
    // Augmented Lagrangian
    {"NLOPT_LN_AUGLAG", NLOPT_LN_AUGLAG},
    {"NLOPT_LN_AUGLAG_EQ", NLOPT_LN_AUGLAG_EQ},
};

nlopt_opt _nlopt_init(const std::string& algorithm, const std::size_t& dim) {
    // Look up algorithm in map; use NLOPT_LN_SBPLX as default fallback
    auto it = ALGORITHM_MAP.find(algorithm);
    nlopt_algorithm algo = (it != ALGORITHM_MAP.end()) ? it->second : NLOPT_LN_SBPLX;
    return nlopt_create(algo, dim);
}
