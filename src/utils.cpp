#include <Rcpp.h>
#include <nloptrAPI.h>
// [[Rcpp::depends(nloptr)]]
#include "utils.h"

nlopt_opt nlopt_init(const std::string& algorithm, const::std::size_t& dim) {
    if (algorithm == "NLOPT_GN_DIRECT") return nlopt_create(NLOPT_GN_DIRECT, dim);
	if (algorithm == "NLOPT_GN_DIRECT_L") return nlopt_create(NLOPT_GN_DIRECT_L, dim);
	if (algorithm == "NLOPT_GN_DIRECT_L_RAND") return nlopt_create(NLOPT_GN_DIRECT_L_RAND, dim);
	if (algorithm == "NLOPT_GN_DIRECT_NOSCAL") return nlopt_create(NLOPT_GN_DIRECT_NOSCAL, dim);
	if (algorithm == "NLOPT_GN_DIRECT_L_NOSCAL") return nlopt_create(NLOPT_GN_DIRECT_L_NOSCAL, dim);
	if (algorithm == "NLOPT_GN_DIRECT_L_RAND_NOSCAL") return nlopt_create(NLOPT_GN_DIRECT_L_RAND_NOSCAL, dim);
	if (algorithm == "NLOPT_GN_ORIG_DIRECT") return nlopt_create(NLOPT_GN_ORIG_DIRECT, dim);
	if (algorithm == "NLOPT_GN_ORIG_DIRECT_L") return nlopt_create(NLOPT_GN_ORIG_DIRECT_L, dim);
	if (algorithm == "NLOPT_GD_STOGO") return nlopt_create(NLOPT_GD_STOGO, dim);
	if (algorithm == "NLOPT_GD_STOGO_RAND") return nlopt_create(NLOPT_GD_STOGO_RAND, dim);
	if (algorithm == "NLOPT_LD_SLSQP") return nlopt_create(NLOPT_LD_SLSQP, dim);
	if (algorithm == "NLOPT_LD_LBFGS_NOCEDAL") return nlopt_create(NLOPT_LD_LBFGS_NOCEDAL, dim);
	if (algorithm == "NLOPT_LD_LBFGS") return nlopt_create(NLOPT_LD_LBFGS, dim);
	if (algorithm == "NLOPT_LN_PRAXIS") return nlopt_create(NLOPT_LN_PRAXIS, dim);
	if (algorithm == "NLOPT_LD_VAR1") return nlopt_create(NLOPT_LD_VAR1, dim);
	if (algorithm == "NLOPT_LD_VAR2") return nlopt_create(NLOPT_LD_VAR2, dim);
	if (algorithm == "NLOPT_LD_TNEWTON") return nlopt_create(NLOPT_LD_TNEWTON, dim);
	if (algorithm == "NLOPT_LD_TNEWTON_RESTART") return nlopt_create(NLOPT_LD_TNEWTON_RESTART, dim);
	if (algorithm == "NLOPT_LD_TNEWTON_PRECOND") return nlopt_create(NLOPT_LD_TNEWTON_PRECOND, dim);
	if (algorithm == "NLOPT_LD_TNEWTON_PRECOND_RESTART") return nlopt_create(NLOPT_LD_TNEWTON_PRECOND_RESTART, dim);
	if (algorithm == "NLOPT_GN_CRS2_LM") return nlopt_create(NLOPT_GN_CRS2_LM, dim);
	if (algorithm == "NLOPT_GN_MLSL") return nlopt_create(NLOPT_GN_MLSL, dim);
	if (algorithm == "NLOPT_GD_MLSL") return nlopt_create(NLOPT_GD_MLSL, dim);
	if (algorithm == "NLOPT_GN_MLSL_LDS") return nlopt_create(NLOPT_GN_MLSL_LDS, dim);
	if (algorithm == "NLOPT_GD_MLSL_LDS") return nlopt_create(NLOPT_GD_MLSL_LDS, dim);
	if (algorithm == "NLOPT_LD_MMA") return nlopt_create(NLOPT_LD_MMA, dim);
	if (algorithm == "NLOPT_LD_CCSAQ") return nlopt_create(NLOPT_LD_CCSAQ, dim);
	if (algorithm == "NLOPT_LN_COBYLA") return nlopt_create(NLOPT_LN_COBYLA, dim);
	if (algorithm == "NLOPT_LN_NEWUOA") return nlopt_create(NLOPT_LN_NEWUOA, dim);
	if (algorithm == "NLOPT_LN_NEWUOA_BOUND") return nlopt_create(NLOPT_LN_NEWUOA_BOUND, dim);
	if (algorithm == "NLOPT_LN_NELDERMEAD") return nlopt_create(NLOPT_LN_NELDERMEAD, dim);
	if (algorithm == "NLOPT_LN_SBPLX") return nlopt_create(NLOPT_LN_SBPLX, dim);
	if (algorithm == "NLOPT_LN_AUGLAG") return nlopt_create(NLOPT_LN_AUGLAG, dim);
	if (algorithm == "NLOPT_LD_AUGLAG") return nlopt_create(NLOPT_LD_AUGLAG, dim);
	if (algorithm == "NLOPT_LN_AUGLAG_EQ") return nlopt_create(NLOPT_LN_AUGLAG_EQ, dim);
	if (algorithm == "NLOPT_LD_AUGLAG_EQ") return nlopt_create(NLOPT_LD_AUGLAG_EQ, dim);
	if (algorithm == "NLOPT_LN_BOBYQA") return nlopt_create(NLOPT_LN_BOBYQA, dim);
	if (algorithm == "NLOPT_GN_ISRES") return nlopt_create(NLOPT_GN_ISRES, dim);
    // default if algorithm does not match any
    return nlopt_create(NLOPT_LN_SBPLX, dim); 
}
