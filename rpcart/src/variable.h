#pragma once
#pragma once

#include "common.h"

namespace pcart {

struct BaseVar {
	string name;
	size_t dataSrcIdx;
};

struct RealVar : public BaseVar {
	typedef double Val;

	double minVal;
	double maxVal;
	size_t maxSubdiv;

	double nu;
	double lambda;
	double barmu;
	double a;

	double parseDataSrcVal(const vector<double>& src) const {
		if(dataSrcIdx >= src.size()) {
			fail("Too short vector in data source");
		}
		double val = src[dataSrcIdx];
		if(val < minVal || val > maxVal) {
			fail("Value ", val, " in data source is outside the allowed range [", minVal, ", ", maxVal, "] for variable ", name);
		}
		return val;
	}
};

struct CatInfo {
	string name;
	double alpha;
};

struct CatVar : public BaseVar {
	typedef size_t Val;

	vector<CatInfo> cats;

	// If enabled, the hyperparameter alpha is scaled by the size of the leaf cell in the data space.
	// Otherwise, the same hyperparameter is used in each leaf cell.
	bool bdeu;

	size_t parseDataSrcVal(const vector<double>& src) const {
		if(dataSrcIdx >= src.size()) {
			fail("Too short vector in data source");
		}
		double srcVal = src[dataSrcIdx];
		size_t val = (size_t)srcVal;
		if(srcVal != (double)val || val >= cats.size()) {
			fail("Invalid value ", srcVal, " for variable ", name, " in data source");
		}
		return val;
	}
};

struct OrdVar : public BaseVar {
	typedef size_t Val;

	vector<CatInfo> cats;
	double kappa;

	size_t parseDataSrcVal(const vector<double>& src) const {
		if(dataSrcIdx >= src.size()) {
			fail("Too short vector in data source");
		}
		double srcVal = src[dataSrcIdx];
		size_t val = (size_t)srcVal;
		if(srcVal != (double)val || val >= cats.size()) {
			fail("Invalid value ", srcVal, " for variable ", name, " in data source");
		}
		return val;
	}
};

typedef shared_ptr<const RealVar> RealVarPtr;
typedef shared_ptr<const CatVar> CatVarPtr;
typedef shared_ptr<const OrdVar> OrdVarPtr;
typedef variant<RealVarPtr, CatVarPtr, OrdVarPtr> VarPtr;

RealVarPtr createRealVar(
	string name,
	size_t dataSrcIdx, // Index of the value of the variable in data points vectors
	double minVal, double maxVal, // The range of allowed values for the variable
	size_t maxSubdiv, // Maximum number of times the variable can be subdivided into half
	// Hyperparameters (nu, lambda, barmu, a), the defaults are possibly ok for N(0, 1) normalized data
	double nu = 1.0,
	double lambda = 1.0,
	double barmu = 0.0,
	double a = 1.0
);

CatVarPtr createCatVar(
	string name,
	size_t dataSrcIdx, // Index of the value of the variable in data points vectors
	const vector<string>& catNames, // catNames[i] = the name of the category i
	double alpha = 0.5 // Hyperparameter alpha for every category, default is Jeffrey's prior
);
CatVarPtr createBDeuCatVar(
	string name,
	size_t dataSrcIdx, // Index of the value of the variable in data points vectors
	const vector<string>& catNames, // catNames[i] = the name of the category i
	double ess = 1.0 // Equivalent sample size
);

OrdVarPtr createOrdVar(
	string name,
	size_t dataSrcIdx, // Index of the value of the variable in data points vectors
	const vector<string>& catNames, // catNames[i] = the name of the category i
	double alpha = 0.5, // Hyperparameter alpha for every category, default is Jeffrey's prior
	double kappa = 0.5
);

}
