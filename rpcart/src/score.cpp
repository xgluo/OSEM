#include "score.h"

#include "bits.h"

namespace pcart {

StructureScoreTerms computeStructureScoreTerms(const vector<VarPtr>& predictors) {
	double c = 0.0;
	vector<size_t> realStateCounts;
	vector<size_t> catStateCounts;
	vector<size_t> ordStateCounts;
	for(const VarPtr& var : predictors) {
		lambdaVisit(var,
			[&](const RealVarPtr& var) {
				c += 1.0;
				realStateCounts.push_back(var->maxSubdiv + 1);
			},
			[&](const CatVarPtr& var) {
				c += (double)bit64(var->cats.size() - 1) - 1.0;
				catStateCounts.push_back(var->cats.size());
			},
			[&](const OrdVarPtr& var) {
				c += (double)var->cats.size() - 1.0;
				ordStateCounts.push_back(var->cats.size());
			}
		);
	}
	double term = 0.0;
	if(c != 0.0) {
		term = -log(4.0 * c);
	}

	// For every cell type (parameterized by the number of splits done in each
	// real variable and number of remaining categories in each categorical
	// variable), compute the logarithm sum of penalties over all subtrees
	// starting from a cell of such type.
	size_t stateCount = 1;
	for(size_t x : realStateCounts) {
		stateCount *= x;
	}
	for(size_t x : catStateCounts) {
		stateCount *= x;
	}
	for(size_t x : ordStateCounts) {
		stateCount *= x;
	}
	vector<double> sums(stateCount);
	for(size_t conf = 0; conf < stateCount; ++conf) {
		double val = term;
		size_t mul = 1;
		for(size_t x : realStateCounts) {
			size_t s = (conf / mul) % x;
			if(s > 0) {
				double sub = sums[conf - mul];
				val = addLog(val, sub + sub);
			}
			mul *= x;
		}
		for(size_t x : catStateCounts) {
			size_t s = (conf / mul) % x;
			for(size_t a = 1; a <= s; ++a) {
				double sub1 = sums[conf - a * mul];
				double sub2 = sums[conf - (s - a + 1) * mul];
				val = addLog(val, sub1 + sub2 + lognCr((double)s, (double)a));
			}
			mul *= x;
		}
		for(size_t x : ordStateCounts) {
			size_t s = (conf / mul) % x;
			for(size_t a = 1; a <= s; ++a) {
				double sub1 = sums[conf - a * mul];
				double sub2 = sums[conf - (s - a + 1) * mul];
				val = addLog(val, sub1 + sub2);
			}
			mul *= x;
		}
		sums[conf] = val;
	}

	return StructureScoreTerms{ term, -sums.back() };
}

}
