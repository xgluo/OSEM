#pragma once

#include "score.h"
#include "variable.h"

namespace pcart {

template <typename T>
struct Leaf;
typedef Leaf<RealVar> RealLeaf;
typedef Leaf<CatVar> CatLeaf;
typedef Leaf<OrdVar> OrdLeaf;

template <typename T>
struct Split;
typedef Split<RealVar> RealSplit;
typedef Split<CatVar> CatSplit;
typedef Split<OrdVar> OrdSplit;

typedef std::variant<RealLeaf, CatLeaf, OrdLeaf, RealSplit, CatSplit, OrdSplit> Tree;
typedef std::unique_ptr<Tree> TreePtr;

struct BaseLeaf {};

template <typename T>
struct Leaf : public BaseLeaf {
	shared_ptr<const T> var;
	LeafStats<T> stats;

	template <typename F>
	Leaf(const shared_ptr<const T>& var, size_t dataCount, F getVal)
		: var(var),
		  stats(*var, dataCount, getVal)
	{}
};

struct BaseSplit {
	unique_ptr<Tree> leftChild;
	unique_ptr<Tree> rightChild;
};

template <>
struct Split<RealVar> : public BaseSplit {
	RealVarPtr var;
	double splitVal;
};

template <>
struct Split<CatVar> : public BaseSplit {
	CatVarPtr var;
	uint64_t leftCatMask;
	uint64_t rightCatMask;
};

template <>
struct Split<OrdVar> : public BaseSplit {
	OrdVarPtr var;
	uint64_t leftCatMask;
	uint64_t rightCatMask;
};

struct TreeResult {
	double dataScore;
	double structureScore;
	TreePtr tree;

	double totalScore() const {
		return dataScore + structureScore;
	}
};

// Returns the optimum tree for given predictor variables, response variable and data.
// In each point vector data[i], the value for each variable is given by the element
// in the index defined by the dataSrcIdx field of the variable. The values for a
// categorical variable with n categories should be represented by integer values
// 0, 1, ..., n-1. Set useStructureScore to false to set all structure scores to zero.
TreeResult optimizeTree(
	const vector<VarPtr>& predictors,
	const VarPtr& response,
	const vector<vector<double>>& data,
	bool useStructureScore = true
);

// Like optimizeTree, but returns the tree corresponding to the full table with
// given predictors (the tree that makes all the possible splits). The only part
// that optimized is the discretization: each real predictor is discretized into
// 2..maxDiscretizationBins quantiles. The structure score of the result is
// always 0.
TreeResult optimizeFullTable(
	const vector<VarPtr>& predictors,
	const VarPtr& response,
	const vector<vector<double>>& data,
	size_t maxDiscretizationBins
);

void iterateTrees(
	const vector<VarPtr>& predictors,
	const VarPtr& response,
	const vector<vector<double>>& data,
	function<void(const TreeResult&)> f
);

void printTree(const TreePtr& tree, ostream& out = cout);

}
