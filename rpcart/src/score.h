#pragma once

#include "variable.h"

namespace pcart {

template <typename T>
struct LeafStats;

template <>
struct LeafStats<RealVar> {
	double avg;
	double stddev;
	size_t dataCount;

	template <typename F>
	LeafStats(const RealVar& var, size_t dataCount, F getVal)
		: avg(0.0),
		  stddev(0.0),
		  dataCount(dataCount)
	{
		for(size_t i = 0; i < dataCount; ++i) {
			avg += getVal(i);
		}
		avg /= (double)dataCount;
		for(size_t i = 0; i < dataCount; ++i) {
			double d = getVal(i) - avg;
			stddev += d * d;
		}
		stddev = sqrt(stddev / (double)dataCount);
	}

	double dataScore(const RealVar& var, double cellSize) const {
		if(dataCount == 0) {
			return 0.0;
		}

		double n = (double)dataCount;

		double mu = avg;
		double vari = stddev * stddev;

		double s = (n - 1.0) * vari;
		double mud = mu - var.barmu;
		double t = (n * var.a / (n + var.a)) * mud * mud;

		double score = 0.0;
		score -= 0.5 * n * log(pi);
		score += 0.5 * var.nu * log(var.lambda * var.nu);
		score += 0.5 * log(var.a);
		score -= 0.5 * log(n + var.a);
		score += lgamma(0.5 * (n + var.nu));
		score -= lgamma(0.5 * var.nu);
		score -= 0.5 * (n + var.nu) * log(s + t + var.nu * var.lambda);

		return score;
	}
};
template <>
struct LeafStats<CatVar> {
	std::vector<size_t> catCount;
	size_t dataCount;

	template <typename F>
	LeafStats(const CatVar& var, size_t dataCount, F getVal)
		: catCount(var.cats.size(), 0),
		  dataCount(dataCount)
	{
		for(size_t i = 0; i < dataCount; ++i) {
			++catCount[getVal(i)];
		}
	}

	double dataScore(const CatVar& var, double cellSize) const {
		if(dataCount == 0) {
			return 0.0;
		}

		double coef = var.bdeu ? cellSize : 1.0;

		double score = 0.0;
		
		double alphaSum = 0.0;
		for(const CatInfo& info : var.cats) {
			score -= lgamma(info.alpha * coef);
			alphaSum += info.alpha * coef;
		}
		score += lgamma(alphaSum);

		for(size_t i = 0; i < var.cats.size(); ++i) {
			score += lgamma((double)catCount[i] + var.cats[i].alpha * coef);
		}
		score -= lgamma((double)dataCount + alphaSum);

		return score;
	}
};
template <>
struct LeafStats<OrdVar> {
	std::vector<size_t> catCount;
	size_t dataCount;

	template <typename F>
	LeafStats(const OrdVar& var, size_t dataCount, F getVal)
		: catCount(var.cats.size(), 0),
		  dataCount(dataCount)
	{
		for(size_t i = 0; i < dataCount; ++i) {
			++catCount[getVal(i)];
		}
	}

	double dataScore(const OrdVar& var, double cellSize) const {
		double alphaSum = 0.0;
		for(const CatInfo& info : var.cats) {
			alphaSum += info.alpha;
		}

		double logKappa = log(var.kappa);

		if(var.cats.size() > 63) {
			fail("More than 63 categories in an ordinal variable not supported");
		}
		double suffixes[64];
		suffixes[var.cats.size()] = 0.0;

		size_t segStart = var.cats.size();
		while(segStart --> 0) {
			suffixes[segStart] = -numeric_limits<double>::infinity();

			double segAlpha = 0.0;
			size_t segCount = 0;
			double base = logKappa;

			for(size_t segEnd = segStart + 1; segEnd <= var.cats.size(); ++segEnd) {
				size_t segSize = segEnd - segStart;
				segAlpha += var.cats[segEnd - 1].alpha;
				segCount += catCount[segEnd - 1];
				base -= lgamma((double)catCount[segEnd - 1] + 1.0);

				double factor = base;
				factor += lgamma((double)segCount + segAlpha);
				factor -= lgamma(segAlpha);
				factor += lgamma((double)segCount + 1.0);
				factor -= (double)segCount * log((double)segSize);

				suffixes[segStart] = addLog(suffixes[segStart], suffixes[segEnd] + factor);
			}
		}

		double score = suffixes[0];
		score += lgamma(alphaSum);
		score -= lgamma((double)dataCount + alphaSum);
		return score;
	}
};

struct StructureScoreTerms {
	double leafPenaltyTerm;
	double normalizerTerm;
};

StructureScoreTerms computeStructureScoreTerms(const vector<VarPtr>& predictors);

}
