#pragma once

#include "variable.h"
#include "bits.h"
#include "tree.h"

namespace pcart {

class Cell {
public:
	Cell() {}

	bool operator==(Cell other) const {
		return repr == other.repr;
	}

private:
	uint64_t repr;

	friend class CellCtx;
	friend class DataSplitter;
	friend struct std::hash<Cell>;
};

class DataSplitter {
public:
	DataSplitter() {}

	template <typename T>
	size_t split(pair<Cell, T>* data, size_t dataCount) const {
		size_t a = 0;
		size_t b = dataCount;
		while(a != b) {
			if(data[a].first.repr & mask) {
				--b;
				swap(data[a], data[b]);
			} else {
				++a;
			}
		}
		return a;
	}

private:
	uint64_t mask;

	friend class CellCtx;
};

class CellCtx {
public:
	CellCtx(const vector<VarPtr>& vars) {
		size_t bitPos = 0;
		for(const VarPtr& var : vars) {
			size_t bitCount = lambdaVisit(var,
				[&](const RealVarPtr& var) {
					return var->maxSubdiv + 1;
				},
				[&](const CatVarPtr& var) {
					return var->cats.size();
				},
				[&](const OrdVarPtr& var) {
					return var->cats.size();
				}
			);
			varInfo_.push_back({var, bitPos, bitCount, 0});
			bitPos += bitCount;
		}
		
		if(bitPos > 64) {
			fail("Storing a cell requires ", bitPos, " bits, but CellCtx supports at most 64");
		}
		
		realMask_ = 0;
		catOrdMask_ = 0;
		root_.repr = 0;
		for(VarInfo& info : varInfo_) {
			info.mask = ones64(info.bitCount) << info.startBit;
			lambdaVisit(info.var,
				[&](const RealVarPtr& var) {
					realMask_ |= info.mask;
					info.putRepr(root_, (uint64_t)-1 << 1);
				},
				[&](const CatVarPtr& var) {
					catOrdMask_ |= info.mask;
					info.putRepr(root_, (uint64_t)-1);
				},
				[&](const OrdVarPtr& var) {
					catOrdMask_ |= info.mask;
					info.putRepr(root_, (uint64_t)-1);
				}
			);
		}
	}

	Cell root() const {
		return root_;
	}

	Cell pointCell(const vector<double>& src) const {
		Cell cell;
		cell.repr = 0;
		for(const VarInfo& info : varInfo_) {
			lambdaVisit(info.var,
				[&](const RealVarPtr& var) {
					double val = var->parseDataSrcVal(src);
					double t = (val - var->minVal) / (var->maxVal - var->minVal);
					t *= (double)bit64(var->maxSubdiv);
					t = floor(t);
					t = max(t, 0.0);
					uint64_t repr = (uint64_t)t;
					repr = min(repr, ones64(var->maxSubdiv));
					info.putRepr(cell, repr);
				},
				[&](const CatVarPtr& var) {
					size_t val = var->parseDataSrcVal(src);
					info.putRepr(cell, bit64(val));
				},
				[&](const OrdVarPtr& var) {
					size_t val = var->parseDataSrcVal(src);
					info.putRepr(cell, bit64(val));
				}
			);
		}
		return cell;
	}

	template <typename F>
	void iterateSplits(Cell cell, F f) const {
		for(const VarInfo& info : varInfo_) {
			uint64_t repr = info.getRepr(cell);
			lambdaVisit(info.var,
				[&](const RealVarPtr& var) {
					uint64_t topBit = bit64(var->maxSubdiv);
					if(repr & topBit) {
						uint64_t leftRepr = (repr ^ topBit) << 1;
						uint64_t rightRepr = leftRepr | (uint64_t)1;

						Cell left = cell;
						info.putRepr(left, leftRepr);
						Cell right = cell;
						info.putRepr(right, rightRepr);

						size_t shift = clz64(~repr << (64 - var->maxSubdiv));
						DataSplitter dataSplitter;
						dataSplitter.mask = bit64(info.startBit + shift);

						auto lazySplit = [&](TreePtr leftChild, TreePtr rightChild) {
							size_t len = var->maxSubdiv - shift;
							double t = (double)(rightRepr & ones64(len)) / (double)bit64(len);
							t = max(min(t, 1.0), 0.0);
							double splitVal = (1.0 - t) * var->minVal + t * var->maxVal;
							
							return make_unique<Tree>(RealSplit{
								move(leftChild), move(rightChild),
								var, splitVal
							});
						};

						f(info.var, left, right, 0.5, 0.5, dataSplitter, lazySplit);
					}
				},
				[&](const CatVarPtr& var) {
					uint64_t bottom = bottomOne64(repr);
					uint64_t mask = repr ^ bottom;
					for(uint64_t sub = mask; sub; sub = (sub - 1) & mask) {
						uint64_t leftCatMask = repr ^ sub;
						uint64_t rightCatMask = sub;

						Cell left = cell;
						info.putRepr(left, leftCatMask);
						Cell right = cell;
						info.putRepr(right, rightCatMask);

						DataSplitter dataSplitter;
						dataSplitter.mask = right.repr & info.mask;

						auto lazySplit = [&](TreePtr leftChild, TreePtr rightChild) {
							return make_unique<Tree>(CatSplit{
								move(leftChild), move(rightChild),
								var, leftCatMask, rightCatMask
							});
						};

						double leftCoef = (double)popcount64(leftCatMask) / (double)popcount64(repr);
						double rightCoef = 1.0 - leftCoef;

						f(info.var, left, right, leftCoef, rightCoef, dataSplitter, lazySplit);
					}
				},
				[&](const OrdVarPtr& var) {
					uint64_t rightCatMask = repr ^ bottomOne64(repr);
					for(; rightCatMask; rightCatMask ^= bottomOne64(rightCatMask)) {
						uint64_t leftCatMask = repr ^ rightCatMask;

						Cell left = cell;
						info.putRepr(left, leftCatMask);
						Cell right = cell;
						info.putRepr(right, rightCatMask);

						DataSplitter dataSplitter;
						dataSplitter.mask = right.repr & info.mask;

						auto lazySplit = [&](TreePtr leftChild, TreePtr rightChild) {
							return make_unique<Tree>(OrdSplit{
								move(leftChild), move(rightChild),
								var, leftCatMask, rightCatMask
							});
						};

						double leftCoef = (double)popcount64(leftCatMask) / (double)popcount64(repr);
						double rightCoef = 1.0 - leftCoef;

						f(info.var, left, right, leftCoef, rightCoef, dataSplitter, lazySplit);
					}
				}
			);
		}
	}

	template <typename T, typename F>
	void iterateDataSplits(
		Cell cell,
		pair<Cell, T>* data,
		size_t dataCount,
		bool allowEmpty,
		bool adaptCell,
		F f
	) const {
		if(allowEmpty && adaptCell) {
			fail("iterateDataSplits called with both allowEmpty and adaptCell set to true, which is not supported");
		}

		if(!allowEmpty && dataCount <= 1) {
			return;
		}

		iterateSplits(cell, [&](
			const VarPtr& var,
			Cell left, Cell right,
			double leftCoef, double rightCoef,
			DataSplitter dataSplitter,
			auto& lazySplit
		) {
			size_t splitPos = dataSplitter.split(data, dataCount);
			if(allowEmpty || (splitPos != 0 && splitPos != dataCount)) {
				if(adaptCell) {
					left = adaptCellToData_(left, data, splitPos);
					right = adaptCellToData_(right, data + splitPos, dataCount - splitPos);
				}
				f(var, left, right, leftCoef, rightCoef, splitPos, lazySplit);
			}
		});
	}

private:
	struct VarInfo {
		VarPtr var;
		size_t startBit;
		size_t bitCount;
		uint64_t mask;

		uint64_t getRepr(Cell cell) const {
			return (cell.repr & mask) >> startBit;
		}
		void putRepr(Cell& cell, uint64_t repr) const {
			cell.repr ^= (cell.repr ^ (repr << startBit)) & mask;
		}
	};

	template <typename T>
	Cell adaptCellToData_(Cell cell, const pair<Cell, T>* data, size_t dataCount) const {
		uint64_t dataOr = 0;
		for(size_t i = 0; i < dataCount; ++i) {
			dataOr |= data[i].first.repr;
		}
		Cell adapted;
		adapted.repr = (realMask_ & cell.repr) | (catOrdMask_ & dataOr);
		return adapted;
	}

	vector<VarInfo> varInfo_;
	uint64_t realMask_;
	uint64_t catOrdMask_;
	Cell root_;
};

}

namespace std {
template <>
struct hash<pcart::Cell> {
	size_t operator()(pcart::Cell cell) const {
		return hash<uint64_t>()(cell.repr);
	}
};
}
