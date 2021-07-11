#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

namespace pcart {

using std::abs;
using std::ceil;
using std::cin;
using std::cerr;
using std::cout;
using std::exp;
using std::floor;
using std::function;
using std::get;
using std::isfinite;
using std::log;
using std::lower_bound;
using std::make_pair;
using std::make_shared;
using std::make_tuple;
using std::make_unique;
using std::max;
using std::min;
using std::move;
using std::numeric_limits;
using std::ostream;
using std::pair;
using std::round;
using std::shared_ptr;
using std::size_t;
using std::sqrt;
using std::string;
using std::stringstream;
using std::tie;
using std::tuple;
using std::uint64_t;
using std::unique_ptr;
using std::unordered_map;
using std::upper_bound;
using std::variant;
using std::vector;
using std::visit;

const double pi = 4.0 * atan(1.0);

inline void stderrPrint() {
	std::cerr << "\n";
}
template <typename F, typename... T>
void stderrPrint(const F& f, const T&... p) {
	std::cerr << f;
	stderrPrint(p...);
}
template <typename... T>
void fail(const T&... p) {
	stderrPrint("FAIL: ", p...);
	abort();
}

template <typename T>
T fromString(const std::string& str) {
	T ret;
	std::stringstream ss(str);
	ss >> ret;
	if (ss.fail()) fail("fromString: Could not convert string '", str, "' to given type.");
	return ret;
}

template<typename... T> struct LambdaVisitor : T... { using T::operator()...; };
template<typename... T> LambdaVisitor(T...)->LambdaVisitor<T...>;

template <typename Variant, typename... Visitor>
constexpr auto lambdaVisit(Variant&& var, Visitor&&... vis) {
	return visit(LambdaVisitor{ vis... }, var);
}

inline double lognCr(double n, double k) {
	return lgamma(n + 1.0) - lgamma(k + 1.0) - lgamma(n - k + 1.0);
}

inline double addLog(double x, double y) {
	if(x < y) {
		std::swap(x, y);
	}
	return x + log(1.0 + exp(y - x));
}

}
