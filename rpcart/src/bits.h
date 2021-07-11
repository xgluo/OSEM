#pragma once

#include "common.h"

#ifdef _MSC_VER // Only for MSVC
#   include <intrin.h>
#	ifdef _M_IX86
#		pragma intrinsic(_BitScanReverse) // topOne64
#		pragma intrinsic(__popcnt) // popcount64
#	else
#		pragma intrinsic(_BitScanReverse64) // topOne64
#		pragma intrinsic(__popcnt64) // popcount64
#	endif
#	pragma warning (disable : 4146) // bottomOne64
#endif

namespace pcart {

inline constexpr uint64_t ones64(size_t count) {
	return (count == 64)
		? (uint64_t)-1
		: (((uint64_t)1 << count) - (uint64_t)1);
}

inline constexpr uint64_t bit64(size_t pos) {
	return (uint64_t)1 << pos;
}

inline int clz64(uint64_t x) {
#ifdef _MSC_VER
	unsigned long idx;
#	ifdef _M_IX86
	if(_BitScanReverse(&idx, (uint32_t)(x >> 32))) {
		idx += 32;
	} else {
		_BitScanReverse(&idx, (uint32_t)x);
	}
#	else
	_BitScanReverse64(&idx, x);
#	endif
	return 63 - idx;
#else
	static_assert(sizeof(unsigned long long) == sizeof(uint64_t), "The case unsigned long long != uint64_t is not implemented");
	return __builtin_clzll(x);
#endif
}

inline constexpr uint64_t bottomOne64(uint64_t x) {
	return x & -x;
}

inline int popcount64(uint64_t x) {
#ifdef _MSC_VER
#	ifdef _M_IX86
	return (int)__popcnt((uint32_t)(x >> 32)) + (int)__popcnt((uint32_t)x);
#	else
	return (int)__popcnt64(x);
#	endif
#else
	static_assert(sizeof(unsigned long long) == sizeof(uint64_t), "The case unsigned long long != uint64_t is not implemented");
	return __builtin_popcountll(x);
#endif
}

}
