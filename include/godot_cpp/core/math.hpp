/**************************************************************************/
/*  math.hpp                                                              */
/**************************************************************************/
/*                         This file is part of:                          */
/*                             GODOT ENGINE                               */
/*                        https://godotengine.org                         */
/**************************************************************************/
/* Copyright (c) 2014-present Godot Engine contributors (see AUTHORS.md). */
/* Copyright (c) 2007-2014 Juan Linietsky, Ariel Manzur.                  */
/*                                                                        */
/* Permission is hereby granted, free of charge, to any person obtaining  */
/* a copy of this software and associated documentation files (the        */
/* "Software"), to deal in the Software without restriction, including    */
/* without limitation the rights to use, copy, modify, merge, publish,    */
/* distribute, sublicense, and/or sell copies of the Software, and to     */
/* permit persons to whom the Software is furnished to do so, subject to  */
/* the following conditions:                                              */
/*                                                                        */
/* The above copyright notice and this permission notice shall be         */
/* included in all copies or substantial portions of the Software.        */
/*                                                                        */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,        */
/* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF     */
/* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. */
/* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY   */
/* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,   */
/* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE      */
/* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                 */
/**************************************************************************/

#ifndef GODOT_MATH_HPP
#define GODOT_MATH_HPP

#include <godot_cpp/core/defs.hpp>

#include <gdextension_interface.h>

#include <cmath>

namespace godot {

#define Math_SQRT12 0.7071067811865475244008443621048490
#define Math_SQRT2 1.4142135623730950488016887242
#define Math_LN2 0.6931471805599453094172321215
#define Math_PI 3.1415926535897932384626433833
#define Math_TAU 6.2831853071795864769252867666
#define Math_E 2.7182818284590452353602874714
#define Math_INF INFINITY
#define Math_NAN NAN

// Make room for our constexpr's below by overriding potential system-specific macros.
#undef ABS
#undef SIGN
#undef MIN
#undef MAX
#undef CLAMP

// Generic ABS function, for math uses please use Math::abs.
template <typename T>
constexpr T ABS(T m_v) {
	return m_v < 0 ? -m_v : m_v;
}

template <typename T>
constexpr const T SIGN(const T m_v) {
	return m_v == 0 ? 0.0f : (m_v < 0 ? -1.0f : +1.0f);
}

template <typename T, typename T2>
constexpr auto MIN(const T m_a, const T2 m_b) {
	return m_a < m_b ? m_a : m_b;
}

template <typename T, typename T2>
constexpr auto MAX(const T m_a, const T2 m_b) {
	return m_a > m_b ? m_a : m_b;
}

template <typename T, typename T2, typename T3>
constexpr auto CLAMP(const T m_a, const T2 m_min, const T3 m_max) {
	return m_a < m_min ? m_min : (m_a > m_max ? m_max : m_a);
}

// Generic swap template.
#ifndef SWAP
#define SWAP(m_x, m_y) __swap_tmpl((m_x), (m_y))
template <typename T>
inline void __swap_tmpl(T &x, T &y) {
	T aux = x;
	x = y;
	y = aux;
}
#endif // SWAP

/* Functions to handle powers of 2 and shifting. */

// Function to find the next power of 2 to an integer.
static _FORCE_INLINE_ unsigned int next_power_of_2(unsigned int x) {
	if (x == 0) {
		return 0;
	}

	--x;
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;

	return ++x;
}

// Function to find the previous power of 2 to an integer.
static _FORCE_INLINE_ unsigned int previous_power_of_2(unsigned int x) {
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	return x - (x >> 1);
}

// Function to find the closest power of 2 to an integer.
static _FORCE_INLINE_ unsigned int closest_power_of_2(unsigned int x) {
	unsigned int nx = next_power_of_2(x);
	unsigned int px = previous_power_of_2(x);
	return (nx - x) > (x - px) ? px : nx;
}

// Get a shift value from a power of 2.
static inline int get_shift_from_power_of_2(unsigned int p_bits) {
	for (unsigned int i = 0; i < 32; i++) {
		if (p_bits == (unsigned int)(1 << i)) {
			return i;
		}
	}

	return -1;
}

template <typename T>
static _FORCE_INLINE_ T nearest_power_of_2_templated(T x) {
	--x;

	// The number of operations on x is the base two logarithm
	// of the number of bits in the type. Add three to account
	// for sizeof(T) being in bytes.
	size_t num = get_shift_from_power_of_2(sizeof(T)) + 3;

	// If the compiler is smart, it unrolls this loop.
	// If it's dumb, this is a bit slow.
	for (size_t i = 0; i < num; i++) {
		x |= x >> (1 << i);
	}

	return ++x;
}

// Function to find the nearest (bigger) power of 2 to an integer.
static inline unsigned int nearest_shift(unsigned int p_number) {
	for (int i = 30; i >= 0; i--) {
		if (p_number & (1 << i)) {
			return i + 1;
		}
	}

	return 0;
}

// constexpr function to find the floored log2 of a number
template <typename T>
constexpr T floor_log2(T x) {
	return x < 2 ? x : 1 + floor_log2(x >> 1);
}

// Get the number of bits needed to represent the number.
// IE, if you pass in 8, you will get 4.
// If you want to know how many bits are needed to store 8 values however, pass in (8 - 1).
template <typename T>
constexpr T get_num_bits(T x) {
	return floor_log2(x);
}

// Swap 16, 32 and 64 bits value for endianness.
#if defined(__GNUC__)
#define BSWAP16(x) __builtin_bswap16(x)
#define BSWAP32(x) __builtin_bswap32(x)
#define BSWAP64(x) __builtin_bswap64(x)
#else
static inline uint16_t BSWAP16(uint16_t x) {
	return (x >> 8) | (x << 8);
}

static inline uint32_t BSWAP32(uint32_t x) {
	return ((x << 24) | ((x << 8) & 0x00FF0000) | ((x >> 8) & 0x0000FF00) | (x >> 24));
}

static inline uint64_t BSWAP64(uint64_t x) {
	x = (x & 0x00000000FFFFFFFF) << 32 | (x & 0xFFFFFFFF00000000) >> 32;
	x = (x & 0x0000FFFF0000FFFF) << 16 | (x & 0xFFFF0000FFFF0000) >> 16;
	x = (x & 0x00FF00FF00FF00FF) << 8 | (x & 0xFF00FF00FF00FF00) >> 8;
	return x;
}
#endif

namespace Math {

// This epsilon should match the one used by Godot for consistency.
// Using `f` when `real_t` is float.
#define CMP_EPSILON 0.00001f
#define CMP_EPSILON2 (CMP_EPSILON * CMP_EPSILON)

// This epsilon is for values related to a unit size (scalar or vector len).
#ifdef PRECISE_MATH_CHECKS
#define UNIT_EPSILON 0.00001
#else
// Tolerate some more floating point error normally.
#define UNIT_EPSILON 0.001
#endif

// Functions reproduced as in Godot's source code `math_funcs.h`.
// Some are overloads to automatically support changing real_t into either double or float in the way Godot does.

inline double fmod(double p_x, double p_y) {
	return ::fmod(p_x, p_y);
}
inline float fmod(float p_x, float p_y) {
	return ::fmodf(p_x, p_y);
}

inline double fposmod(double p_x, double p_y) {
	double value = Math::fmod(p_x, p_y);
	if ((value < 0 && p_y > 0) || (value > 0 && p_y < 0)) {
		value += p_y;
	}
	value += 0.0;
	return value;
}
inline float fposmod(float p_x, float p_y) {
	float value = Math::fmod(p_x, p_y);
	if ((value < 0 && p_y > 0) || (value > 0 && p_y < 0)) {
		value += p_y;
	}
	value += 0.0f;
	return value;
}

inline float fposmodp(float p_x, float p_y) {
	float value = Math::fmod(p_x, p_y);
	if (value < 0) {
		value += p_y;
	}
	value += 0.0f;
	return value;
}
inline double fposmodp(double p_x, double p_y) {
	double value = Math::fmod(p_x, p_y);
	if (value < 0) {
		value += p_y;
	}
	value += 0.0;
	return value;
}

inline int64_t posmod(int64_t p_x, int64_t p_y) {
	int64_t value = p_x % p_y;
	if ((value < 0 && p_y > 0) || (value > 0 && p_y < 0)) {
		value += p_y;
	}
	return value;
}

/*template trig func*/
template<typename T, uint8_t flag>
inline T t_sine90(const T x) {
 const	T x_squared = -x * x;
	T accumulation = x;
	T iteration = x;
	switch(flag) {
	 case 1:
                            //  1 / (3 * 2)   
	iteration *= x_squared * 0.1666666666666666574148081281236954964697;
	accumulation += iteration;
                            //  1 / (5 * 4)   
	iteration *= x_squared * 0.0500000000000000027755575615628913510591;
	accumulation += iteration;
                            //  1 / (7 * 6)   
	iteration *= x_squared * 0.0238095238095238082021154468748136423528;
	accumulation += iteration;
                            //  1 / (9 * 8)   
	iteration *= x_squared * 0.0138888888888888881179006773436412913725;
	accumulation += iteration;
                            //  1 / (11 * 10)   
	iteration *= x_squared * 0.0090909090909090904675249333877218305133;
	accumulation += iteration;
	 break;
	 case 2:
                            //  1 / (3 * 2)   
	iteration *= x_squared * 0.1666666666666666574148081281236954964697;
	accumulation += iteration;
                            //  1 / (5 * 4)   
	iteration *= x_squared * 0.0500000000000000027755575615628913510591;
	accumulation += iteration;
	 break;
	 default:
                            //  1 / (3 * 2)   
	iteration *= x_squared * 0.1666666666666666574148081281236954964697;
	accumulation += iteration;
                            //  1 / (5 * 4)   
	iteration *= x_squared * 0.0500000000000000027755575615628913510591;
	accumulation += iteration;
                            //  1 / (7 * 6)   
	iteration *= x_squared * 0.0238095238095238082021154468748136423528;
	accumulation += iteration;
                            //  1 / (9 * 8)   
	iteration *= x_squared * 0.0138888888888888881179006773436412913725;
	accumulation += iteration;
                            //  1 / (11 * 10)   
	iteration *= x_squared * 0.0090909090909090904675249333877218305133;
	accumulation += iteration;
                            //  1 / (13 * 12)   
	iteration *= x_squared * 0.0064102564102564100340098107722042186651;
	accumulation += iteration;
                            //  1 / (15 * 14)   
	iteration *= x_squared * 0.0047619047619047623343124797656855662353;
	accumulation += iteration;
                            //  1 / (17 * 16)   
	iteration *= x_squared * 0.0036764705882352940666257801183292031055;
	accumulation += iteration;
                            //  1 / (19 * 18)   
	iteration *= x_squared * 0.0029239766081871343406106689144507981837;
	accumulation += iteration;
	}
	
 return accumulation;
}


template<typename T, uint8_t flag>
inline T t_cosine90(const T x) {
 const	T x_squared = -x * x;
	T accumulation = 1.0;
	T iteration = 1.0;
	switch(flag) {
	 case 1:
                            //  1 / (2 * 1)   
	iteration *= x_squared * 0.5000000000000000000000000000000000000000;
	accumulation += iteration;
                            //  1 / (4 * 3)   
	iteration *= x_squared * 0.0833333333333333333333333333333333293210;
	accumulation += iteration;
                            //  1 / (6 * 5)   
	iteration *= x_squared * 0.0333333333333333333333333333333333329321;
	accumulation += iteration;
                            //  1 / (8 * 7)   
	iteration *= x_squared * 0.0178571428571428571428571428571428562831;
	accumulation += iteration;
                            //  1 / (10 * 9)   
	iteration *= x_squared * 0.0111111111111111111111111111111111114789;
	accumulation += iteration;
	 break;
	 case 2:
                            //  1 / (2 * 1)   
	iteration *= x_squared * 0.5000000000000000000000000000000000000000;
	accumulation += iteration;
                            //  1 / (4 * 3)   
	iteration *= x_squared * 0.0833333333333333333333333333333333293210;
	accumulation += iteration;
	 break;
	 default:
                            //  1 / (2 * 1)   
	iteration *= x_squared * 0.5000000000000000000000000000000000000000;
	accumulation += iteration;
                            //  1 / (4 * 3)   
	iteration *= x_squared * 0.0833333333333333333333333333333333293210;
	accumulation += iteration;
                            //  1 / (6 * 5)   
	iteration *= x_squared * 0.0333333333333333333333333333333333329321;
	accumulation += iteration;
                            //  1 / (8 * 7)   
	iteration *= x_squared * 0.0178571428571428571428571428571428562831;
	accumulation += iteration;
                            //  1 / (10 * 9)   
	iteration *= x_squared * 0.0111111111111111111111111111111111114789;
	accumulation += iteration;
                            //  1 / (12 * 11)   
	iteration *= x_squared * 0.0075757575757575757575757575757575759400;
	accumulation += iteration;
                            //  1 / (14 * 13)   
	iteration *= x_squared * 0.0054945054945054945054945054945054948775;
	accumulation += iteration;
                            //  1 / (16 * 15)   
	iteration *= x_squared * 0.0041666666666666666666666666666666666165;
	accumulation += iteration;
                            //  1 / (18 * 17)   
	iteration *= x_squared * 0.0032679738562091503267973856209150326773;
	accumulation += iteration;
	}
 return accumulation;
}


//approximating longest curve in arc tangent from atan(tan(66°)) towards pi/2
//using linear interpolation and bezier curve
template<typename T>
inline T __atan66(const T x) noexcept {
/*
	const double _tan_89 = 57.28996163075914438423;
	
	const double _tan_66 = 2.24603677390421641036;
	const double _atan_66 = 1.15191730631625754988;
	
	const double _tan_76 = 4.01078093353584552716;
	const double _atan_76 = 1.32645023151569052544;
	
	const double _tan_83 = 8.14434642797459318331;
	const double _atan_83 = 1.44862327915529354172;
*/
	constexpr T _pi_half = Math_PI/2;
  // x >= _tan_66 && x < _tan_76
	if(x >= 2.24603677390421641036 && x < 4.01078093353584552716) {
		const T t = (x-2.24603677390421641036) / (57.28996163075914438423-2.24603677390421641036) * 0.00001;
		const T mp = ((1.15191730631625754988+_pi_half) * 0.5) + 299000.0;
		const T ma = 1.15191730631625754988 + t * (mp-1.15191730631625754988);
		const T mb = mp + t * (_pi_half-mp);
		return ma + t * (mb-ma);
		//x >= _tan76 && x < _tan_83
	} else if(x >= 4.01078093353584552716 && x < 8.14434642797459318331) {
		
		const T t = (x-4.01078093353584552716) / (57.28996163075914438423-4.01078093353584552716) * 0.00001;
		const T mp = ((1.32645023151569052544+_pi_half) * 0.5) + 100000.0;
		const T ma = 1.32645023151569052544 + t * (mp-1.32645023151569052544);
		const T mb = mp + t * (_pi_half-mp);
		return ma + t * (mb-ma);
		//x >= _tan_83 && x < _tan_89
	} else if(x >= 8.14434642797459318331 && x < 57.28996163075914438423) {
		const T t = (x-8.14434642797459318331) / (57.28996163075914438423-8.14434642797459318331) * 0.00001;
		const T mp = ((1.44862327915529354172+_pi_half) * 0.5) + 10000.0;
		const T ma = 1.44862327915529354172 + t * (mp-1.44862327915529354172);
		const T mb = mp + t * (_pi_half-mp);
		return ma + t * (mb-ma);
		//x >= _tan_89 && x < _tan_90
	} else if(x >= 57.28996163075914438423 && x < 16331239353195370.0) {
		const T t = (x-57.28996163075914438423) / (16331239353195370.0-57.28996163075914438423) * 0.00001;
		const T mp = ((1.55334303427495323824+_pi_half) * 0.5) + 10000.0;
		const T ma = 1.55334303427495323824 + t * (mp-1.55334303427495323824);
		const T mb = mp + t * (_pi_half-mp);
		return ma + t * (mb-ma);
	} else 
	 return _pi_half;
	return 0;
}


template<typename T, uint8_t flag>
inline T t_sin(const T theta) {
	constexpr T _pi = Math_PI;
	constexpr T _pi_half = Math_PI/2;
	constexpr T _pi_270 = Math_PI+_pi_half;
	constexpr T _2_pi = Math_PI*2;
	constexpr T _inv_2_pi = 1.0 / _2_pi;
	
	const T x = ((theta) - (int)((theta) * _inv_2_pi) * (_2_pi));

 if(x <= _pi_half) //0 - 90 deg
  return t_sine90<T, flag>(x);
 else if(x > _pi_half && x <= _pi) //90 - 180 deg
  return t_sine90<T, flag>(_pi_half-(x-_pi_half));
 else if(x > _pi && x <= _pi_270) //180 - 270 deg
  return -t_sine90<T, flag>(x-_pi);
 else //270 - 360 deg
  return -t_sine90<T, flag>(_pi_half-(x-_pi_270));
 return 0;
}


template<typename T, uint8_t flag>
inline T t_cos(const T theta) {
	constexpr T _pi = Math_PI;
	constexpr T _pi_half = Math_PI/2;
	constexpr T _pi_270 = Math_PI+_pi_half;
	constexpr T _2_pi = Math_PI*2;
	constexpr T _inv_2_pi = 1.0 / _2_pi;

	const T x = ((theta) - (int)((theta) * _inv_2_pi) * (_2_pi));

 if(x <= _pi_half) //0 - 90 deg
  return  t_cosine90<T, flag>(x);
 else if(x > _pi_half && x <= _pi) //90 - 180 deg
  return -t_cosine90<T, flag>(_pi_half-(x-_pi_half));
 else if(x > _pi && x <= _pi_270) //180 - 270 deg
  return -t_cosine90<T, flag>(x-_pi);
 else //270 - 360 deg
  return t_cosine90<T, flag>(_pi_half-(x-_pi_270));
 return 0;
}


template<typename T, uint8_t flag>
inline T t_tan(const T theta) {
	constexpr T _pi = Math_PI;
	constexpr T _pi_half = Math_PI/2;
	constexpr T _pi_270 = Math_PI+_pi_half;
	constexpr T _2_pi = Math_PI*2;
	constexpr T _inv_2_pi = 1.0 / _2_pi;
	
	const T x = ((theta) - (int)((theta) * _inv_2_pi) * (_2_pi));

 const T x_sub_pi_half = x-_pi_half;
 const T x_sub_pi = x-_pi;
 const T x_sub_pi_270 = x-_pi_270;
 
 if(x <= _pi_half) //0 - 90 deg
  return t_sine90<T, flag>(x) / t_cosine90<T, flag>(x);
 else if(x > _pi_half && x <= _pi) //90 - 180 deg
  return t_sine90<T, flag>(_pi_half-x_sub_pi_half) / -t_cosine90<T, flag>(_pi_half-x_sub_pi_half);
 else if(x > _pi && x <= _pi_270) //180 - 270 deg
  return -t_sine90<T, flag>(x_sub_pi) / -t_cosine90<T, flag>(x_sub_pi);
 else //270 - 360 deg
  return -t_sine90<T, flag>(_pi_half-x_sub_pi_270) / t_cosine90<T, flag>(_pi_half-x_sub_pi_270);
 return 0;
}

template<typename T, uint8_t flag>
inline T t_asin(const T ix) {
 constexpr T _pi = Math_PI;
	constexpr T _pi_half = Math_PI/2;
 constexpr T _inv_pi = 1.0 / _pi;
 T x = ((ix) - (int)((ix) * _inv_pi) * (_pi)); //-180 - 180 deg
 
 if(x > _pi_half)
  x -= _pi_half;
 else if (x <= -_pi_half)
  x += _pi_half;

 const	T x_squared = x * x;
	T accumulation = x;
	T iteration = x;
	T final_iter;
	
 switch(flag) {
 	case 1:
	iteration *= x_squared;
	              //(1/2)                 1 / 3   
	accumulation += 0.5 * (iteration * 0.3333333333333333333333333333333333172839);

	iteration *= x_squared;
	            //((1*3)/(2*4))              1 / 5   
	accumulation += 0.375 * (iteration * 0.2000000000000000000000000000000000096296);

	iteration *= x_squared;
	          //((1*3*5)/(2*4*6))             1 / 7   
	accumulation += 0.3125 * (iteration * 0.1428571428571428571428571428571428502645);

	iteration *= x_squared;
	           //((1*3*5*7)/(2*4*6*8))         1 / 9   
	accumulation += 0.2734375 * (iteration * 0.1111111111111111111111111111111111057613);

	iteration *= x_squared;  
	           //((1*3*5*7*9)/(2*4*6*8*10))         1 / 11   
	accumulation += 0.24609375 * (iteration * 0.0909090909090909090909090909090909112795);

	iteration *= x_squared;
	    //((1*3*5*7*9*11)/(2*4*6*8*10*12))          1 / 13   
	accumulation += 0.2255859375 * (iteration * 0.0769230769230769230769230769230769267806);

	iteration *= x_squared;
	//((1*3*5*7*9*11*13)/(2*4*6*8*10*12*14))               1 / 15   
	accumulation += 0.20947265625 * (iteration * 0.0666666666666666666666666666666666658642);

	iteration *= x_squared;
	//((1*3*5*7*9*11*13*15)/(2*4*6*8*10*12*14*16))          1 / 17   
	accumulation += 0.196380615234375 * (iteration * 0.0588235294117647058823529411764705889434);

 final_iter = accumulation;
 for(char i = 0; i < 6; i++)
	 final_iter = final_iter - (t_sine90<T, flag>(final_iter) - x) / t_cosine90<T, flag>(final_iter);
 	break;
 	case 2:
	iteration *= x_squared;
	              //(1/2)                 1 / 3   
	accumulation += 0.5 * (iteration * 0.3333333333333333333333333333333333172839);

	iteration *= x_squared;
	            //((1*3)/(2*4))              1 / 5   
	accumulation += 0.375 * (iteration * 0.2000000000000000000000000000000000096296);

	iteration *= x_squared;
	          //((1*3*5)/(2*4*6))             1 / 7   
	accumulation += 0.3125 * (iteration * 0.1428571428571428571428571428571428502645);

	iteration *= x_squared;
	           //((1*3*5*7)/(2*4*6*8))         1 / 9   
	accumulation += 0.2734375 * (iteration * 0.1111111111111111111111111111111111057613);
 	
 final_iter = accumulation;
 for(char i = 0; i < 2; i++)
	 final_iter = final_iter - (t_sine90<T, flag>(final_iter) - x) / t_cosine90<T, flag>(final_iter);

 	break;
 	default:
	iteration *= x_squared;
	              //(1/2)                 1 / 3   
	accumulation += 0.5 * (iteration * 0.3333333333333333333333333333333333172839);

	iteration *= x_squared;
	            //((1*3)/(2*4))              1 / 5   
	accumulation += 0.375 * (iteration * 0.2000000000000000000000000000000000096296);

	iteration *= x_squared;
	          //((1*3*5)/(2*4*6))             1 / 7   
	accumulation += 0.3125 * (iteration * 0.1428571428571428571428571428571428502645);

	iteration *= x_squared;
	           //((1*3*5*7)/(2*4*6*8))         1 / 9   
	accumulation += 0.2734375 * (iteration * 0.1111111111111111111111111111111111057613);

	iteration *= x_squared;  
	           //((1*3*5*7*9)/(2*4*6*8*10))         1 / 11   
	accumulation += 0.24609375 * (iteration * 0.0909090909090909090909090909090909112795);

	iteration *= x_squared;
	    //((1*3*5*7*9*11)/(2*4*6*8*10*12))          1 / 13   
	accumulation += 0.2255859375 * (iteration * 0.0769230769230769230769230769230769267806);

	iteration *= x_squared;
	//((1*3*5*7*9*11*13)/(2*4*6*8*10*12*14))               1 / 15   
	accumulation += 0.20947265625 * (iteration * 0.0666666666666666666666666666666666658642);

	iteration *= x_squared;
	//((1*3*5*7*9*11*13*15)/(2*4*6*8*10*12*14*16))          1 / 17   
	accumulation += 0.196380615234375 * (iteration * 0.0588235294117647058823529411764705889434);

	iteration *= x_squared;
	//((1*3*5*7*9*11*13*15*17)/(2*4*6*8*10*12*14*16*18))            1 / 19   
	accumulation += 0.1854705810546875 * (iteration * 0.0526315789473684210526315789473684213694);

 iteration *= x_squared;
 //((1*3*5*7*9*11*13*15*17*19)/(2*4*6*8*10*12*14*16*18*20))          1 / 21   
 accumulation += 0.176197052001953125 * (iteration * 0.0476190476190476190476190476190476167548);

 iteration *= x_squared;
 //((1*3*5*7*9*11*13*15*17*19*21)/(2*4*6*8*10*12*14*16*18*20*22))               1 / 23   
 accumulation += 0.1681880950927734375 * (iteration * 0.0434782608695652173913043478260869591385);
 
 iteration *= x_squared;
 //((1*3*5*7*9*11*13*15*17*19*21*23)/(2*4*6*8*10*12*14*16*18*20*22*24))           1 / 25   
 accumulation += 0.1611802577972412109375 * (iteration * 0.0400000000000000000000000000000000007222);

 iteration *= x_squared;
 //((1*3*5*7*9*11*13*15*17*19*21*23*25)/(2*4*6*8*10*12*14*16*18*20*22*24*26))             1 / 27   
 accumulation += 0.15498101711273193359375 * (iteration * 0.0370370370370370370370370370370370372599);

 iteration *= x_squared;                                    //  1 / 29   
 accumulation += 0.1494459807872772216796875 * (iteration * 0.0344827586206896551724137931034482752395);

 iteration *= x_squared;                                     //  1 / 31   
 accumulation += 0.14446444809436798095703125 * (iteration * 0.0322580645161290322580645161290322572879);

 iteration *= x_squared;                                          //  1 / 33   
 accumulation += 0.1399499340914189815521240234375 * (iteration * 0.0303030303030303030303030303030303037598);
 
 final_iter = accumulation;
 for(char i = 0; i < 8; i++)
	 final_iter = final_iter - (t_sine90<T, flag>(final_iter) - x) / t_cosine90<T, flag>(final_iter);

 }
	
	return final_iter;
}


template<typename T, uint8_t flag>
inline T t_acos(const T ix) {
	constexpr T _pi_half = Math_PI/2;
	constexpr T _inv_pi_half = 1.0 / _pi_half;
 const T x = ((ix) - (int)((ix) * _inv_pi_half) * (_pi_half)); //-90 - 90 deg

 const T x_squared = x * x;
 T iteration = x;
 T accumulation = _pi_half;
 T final_iter;
 
 //floating point error correction between cos(1°) - cos(41°)
 if(x <= 0.99984769515639126958 && x >= 0.75470958022277201405)
  return _pi_half - t_asin<T, flag>(x);
 
 switch(flag) {
  case 1:
 iteration *= x_squared;
                //(2 - 1) * (2 - 3)
 accumulation -= -1.0 * iteration * 0.5000000000000000000000000000000000000000; //1 / 2

 iteration *= x_squared;
                //(2 * 2 - 1) * (2 * 2 - 3)
 accumulation -= 3.0 * iteration * 0.2500000000000000000000000000000000000000; //1 / (2 * 2)

 iteration *= x_squared;
                //(2 * 3 - 1) * (2 * 3 - 3)
 accumulation -= 15.0 * iteration * 0.1666666666666666666666666666666666586420; //1 / (2 * 3)
 
 iteration *= x_squared;
                //(2 * 4 - 1) * (2 * 4 - 3)
 accumulation -= 35.0 * iteration * 0.1250000000000000000000000000000000000000; //1 / (2 * 4)
 
 iteration *= x_squared;
                //(2 * 5 - 1) * (2 * 5 - 3)
 accumulation -= 63.0 * iteration * 0.1000000000000000000000000000000000048148; //1 / (2 * 5)
 
 iteration *= x_squared;
                //(2 * 6 - 1) * (2 * 6 - 3)
 accumulation -= 99.0 * iteration * 0.0833333333333333333333333333333333293210; //1 / (2 * 6)
 
 iteration *= x_squared;
                //(2 * 7 - 1) * (2 * 7 - 3)
 accumulation -= 143.0 * iteration * 0.0714285714285714285714285714285714251323; //1 / (2 * 7)
 
 iteration *= x_squared;
                //(2 * 8 - 1) * (2 * 8 - 3)
 accumulation -= 195.0 * iteration * 0.0625000000000000000000000000000000000000; //1 / (2 * 8)
 
 final_iter = accumulation;
	for(char i = 0; i < 6; i++)
	 final_iter = final_iter - (t_cosine90<T, flag>(final_iter) - x) / -t_sine90<T, flag>(final_iter);
  break;
  case 2:
 iteration *= x_squared;
                //(2 - 1) * (2 - 3)
 accumulation -= -1.0 * iteration * 0.5000000000000000000000000000000000000000; //1 / 2

 iteration *= x_squared;
                //(2 * 2 - 1) * (2 * 2 - 3)
 accumulation -= 3.0 * iteration * 0.2500000000000000000000000000000000000000; //1 / (2 * 2)

 iteration *= x_squared;
                //(2 * 3 - 1) * (2 * 3 - 3)
 accumulation -= 15.0 * iteration * 0.1666666666666666666666666666666666586420; //1 / (2 * 3)
 
 iteration *= x_squared;
                //(2 * 4 - 1) * (2 * 4 - 3)
 accumulation -= 35.0 * iteration * 0.1250000000000000000000000000000000000000; //1 / (2 * 4)
 
 final_iter = accumulation;
	for(char i = 0; i < 2; i++)
	 final_iter = final_iter - (t_cosine90<T, flag>(final_iter) - x) / -t_sine90<T, flag>(final_iter);
  break;
  default:
 iteration *= x_squared;
                //(2 - 1) * (2 - 3)
 accumulation -= -1.0 * iteration * 0.5000000000000000000000000000000000000000; //1 / 2

 iteration *= x_squared;
                //(2 * 2 - 1) * (2 * 2 - 3)
 accumulation -= 3.0 * iteration * 0.2500000000000000000000000000000000000000; //1 / (2 * 2)

 iteration *= x_squared;
                //(2 * 3 - 1) * (2 * 3 - 3)
 accumulation -= 15.0 * iteration * 0.1666666666666666666666666666666666586420; //1 / (2 * 3)
 
 iteration *= x_squared;
                //(2 * 4 - 1) * (2 * 4 - 3)
 accumulation -= 35.0 * iteration * 0.1250000000000000000000000000000000000000; //1 / (2 * 4)
 
 iteration *= x_squared;
                //(2 * 5 - 1) * (2 * 5 - 3)
 accumulation -= 63.0 * iteration * 0.1000000000000000000000000000000000048148; //1 / (2 * 5)
 
 iteration *= x_squared;
                //(2 * 6 - 1) * (2 * 6 - 3)
 accumulation -= 99.0 * iteration * 0.0833333333333333333333333333333333293210; //1 / (2 * 6)
 
 iteration *= x_squared;
                //(2 * 7 - 1) * (2 * 7 - 3)
 accumulation -= 143.0 * iteration * 0.0714285714285714285714285714285714251323; //1 / (2 * 7)
 
 iteration *= x_squared;
                //(2 * 8 - 1) * (2 * 8 - 3)
 accumulation -= 195.0 * iteration * 0.0625000000000000000000000000000000000000; //1 / (2 * 8)
 
 iteration *= x_squared;
                //(2 * 9 - 1) * (2 * 9 - 3)
 accumulation -= 255.0 * iteration * 0.0555555555555555555555555555555555528807; //1 / (2 * 9)
 
 iteration *= x_squared;
                //(2 * 10 - 1) * (2 * 10 - 3)
 accumulation -= 323.0 * iteration * 0.0500000000000000000000000000000000024074; //1 / (2 * 10)
 
 iteration *= x_squared;
                //(2 * 11 - 1) * (2 * 11 - 3)
 accumulation -= 399.0 * iteration * 0.0454545454545454545454545454545454556397; //1 / (2 * 11)
 
 iteration *= x_squared;
                //(2 * 12 - 1) * (2 * 12 - 3)
 accumulation -= 483.0 * iteration * 0.0416666666666666666666666666666666646605; //1 / (2 * 12)
 
 iteration *= x_squared;
                //(2 * 13 - 1) * (2 * 13 - 3)
 accumulation -= 575.0 * iteration * 0.0384615384615384615384615384615384633903; //1 / (2 * 13)
 
 iteration *= x_squared;
                //(2 * 14 - 1) * (2 * 14 - 3)
 accumulation -= 675.0 * iteration * 0.0357142857142857142857142857142857125661; //1 / (2 * 14)
 
 iteration *= x_squared;
                //(2 * 15 - 1) * (2 * 15 - 3)
 accumulation -= 783.0 * iteration * 0.0333333333333333333333333333333333329321; //1 / (2 * 15)
 
 iteration *= x_squared;
                //(2 * 16 - 1) * (2 * 16 - 3)
 accumulation -= 899.0 * iteration * 0.0312500000000000000000000000000000000000; //1 / (2 * 16)

 final_iter = accumulation;
	for(char i = 0; i < 8; i++)
	 final_iter = final_iter - (t_cosine90<T, flag>(final_iter) - x) / -t_sine90<T, flag>(final_iter);
 }

 return final_iter;
}


template<typename T, uint8_t flag>
inline T t_atan(const T ix) {
	constexpr T _pi = Math_PI;
	constexpr T _pi_half = Math_PI/2;
	constexpr T _inv_pi_half = 1.0 / _pi_half;
 T x = ix;

//floating point error correction when x is between tan(66°) - tan(89°)
 if(x >= 2.24603677390421641036)
  return __atan66<T>(x);
 else if(x <= -2.24603677390421641036)
  return -__atan66<T>(-x);

 	const T x_squared = -x * x;
	T accumulation = x;
	T iteration = x;
 T tangent;
 T final_iter;
 
 //floating point correction between tan(67°) - tan(89°)
 if(x >= 2.35585236582375312508 && x <= 57.28996163075914438423) {
                   //tan(67°)                 tan(89°)
  const T p = (x - 2.35585236582375312508) / (57.28996163075914438423 - 2.35585236582375312508);
  const T max = 1.55334303427495323824; //atan(tan(89°))
  const T min = 1.16937059883620086964; //atan(tan(67°))
  return min + p * (max-min); //lerp
 }
 
	switch(flag) {
	 case 1:
                           //  1 / (3 * 2)   
 iteration *= x_squared * 0.1666666666666666574148081281236954964697;
 accumulation += iteration;
                           //  1 / (5 * 4)   
 iteration *= x_squared * 0.0500000000000000027755575615628913510591;
 accumulation += iteration;
                           //  1 / (7 * 6)   
 iteration *= x_squared * 0.0238095238095238082021154468748136423528;
 accumulation += iteration;
                           //  1 / (9 * 8)   
 iteration *= x_squared * 0.0138888888888888881179006773436412913725;
 accumulation += iteration;
                           //  1 / (11 * 10)   
 iteration *= x_squared * 0.0090909090909090904675249333877218305133;
 accumulation += iteration;
                           //  1 / (13 * 12)   
 iteration *= x_squared * 0.0064102564102564100340098107722042186651;
 accumulation += iteration;
                           //  1 / (15 * 14)   
 iteration *= x_squared * 0.0047619047619047623343124797656855662353;
 accumulation += iteration;
                           //  1 / (17 * 16)   
 iteration *= x_squared * 0.0036764705882352940666257801183292031055;
 accumulation += iteration;
 
 final_iter = accumulation;
 for(char i = 0; i < 6; i++) {
  tangent = t_tan<T, flag>(final_iter);
  final_iter = final_iter - (tangent - x) / (1.0 + tangent * tangent);
 }
	 break;
	 case 2:
                           //  1 / (3 * 2)   
 iteration *= x_squared * 0.1666666666666666574148081281236954964697;
 accumulation += iteration;
                           //  1 / (5 * 4)   
 iteration *= x_squared * 0.0500000000000000027755575615628913510591;
 accumulation += iteration;
                           //  1 / (7 * 6)   
 iteration *= x_squared * 0.0238095238095238082021154468748136423528;
 accumulation += iteration;
                           //  1 / (9 * 8)   
 iteration *= x_squared * 0.0138888888888888881179006773436412913725;
 accumulation += iteration;
 
 final_iter = accumulation;
 for(char i = 0; i < 2; i++) {
  tangent = t_tan<T, flag>(final_iter);
  final_iter = final_iter - (tangent - x) / (1.0 + tangent * tangent);
 }
	 break;
	 default:
                           //  1 / (3 * 2)   
 iteration *= x_squared * 0.1666666666666666574148081281236954964697;
 accumulation += iteration;
                           //  1 / (5 * 4)   
 iteration *= x_squared * 0.0500000000000000027755575615628913510591;
 accumulation += iteration;
                           //  1 / (7 * 6)   
 iteration *= x_squared * 0.0238095238095238082021154468748136423528;
 accumulation += iteration;
                           //  1 / (9 * 8)   
 iteration *= x_squared * 0.0138888888888888881179006773436412913725;
 accumulation += iteration;
                           //  1 / (11 * 10)   
 iteration *= x_squared * 0.0090909090909090904675249333877218305133;
 accumulation += iteration;
                           //  1 / (13 * 12)   
 iteration *= x_squared * 0.0064102564102564100340098107722042186651;
 accumulation += iteration;
                           //  1 / (15 * 14)   
 iteration *= x_squared * 0.0047619047619047623343124797656855662353;
 accumulation += iteration;
                           //  1 / (17 * 16)   
 iteration *= x_squared * 0.0036764705882352940666257801183292031055;
 accumulation += iteration;
                           //  1 / (19 * 18)   
 iteration *= x_squared * 0.0029239766081871343406106689144507981837;
 accumulation += iteration;
                           //  1 / (21 * 20)   
 iteration *= x_squared * 0.0023809523809523809523809523809523811387;
 accumulation += iteration;
                           //  1 / (23 * 22)   
 iteration *= x_squared * 0.0019762845849802371541501976284584980059;
 accumulation += iteration;
                           //  1 / (25 * 24)   
 iteration *= x_squared * 0.0016666666666666666666666666666666667595;
 accumulation += iteration;
                           //  1 / (27 * 26)   
 iteration *= x_squared * 0.0014245014245014245014245014245014244377;
 accumulation += iteration;
                           //  1 / (29 * 28)   
 iteration *= x_squared * 0.0012315270935960591133004926108374384551;
 accumulation += iteration;
                           //  1 / (31 * 30)   
 iteration *= x_squared * 0.0010752688172043010752688172043010752680;
 accumulation += iteration;
 
 final_iter = accumulation;
 for(char i = 0; i < 8; i++) {
  tangent = t_tan<T, flag>(final_iter);
  final_iter = final_iter - (tangent - x) / (1.0 + tangent * tangent);
 }
	}

	return final_iter;
}


template<typename T, uint8_t flag>
inline T t_atan2(const T y, const T x) {
	constexpr T _pi = Math_PI;
	constexpr T _pi_half = Math_PI/2;
 if(x > 0)
  return t_atan<T, flag>(y / x);
 else if(x < 0) {
   if(y >= 0)
     return t_atan<T, flag>(y / x) + _pi;
    else
     return t_atan<T, flag>(y / x) - _pi;
  } else {
   if(y > 0)
    return _pi_half;
   else if(y < 0)
    return -_pi_half;
  }
 return 0;
}
/*end template func*/

inline double floor(double p_x) {
	return ::floor(p_x);
}
inline float floor(float p_x) {
	return ::floorf(p_x);
}

inline double ceil(double p_x) {
	return ::ceil(p_x);
}
inline float ceil(float p_x) {
	return ::ceilf(p_x);
}

inline double exp(double p_x) {
	return ::exp(p_x);
}
inline float exp(float p_x) {
	return ::expf(p_x);
}

inline double sin(double p_x) {
#if __TRIG_FUNC >= 0
 return t_sin<double, __TRIG_FUNC>(p_x);
#else 
	return ::sin(p_x);
#endif
}
inline float sin(float p_x) {
#if __TRIG_FUNC >= 0
 return t_sin<float, __TRIG_FUNC>(p_x);
#else
	return ::sinf(p_x);
#endif
}

inline double cos(double p_x) {
#if __TRIG_FUNC >= 0
 return t_cos<double, __TRIG_FUNC>(p_x);
#else 
	return ::cos(p_x);
#endif
}
inline float cos(float p_x) {
#if __TRIG_FUNC >= 0
 return t_cos<float, __TRIG_FUNC>(p_x);
#else
	return ::cosf(p_x);
#endif
}

inline double tan(double p_x) {
#if __TRIG_FUNC >= 0
 return t_tan<double, __TRIG_FUNC>(p_x);
#else
	return ::tan(p_x);
#endif
}
inline float tan(float p_x) {
#if __TRIG_FUNC >= 0
 return t_tan<float, __TRIG_FUNC>(p_x);
#else
	return ::tanf(p_x);
#endif
}

inline double sinh(double p_x) {
	return ::sinh(p_x);
}
inline float sinh(float p_x) {
	return ::sinhf(p_x);
}

inline float sinc(float p_x) {
	return p_x == 0 ? 1 : ::sin(p_x) / p_x;
}
inline double sinc(double p_x) {
	return p_x == 0 ? 1 : ::sin(p_x) / p_x;
}

inline float sincn(float p_x) {
	return (float)sinc(Math_PI * p_x);
}
inline double sincn(double p_x) {
	return sinc(Math_PI * p_x);
}

inline double cosh(double p_x) {
	return ::cosh(p_x);
}
inline float cosh(float p_x) {
	return ::coshf(p_x);
}

inline double tanh(double p_x) {
	return ::tanh(p_x);
}
inline float tanh(float p_x) {
	return ::tanhf(p_x);
}

inline double asin(double p_x) {
#if __TRIG_FUNC >= 0
 return t_asin<double, __TRIG_FUNC>(p_x);
#else
	return ::asin(p_x);
#endif
}
inline float asin(float p_x) {
#if __TRIG_FUNC >= 0
 return t_asin<float, __TRIG_FUNC>(p_x);
#else
	return ::asinf(p_x);
#endif
}

inline double acos(double p_x) {
#if __TRIG_FUNC >= 0
 return t_acos<double, __TRIG_FUNC>(p_x);
#else
	return ::acos(p_x);
#endif
}
inline float acos(float p_x) {
#if __TRIG_FUNC >= 0
 return t_acos<float, __TRIG_FUNC>(p_x);
#else
	return ::acosf(p_x);
#endif
}

inline double atan(double p_x) {
#if __TRIG_FUNC >= 0
 return t_atan<double, __TRIG_FUNC>(p_x);
#else
	return ::atan(p_x);
#endif
}
inline float atan(float p_x) {
#if __TRIG_FUNC >= 0
 return t_atan<float, __TRIG_FUNC>(p_x);
#else
	return ::atanf(p_x);
#endif
}

inline double atan2(double p_y, double p_x) {
#if __TRIG_FUNC >= 0
 return t_atan2<double, __TRIG_FUNC>(p_y, p_x);
#else
	return ::atan2(p_y, p_x);
#endif
}
inline float atan2(float p_y, float p_x) {
#if __TRIG_FUNC >= 0
 return t_atan2<float, __TRIG_FUNC>(p_y, p_x);
#else 
	return ::atan2f(p_y, p_x);
#endif
}

inline double sqrt(double p_x) {
	return ::sqrt(p_x);
}
inline float sqrt(float p_x) {
	return ::sqrtf(p_x);
}

inline double pow(double p_x, double p_y) {
	return ::pow(p_x, p_y);
}
inline float pow(float p_x, float p_y) {
	return ::powf(p_x, p_y);
}

inline double log(double p_x) {
	return ::log(p_x);
}
inline float log(float p_x) {
	return ::logf(p_x);
}

inline float lerp(float minv, float maxv, float t) {
	return minv + t * (maxv - minv);
}
inline double lerp(double minv, double maxv, double t) {
	return minv + t * (maxv - minv);
}

inline double lerp_angle(double p_from, double p_to, double p_weight) {
	double difference = fmod(p_to - p_from, Math_TAU);
	double distance = fmod(2.0 * difference, Math_TAU) - difference;
	return p_from + distance * p_weight;
}
inline float lerp_angle(float p_from, float p_to, float p_weight) {
	float difference = fmod(p_to - p_from, (float)Math_TAU);
	float distance = fmod(2.0f * difference, (float)Math_TAU) - difference;
	return p_from + distance * p_weight;
}

inline double cubic_interpolate(double p_from, double p_to, double p_pre, double p_post, double p_weight) {
	return 0.5 *
			((p_from * 2.0) +
					(-p_pre + p_to) * p_weight +
					(2.0 * p_pre - 5.0 * p_from + 4.0 * p_to - p_post) * (p_weight * p_weight) +
					(-p_pre + 3.0 * p_from - 3.0 * p_to + p_post) * (p_weight * p_weight * p_weight));
}

inline float cubic_interpolate(float p_from, float p_to, float p_pre, float p_post, float p_weight) {
	return 0.5f *
			((p_from * 2.0f) +
					(-p_pre + p_to) * p_weight +
					(2.0f * p_pre - 5.0f * p_from + 4.0f * p_to - p_post) * (p_weight * p_weight) +
					(-p_pre + 3.0f * p_from - 3.0f * p_to + p_post) * (p_weight * p_weight * p_weight));
}

inline double cubic_interpolate_angle(double p_from, double p_to, double p_pre, double p_post, double p_weight) {
	double from_rot = fmod(p_from, Math_TAU);

	double pre_diff = fmod(p_pre - from_rot, Math_TAU);
	double pre_rot = from_rot + fmod(2.0 * pre_diff, Math_TAU) - pre_diff;

	double to_diff = fmod(p_to - from_rot, Math_TAU);
	double to_rot = from_rot + fmod(2.0 * to_diff, Math_TAU) - to_diff;

	double post_diff = fmod(p_post - to_rot, Math_TAU);
	double post_rot = to_rot + fmod(2.0 * post_diff, Math_TAU) - post_diff;

	return cubic_interpolate(from_rot, to_rot, pre_rot, post_rot, p_weight);
}

inline float cubic_interpolate_angle(float p_from, float p_to, float p_pre, float p_post, float p_weight) {
	float from_rot = fmod(p_from, (float)Math_TAU);

	float pre_diff = fmod(p_pre - from_rot, (float)Math_TAU);
	float pre_rot = from_rot + fmod(2.0f * pre_diff, (float)Math_TAU) - pre_diff;

	float to_diff = fmod(p_to - from_rot, (float)Math_TAU);
	float to_rot = from_rot + fmod(2.0f * to_diff, (float)Math_TAU) - to_diff;

	float post_diff = fmod(p_post - to_rot, (float)Math_TAU);
	float post_rot = to_rot + fmod(2.0f * post_diff, (float)Math_TAU) - post_diff;

	return cubic_interpolate(from_rot, to_rot, pre_rot, post_rot, p_weight);
}

inline double cubic_interpolate_in_time(double p_from, double p_to, double p_pre, double p_post, double p_weight,
		double p_to_t, double p_pre_t, double p_post_t) {
	/* Barry-Goldman method */
	double t = Math::lerp(0.0, p_to_t, p_weight);
	double a1 = Math::lerp(p_pre, p_from, p_pre_t == 0 ? 0.0 : (t - p_pre_t) / -p_pre_t);
	double a2 = Math::lerp(p_from, p_to, p_to_t == 0 ? 0.5 : t / p_to_t);
	double a3 = Math::lerp(p_to, p_post, p_post_t - p_to_t == 0 ? 1.0 : (t - p_to_t) / (p_post_t - p_to_t));
	double b1 = Math::lerp(a1, a2, p_to_t - p_pre_t == 0 ? 0.0 : (t - p_pre_t) / (p_to_t - p_pre_t));
	double b2 = Math::lerp(a2, a3, p_post_t == 0 ? 1.0 : t / p_post_t);
	return Math::lerp(b1, b2, p_to_t == 0 ? 0.5 : t / p_to_t);
}

inline float cubic_interpolate_in_time(float p_from, float p_to, float p_pre, float p_post, float p_weight,
		float p_to_t, float p_pre_t, float p_post_t) {
	/* Barry-Goldman method */
	float t = Math::lerp(0.0f, p_to_t, p_weight);
	float a1 = Math::lerp(p_pre, p_from, p_pre_t == 0 ? 0.0f : (t - p_pre_t) / -p_pre_t);
	float a2 = Math::lerp(p_from, p_to, p_to_t == 0 ? 0.5f : t / p_to_t);
	float a3 = Math::lerp(p_to, p_post, p_post_t - p_to_t == 0 ? 1.0f : (t - p_to_t) / (p_post_t - p_to_t));
	float b1 = Math::lerp(a1, a2, p_to_t - p_pre_t == 0 ? 0.0f : (t - p_pre_t) / (p_to_t - p_pre_t));
	float b2 = Math::lerp(a2, a3, p_post_t == 0 ? 1.0f : t / p_post_t);
	return Math::lerp(b1, b2, p_to_t == 0 ? 0.5f : t / p_to_t);
}

inline double cubic_interpolate_angle_in_time(double p_from, double p_to, double p_pre, double p_post, double p_weight,
		double p_to_t, double p_pre_t, double p_post_t) {
	double from_rot = fmod(p_from, Math_TAU);

	double pre_diff = fmod(p_pre - from_rot, Math_TAU);
	double pre_rot = from_rot + fmod(2.0 * pre_diff, Math_TAU) - pre_diff;

	double to_diff = fmod(p_to - from_rot, Math_TAU);
	double to_rot = from_rot + fmod(2.0 * to_diff, Math_TAU) - to_diff;

	double post_diff = fmod(p_post - to_rot, Math_TAU);
	double post_rot = to_rot + fmod(2.0 * post_diff, Math_TAU) - post_diff;

	return cubic_interpolate_in_time(from_rot, to_rot, pre_rot, post_rot, p_weight, p_to_t, p_pre_t, p_post_t);
}

inline float cubic_interpolate_angle_in_time(float p_from, float p_to, float p_pre, float p_post, float p_weight,
		float p_to_t, float p_pre_t, float p_post_t) {
	float from_rot = fmod(p_from, (float)Math_TAU);

	float pre_diff = fmod(p_pre - from_rot, (float)Math_TAU);
	float pre_rot = from_rot + fmod(2.0f * pre_diff, (float)Math_TAU) - pre_diff;

	float to_diff = fmod(p_to - from_rot, (float)Math_TAU);
	float to_rot = from_rot + fmod(2.0f * to_diff, (float)Math_TAU) - to_diff;

	float post_diff = fmod(p_post - to_rot, (float)Math_TAU);
	float post_rot = to_rot + fmod(2.0f * post_diff, (float)Math_TAU) - post_diff;

	return cubic_interpolate_in_time(from_rot, to_rot, pre_rot, post_rot, p_weight, p_to_t, p_pre_t, p_post_t);
}

inline double bezier_interpolate(double p_start, double p_control_1, double p_control_2, double p_end, double p_t) {
	/* Formula from Wikipedia article on Bezier curves. */
	double omt = (1.0 - p_t);
	double omt2 = omt * omt;
	double omt3 = omt2 * omt;
	double t2 = p_t * p_t;
	double t3 = t2 * p_t;

	return p_start * omt3 + p_control_1 * omt2 * p_t * 3.0 + p_control_2 * omt * t2 * 3.0 + p_end * t3;
}

inline float bezier_interpolate(float p_start, float p_control_1, float p_control_2, float p_end, float p_t) {
	/* Formula from Wikipedia article on Bezier curves. */
	float omt = (1.0f - p_t);
	float omt2 = omt * omt;
	float omt3 = omt2 * omt;
	float t2 = p_t * p_t;
	float t3 = t2 * p_t;

	return p_start * omt3 + p_control_1 * omt2 * p_t * 3.0f + p_control_2 * omt * t2 * 3.0f + p_end * t3;
}

template <typename T>
inline T clamp(T x, T minv, T maxv) {
	if (x < minv) {
		return minv;
	}
	if (x > maxv) {
		return maxv;
	}
	return x;
}

template <typename T>
inline T min(T a, T b) {
	return a < b ? a : b;
}

template <typename T>
inline T max(T a, T b) {
	return a > b ? a : b;
}

template <typename T>
inline T sign(T x) {
	return static_cast<T>(SIGN(x));
}

template <typename T>
inline T abs(T x) {
	return std::abs(x);
}

inline double deg_to_rad(double p_y) {
	return p_y * Math_PI / 180.0;
}
inline float deg_to_rad(float p_y) {
	return p_y * static_cast<float>(Math_PI) / 180.f;
}

inline double rad_to_deg(double p_y) {
	return p_y * 180.0 / Math_PI;
}
inline float rad_to_deg(float p_y) {
	return p_y * 180.f / static_cast<float>(Math_PI);
}

inline double inverse_lerp(double p_from, double p_to, double p_value) {
	return (p_value - p_from) / (p_to - p_from);
}
inline float inverse_lerp(float p_from, float p_to, float p_value) {
	return (p_value - p_from) / (p_to - p_from);
}

inline double remap(double p_value, double p_istart, double p_istop, double p_ostart, double p_ostop) {
	return Math::lerp(p_ostart, p_ostop, Math::inverse_lerp(p_istart, p_istop, p_value));
}
inline float remap(float p_value, float p_istart, float p_istop, float p_ostart, float p_ostop) {
	return Math::lerp(p_ostart, p_ostop, Math::inverse_lerp(p_istart, p_istop, p_value));
}

inline bool is_nan(float p_val) {
	return std::isnan(p_val);
}

inline bool is_nan(double p_val) {
	return std::isnan(p_val);
}

inline bool is_inf(float p_val) {
	return std::isinf(p_val);
}

inline bool is_inf(double p_val) {
	return std::isinf(p_val);
}

inline bool is_finite(float p_val) {
	return std::isfinite(p_val);
}

inline bool is_finite(double p_val) {
	return std::isfinite(p_val);
}

inline bool is_equal_approx(float a, float b) {
	// Check for exact equality first, required to handle "infinity" values.
	if (a == b) {
		return true;
	}
	// Then check for approximate equality.
	float tolerance = (float)CMP_EPSILON * abs(a);
	if (tolerance < (float)CMP_EPSILON) {
		tolerance = (float)CMP_EPSILON;
	}
	return abs(a - b) < tolerance;
}

inline bool is_equal_approx(float a, float b, float tolerance) {
	// Check for exact equality first, required to handle "infinity" values.
	if (a == b) {
		return true;
	}
	// Then check for approximate equality.
	return abs(a - b) < tolerance;
}

inline bool is_zero_approx(float s) {
	return abs(s) < (float)CMP_EPSILON;
}

inline bool is_equal_approx(double a, double b) {
	// Check for exact equality first, required to handle "infinity" values.
	if (a == b) {
		return true;
	}
	// Then check for approximate equality.
	double tolerance = CMP_EPSILON * abs(a);
	if (tolerance < CMP_EPSILON) {
		tolerance = CMP_EPSILON;
	}
	return abs(a - b) < tolerance;
}

inline bool is_equal_approx(double a, double b, double tolerance) {
	// Check for exact equality first, required to handle "infinity" values.
	if (a == b) {
		return true;
	}
	// Then check for approximate equality.
	return abs(a - b) < tolerance;
}

inline bool is_zero_approx(double s) {
	return abs(s) < CMP_EPSILON;
}

inline float absf(float g) {
	union {
		float f;
		uint32_t i;
	} u;

	u.f = g;
	u.i &= 2147483647u;
	return u.f;
}

inline double absd(double g) {
	union {
		double d;
		uint64_t i;
	} u;
	u.d = g;
	u.i &= (uint64_t)9223372036854775807ull;
	return u.d;
}

inline double smoothstep(double p_from, double p_to, double p_weight) {
	if (is_equal_approx(static_cast<real_t>(p_from), static_cast<real_t>(p_to))) {
		return p_from;
	}
	double x = clamp((p_weight - p_from) / (p_to - p_from), 0.0, 1.0);
	return x * x * (3.0 - 2.0 * x);
}
inline float smoothstep(float p_from, float p_to, float p_weight) {
	if (is_equal_approx(p_from, p_to)) {
		return p_from;
	}
	float x = clamp((p_weight - p_from) / (p_to - p_from), 0.0f, 1.0f);
	return x * x * (3.0f - 2.0f * x);
}

inline double move_toward(double p_from, double p_to, double p_delta) {
	return std::abs(p_to - p_from) <= p_delta ? p_to : p_from + sign(p_to - p_from) * p_delta;
}

inline float move_toward(float p_from, float p_to, float p_delta) {
	return std::abs(p_to - p_from) <= p_delta ? p_to : p_from + sign(p_to - p_from) * p_delta;
}

inline double linear2db(double p_linear) {
	return log(p_linear) * 8.6858896380650365530225783783321;
}
inline float linear2db(float p_linear) {
	return log(p_linear) * 8.6858896380650365530225783783321f;
}

inline double db2linear(double p_db) {
	return exp(p_db * 0.11512925464970228420089957273422);
}
inline float db2linear(float p_db) {
	return exp(p_db * 0.11512925464970228420089957273422f);
}

inline double round(double p_val) {
	return (p_val >= 0) ? floor(p_val + 0.5) : -floor(-p_val + 0.5);
}
inline float round(float p_val) {
	return (p_val >= 0) ? floor(p_val + 0.5f) : -floor(-p_val + 0.5f);
}

inline int64_t wrapi(int64_t value, int64_t min, int64_t max) {
	int64_t range = max - min;
	return range == 0 ? min : min + ((((value - min) % range) + range) % range);
}

inline float wrapf(real_t value, real_t min, real_t max) {
	const real_t range = max - min;
	return is_zero_approx(range) ? min : value - (range * floor((value - min) / range));
}

inline float fract(float value) {
	return value - floor(value);
}

inline double fract(double value) {
	return value - floor(value);
}

inline float pingpong(float value, float length) {
	return (length != 0.0f) ? abs(fract((value - length) / (length * 2.0f)) * length * 2.0f - length) : 0.0f;
}

inline double pingpong(double value, double length) {
	return (length != 0.0) ? abs(fract((value - length) / (length * 2.0)) * length * 2.0 - length) : 0.0;
}

// This function should be as fast as possible and rounding mode should not matter.
inline int fast_ftoi(float a) {
	static int b;

#if (defined(_WIN32_WINNT) && _WIN32_WINNT >= 0x0603) || WINAPI_FAMILY == WINAPI_FAMILY_PHONE_APP // windows 8 phone?
	b = (int)((a > 0.0) ? (a + 0.5) : (a - 0.5));

#elif defined(_MSC_VER) && _MSC_VER < 1800
	__asm fld a __asm fistp b
	/*#elif defined( __GNUC__ ) && ( defined( __i386__ ) || defined( __x86_64__ ) )
	// use AT&T inline assembly style, document that
	// we use memory as output (=m) and input (m)
	__asm__ __volatile__ (
	"flds %1        \n\t"
	"fistpl %0      \n\t"
	: "=m" (b)
	: "m" (a));*/

#else
	b = lrintf(a); // assuming everything but msvc 2012 or earlier has lrint
#endif
	return b;
}

inline double snapped(double p_value, double p_step) {
	if (p_step != 0) {
		p_value = Math::floor(p_value / p_step + 0.5) * p_step;
	}
	return p_value;
}

inline float snap_scalar(float p_offset, float p_step, float p_target) {
	return p_step != 0 ? Math::snapped(p_target - p_offset, p_step) + p_offset : p_target;
}

inline float snap_scalar_separation(float p_offset, float p_step, float p_target, float p_separation) {
	if (p_step != 0) {
		float a = Math::snapped(p_target - p_offset, p_step + p_separation) + p_offset;
		float b = a;
		if (p_target >= 0) {
			b -= p_separation;
		} else {
			b += p_step;
		}
		return (Math::abs(p_target - a) < Math::abs(p_target - b)) ? a : b;
	}
	return p_target;
}

} // namespace Math
} // namespace godot

#endif // GODOT_MATH_HPP
