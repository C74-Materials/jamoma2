/// @file
///	@ingroup 	minlib
/// @author		Timothy Place, Nils Peters, Tristan Matthews
///	@copyright	Copyright (c) 2017, Cycling '74
///	@license	Usage of this file and its contents is governed by the MIT License

#pragma once


namespace c74 {
namespace min {
namespace lib {


	///	Determine if a value is a power-of-two. Only works for ints.
	// TODO: static_assert is_integral

	template<class T>
	bool is_power_of_two(T value) {
		return (value > 0) && ((value & (value-1)) == 0);
	}


	///	Limit input to power-of-two values.
	/// Non-power-of-two values are increased to the next-highest power-of-two upon return.
	/// Only works for ints up to 32-bits.
	/// @seealso http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
	// TODO: static_assert correct type

	template<class T>
	void limit_to_power_of_two(T value) {
		value--;
		value |= value >> 1;
		value |= value >> 2;
		value |= value >> 4;
		value |= value >> 8;
		value |= value >> 16;
		++value;
		return value;
	}


	/// This routine wrapping around the range as much as necessary

	template<class T>
	T wrap(T input, const T low_bound, const T high_bound) {
		if ((input >= low_bound) && (input < high_bound))
			return input; //nothing to wrap
		else if (input - low_bound >= 0)
			return (fmod((double)input  - low_bound, fabs((double)low_bound - high_bound)) + low_bound);
		else
			return (-1.0 * fmod(-1.0 * (input - low_bound), fabs((double)low_bound - high_bound)) + high_bound);
	}


	/// A fast routine for wrapping around the range once.
	/// This is faster than doing an expensive module, where you know the range of the input
	/// will not equal or exceed twice the range.

	template<class T>
	T wrap_once(T input, const T low_bound, const T high_bound) {
		if ((input >= low_bound) && (input < high_bound))
			return input;
		else if (input >= high_bound)
			return ((low_bound - 1) + (input - high_bound));
		else
			return ((high_bound + 1) - (low_bound - input));
	}


	/// This routine folds numbers into the data range

	template<typename T>
	T fold(T input, const T low_bound, const T high_bound) {
		double foldRange;

		if ((input >= low_bound) && (input <= high_bound))
			return input; //nothing to fold
		else {
			foldRange = 2 * fabs((double)low_bound - high_bound);
			return fabs(remainder(input - low_bound, foldRange)) + low_bound;
		}
	}


	/// A utility for scaling one range of values onto another range of values.
	template<class T>
	static T scale(T value, T inlow, T inhigh, T outlow, T outhigh) {
		double inscale, outdiff;

		inscale = 1 / (inhigh - inlow);
		outdiff = outhigh - outlow;

		value = (value - inlow) * inscale;
		value = (value * outdiff) + outlow;
		return(value);
	}


	/// Defines several functions for constraining values within specified boundaries and preventing unwanted values.
	/// A variety of behaviors are offered, including clipping, wrapping and folding.
	/// Exercise caution when using the functions defined here with unsigned values.
	/// Negative, signed integers have the potential to become very large numbers when casting to unsigned integers.
	/// This can cause errors during a boundary check, such as values clipping to the high boundary instead of the
	/// low boundary or numerous iterations of loop to bring a wrapped value back into the acceptable range.

	namespace limit {

		template <typename T>
		class base {
		public:
			base() = delete;
		};


		template <typename T>
		class none : public base<T> {
		public:
			static T apply(T input, T low, T high) {
				return input;
			}

			T operator()(T input, T low, T high) {
				return apply(input, low, high);
			}
		};


		template <typename T>
		class clip : public base<T> {
		public:
			static T apply(T input, T low, T high) {
				return max::clamp<T>(input, low, high);
			}

			T operator()(T input, T low, T high) {
				return apply(input, low, high);
			}
		};


		template <typename T>
		class wrap : public base<T> {
		public:
			static T apply(T input, T low, T high) {
				return lib::wrap(input, low, high);
			}

			T operator()(T input, T low, T high) {
				return apply(input, low, high);
			}
		};

		template <typename T>
		class wrap_once : public base<T> {
		public:
			static T apply(T input, T low, T high) {
				return lib::wrap_once(input, low, high);
			}

			T operator()(T input, T low, T high) {
				return apply(input, low, high);
			}
		};

		template <typename T>
		class fold : public base<T> {
		public:
			static T apply(T input, T low, T high) {
				return lib::fold(input, low, high);
			}

			T operator()(T input, T low, T high) {
				return apply(input, low, high);
			}
		};

	} // namespace Limit


}}}  // namespace c74::min::lib
