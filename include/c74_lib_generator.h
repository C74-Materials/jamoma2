/// @file
///	@ingroup 	minlib
/// @author		Timothy Place, Nathan Wolek
///	@copyright	Copyright (c) 2017, Cycling '74
///	@license	Usage of this file and its contents is governed by the MIT License

#pragma once


namespace c74 {
namespace min {
    /** Defines several functions for use with <a href="http://en.cppreference.com/w/cpp/algorithm/generate">std::generate</a> to fill vectors with common shapes used in computer sound.
     */
	namespace Generator {

        /** Generates a line from -1 to 1 with consistent slope
         @param T       render output as this datatype. algorithm was designed to assume the use of floating point.
         @param size    size of the target vector
         */
        template <typename T>
        class Ramp {
        public:
            Ramp (int size)
            : mCycleSize(size)
            {
                //TODO: we need way to protect against zero. static_assert did not work.
            }

            T operator()() {
                ++mCurrent;
                return ( T(mCurrent) * 2.0 / mCycleSize) - 1.0;
            }

        private:
            int mCurrent = -1;
            int mCycleSize; // required by constructor
        };

        /** Generates an ideal sawtooth waveform from -1 to 1. Not anti-aliased.
         @param T       render output as this datatype. algorithm was designed to assume the use of floating point.
         @param size    size of the target vector
         */
        template <typename T>
        using Sawtooth     = Generator::Ramp<T>;

        /** Generates a line from 0 to 1 with consistent slope
         @param T       render output as this datatype. algorithm was designed to assume the use of floating point.
         @param size    size of the target vector
         */
		template <typename T>
		class UnipolarRamp {
		public:
			UnipolarRamp (int size)
			: mCycleSize(size)
			{
                //TODO: we need way to protect against zero. static_assert did not work.
            }

			T operator()() {
				++mCurrent;
				return T(mCurrent) / mCycleSize;
			}

		private:
			int mCurrent = -1;
			int mCycleSize; // required by constructor
		};

        /** Generates an ideal sawtooth waveform from 0 to 1. Not anti-aliased.
         @param T       render output as this datatype. algorithm was designed to assume the use of floating point.
         @param size    size of the target vector
         */
        template <typename T>
        using UnipolarSawtooth     = Generator::UnipolarRamp<T>;

        /** Generates a sine wave constrained between -1 to 1
         @param T       render output as this datatype. algorithm was designed to assume the use of floating point.
         @param size    size of the target vector
         */
		template <typename T>
		class Sine {
		public:
			Sine (int size)
			: mCycleSize(size)
            {
                //TODO: we need way to protect against zero. static_assert did not work.
            }

			T operator()() {
				++mCurrent;
				auto output = std::sin(mCurrent * kTwoPi / mCycleSize);
				return T(output);
			}

		private:
			int mCurrent = -1;
			int mCycleSize; // required by constructor
		};

        /** Generates a sine wave constrained between 0 to 1
         @param T       render output as this datatype. algorithm was designed to assume the use of floating point.
         @param size    size of the target vector
         */
        template <typename T>
        class UnipolarSine {
        public:
            UnipolarSine (int size)
            : mCycleSize(size)
            {
                //TODO: we need way to protect against zero. static_assert did not work.
            }

            T operator()() {
                ++mCurrent;
                auto output = 0.5 * std::sin(mCurrent * kTwoPi / mCycleSize);
                return T(output) + 0.5;
            }

        private:
            int mCurrent = -1;
            int mCycleSize; // required by constructor
        };

        /** Generates a cosine wave constrained between -1 to 1
         @param T       render output as this datatype. algorithm was designed to assume the use of floating point.
         @param size    size of the target vector
         */
        template <typename T>
        class Cosine {
        public:
            Cosine (int size)
            : mCycleSize(size)
            {
                //TODO: we need way to protect against zero. static_assert did not work.
            }

            T operator()() {
                ++mCurrent;
                auto output = std::cos(mCurrent * kTwoPi / mCycleSize);
                return T(output);
            }

        private:
            int mCurrent = -1;
            int mCycleSize; // required by constructor
        };

        /** Generates a cosine wave constrained between 0 to 1
         @param T       render output as this datatype. algorithm was designed to assume the use of floating point.
         @param size    size of the target vector
         */
        template <typename T>
        class UnipolarCosine {
        public:
            UnipolarCosine (int size)
            : mCycleSize(size)
            {
                //TODO: we need way to protect against zero. static_assert did not work.
            }

            T operator()() {
                ++mCurrent;
                auto output = 0.5 + 0.5 * std::cos(mCurrent * kTwoPi / mCycleSize);
                return T(output);
            }

        private:
            int mCurrent = -1;
            int mCycleSize; // required by constructor
        };


        /** Generates a triangle wave constrained between -1 to 1
         @param T       render output as this datatype. algorithm was designed to assume the use of floating point.
         @param size    size of the target vector
         */
		template <typename T>
		class Triangle {
		public:
			Triangle (int size)
			: mCycleSize(size)
            {
                //TODO: we need way to protect against zero. static_assert did not work.
            }

			T operator()() {
				T out = 0.0;
				++mCurrent;

				if (mCurrent <= mCycleSize/4)
					out = 4.0 * mCurrent / mCycleSize;
				else if (mCurrent >= 3 * mCycleSize / 4)
					out = -4.0 + 4.0 * mCurrent / mCycleSize;
				else
					out = 2.0 - 4.0 * mCurrent / mCycleSize;
				return out;
			}

		private:
			int mCurrent = -1;
			int mCycleSize; // required by constructor
		};

        /** Generates a triangle wave constrained between 0 to 1
         @param T       render output as this datatype. algorithm was designed to assume the use of floating point.
         @param size    size of the target vector
         */
        template <typename T>
        class UnipolarTriangle {
        public:
            UnipolarTriangle (int size)
            : mCycleSize(size)
            {
                //TODO: we need way to protect against zero. static_assert did not work.
            }

            T operator()() {
                T out = 0.0;
                ++mCurrent;

                if (mCurrent <= mCycleSize/4)
                    out = 2.0 * mCurrent / mCycleSize;
                else if (mCurrent >= 3 * mCycleSize / 4)
                    out = -2.0 + 2.0 * mCurrent / mCycleSize;
                else
                    out = 1.0 - 2.0 * mCurrent / mCycleSize;
                return 0.5 + out;
            }

        private:
            int mCurrent = -1;
            int mCycleSize; // required by constructor
        };

	} // namespace Generator
}}  // namespace c74::min
