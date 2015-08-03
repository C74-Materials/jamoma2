/** @file
	
	@ingroup jamoma2
	
	@brief Apply a fourth-order low-pass filter to samples.
	
	@author Timothy Place, Nathan Wolek
	
	@copyright Copyright © 2015 by Jamoma authors and contributors @n
	This code is licensed under the terms of the "BSD 3-Clause License" @n
	https://github.com/jamoma/jamoma2/blob/master/LICENSE.md @n
 */

#pragma once

#include "JamomaAudioObject.h"

namespace Jamoma {

	
	/**	This AudioObject applies a resonant <a href="https://en.wikipedia.org/wiki/Low-pass_filter">low-pass filter</a> to a Sample or SampleBundle.
		A fourth-order algorithm based on the Moog VCF is used to acheive a relatively flat passband.
        Subtle modifications from <a href="http://www.musicdsp.org/archive.php?classid=3#26">the original source</a> to improve control of the frequency cutoff.

		@warning The cutoff frequency is only stable below about one-fourth the sampling rate.
        In addition, it will tend to sound flat below 1000 Hz and go sharp above 1000 Hz.
	 */
	class LowpassFourPole : public AudioObject {
		
		// TODO: change the below to use 2 bundles of 4 frames each
		
		AdaptingSampleBundle		mX1 = {this};					///< previous input sample for each channel
		AdaptingSampleBundle		mX2 = {this};					///< 2nd previous input sample for each channel
		AdaptingSampleBundle		mX3 = {this};					///< 3rd previous input sample for each channel
		AdaptingSampleBundle		mX4 = {this};					///< 4th previous input sample for each channel
		AdaptingSampleBundle		mY1 = {this};					///< previous output sample for each channel
		AdaptingSampleBundle		mY2 = {this};					///< 2nd previous output sample for each channel
		AdaptingSampleBundle		mY3 = {this};					///< 3rd previous output sample for each channel
		AdaptingSampleBundle		mY4 = {this};					///< 4th previous output sample for each channel

		double		mDeciResonance;			///< attrResonance * 0.1
		double		mCoefficientF;			///< filter coefficient
		double		mCoefficientSquaredF;   ///< mCoefficientF * mCoefficientF
		double		mOneMinusCoefficientF;	///< 1 - mCoefficientF
		double		mCoefficientFB;			///< filter coefficient
		double		mCoefficientG;			///< filter coefficient
		
	public:
		static constexpr Classname classname = { "lowpass.4" };
		static constexpr auto tags = { "dspEffectsLib", "audio", "processor", "filter", "lowpass" };
		
		
		/** Filter cutoff frequency.
            Controls the boundary between low frequency pass band and high frequency stop band.
            The frequency setting will tend to sound flat below 1000 Hz and goes sharp above 1000 Hz.
            In addition, the cutoff frequency is only stable below about one-fourth the sampling rate.
         */
		Parameter<double, RangeLimit::clip>	frequency = {	this,
															"frequency",
															1000.0,
															Range<double>(20.0, sampleRate * 0.475),
															[this]{
																//double	radians = hertzToRadians(frequency);
																double  fNormalizedToHalfNyquist = frequency * 4 / sampleRate;
                                                                mCoefficientF = fNormalizedToHalfNyquist * 1.4716;
																mCoefficientSquaredF = mCoefficientF * mCoefficientF;
																mOneMinusCoefficientF = 1.0 - mCoefficientF;
																calculateCoefficients();
															}
		};
		// addAttributeProperty(Frequency,			description,	TT("Cutoff Frequency in Hertz"));
		
		
		/** Filter resonance.
            Controls the prominence of the cutoff frequency.
            The usable range for this parameter is between 0.0 and 40.0.
         */
		Parameter<double>	q = { this, "q", 1.0,
			[this]{
				mDeciResonance = q * 0.1;
				calculateCoefficients();
			}
		};
		// addAttributeProperty(Q,			range,			TTValue(0.01, 100.0));
		// addAttributeProperty(Q,			rangeChecking,	TT("cliplow"));
		// addAttributeProperty(Q,			description,	TT("Strength of Resonance Near the Cutoff Frequency"));
		
		
		/**	This algorithm uses an IIR filter, meaning that it relies on feedback.  If the filter should
		 *	not be producing any signal (such as turning audio off and then back on in a host) or if the
		 *	feedback has become corrupted (such as might happen if a NaN is fed in) then it may be
		 *	neccesary to clear the filter by calling this method.
		 *	@return Returns a TTErr error code.												*/

		Message			clear = { "clear",
									Synopsis("Reset the Filter History"),
									[this]{
										mX1.fill(0.0);
										mX2.fill(0.0);
										mX3.fill(0.0);
										mX4.fill(0.0);
										mY1.fill(0.0);
										mY2.fill(0.0);
										mY3.fill(0.0);
										mY4.fill(0.0);
									}
		};
		
		
		void calculateCoefficients()
		{
			mCoefficientFB = mDeciResonance * (1.0 - 0.4 * mCoefficientSquaredF);
			mCoefficientG = /*TTAntiDenormal*/(0.35013 * (mCoefficientSquaredF * mCoefficientSquaredF));
		}

		
		Sample operator()(Sample x)
		{
			return (*this)(x, 0);
		}
		
		
		Sample operator()(Sample x, int channel)
		{
			// the original source of this algorithm appears to be here:
            // http://musicdsp.org/archive.php?classid=3#26
            
            Sample y = x;
			
			y -= mY4[channel][0] * mCoefficientFB;
			//TTZeroDenormal(y);
			y *= mCoefficientG;
			
			mY1[channel][0] = y + 0.3 * mX1[channel][0] + mOneMinusCoefficientF * mY1[channel][0];
			//TTZeroDenormal(mY1[channel]);
			mX1[channel][0] = y;
			mY2[channel][0] = mY1[channel][0] + 0.3 * mX2[channel][0] + mOneMinusCoefficientF * mY2[channel][0];
			//TTZeroDenormal(mY2[channel]);
			mX2[channel][0] = mY1[channel][0];
			mY3[channel][0] = mY2[channel][0] + 0.3 * mX3[channel][0] + mOneMinusCoefficientF * mY3[channel][0];
			//TTZeroDenormal(mY3[channel] );
			mX3[channel][0] = mY2[channel][0];
			mY4[channel][0] = mY3[channel][0] + 0.3 * mX4[channel][0] + mOneMinusCoefficientF * mY4[channel][0];
			//TTZeroDenormal(mY4[channel]);
			mX4[channel][0] = mY3[channel][0];
			y = mY4[channel][0];
			
			return y;
		}
		
		
		SharedSampleBundleGroup operator()(const SampleBundle& x)
		{
			auto out = mOutput;
			out.adapt(x);
			channelCount = (int)x.channelCount();
			
			for (int channel=0; channel < x.channelCount(); channel++) {
				for	(int i=0; i < x.frameCount(); i++)
					out[0][channel][i] = (*this)(x[channel][i], channel);
			}
			return out;
		}
		
	};
	

} // namespace Jamoma
