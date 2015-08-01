/** @file
	
	@ingroup jamoma2
	
	@brief Apply a first-order low-pass filter to samples.
	
	@author Nathan Wolek
	
	@copyright Copyright © 2015 by Jamoma authors and contributors @n
	This code is licensed under the terms of the "BSD 3-Clause License" @n
	https://github.com/jamoma/jamoma2/blob/master/LICENSE.md @n
 */

#pragma once

#include "JamomaAudioObject.h"

namespace Jamoma {
    
    /**	This AudioObject applies a basic <a href="https://en.wikipedia.org/wiki/Low-pass_filter">low-pass filter</a> to a Sample or SampleBundle.
     */
    class LowpassOnePole : public AudioObject {
        
        AdaptingSampleBundle		mY1 = {this};					///< previous output sample for each channel
        double		mCoefficientF;			///< filter coefficient
        double		mOneMinusCoefficientF;	///< 1 - mCoefficientF
        
    public:
        static constexpr Classname classname = { "lowpass.1" };
        static constexpr auto tags = { "dspEffectsLib", "audio", "processor", "filter", "lowpass" };
        
        /** Filter coefficient.
         */
        Parameter<double, RangeLimit::clip> coefficient = { this,
                                                            "coefficient",
                                                            0.5,
                                                            Range<double>(0.0, 1.0),
                                                            [this]{
                                                                mCoefficientF = coefficient;
                                                                mOneMinusCoefficientF = 1 - coefficient;
                                                            }
        };
        
        Sample operator()(Sample x, int channel)
        {
            
        };
        
    };
    
    
}  // namespace Jamoma