/** @file
	@ingroup jamoma2
 
	@brief Unit test for the Delay class
 
	@author Timothy Place
 
	@copyright Copyright © 2000-2015 by Jamoma authors and contributors @n
	This code is licensed under the terms of the "BSD 3-Clause License" @n
	https://github.com/jamoma/jamoma2/blob/master/LICENSE.md @n
 */

#include "Jamoma.h"


class DelayTest {
	
	Jamoma::UnitTest<DelayTest>*	mTest;
	
public:
	DelayTest(Jamoma::UnitTest<DelayTest>* test)
	: mTest(test)
	{
		testDelayGreaterThanOneVectorSize();
	}

	
	void testDelayGreaterThanOneVectorSize()
	{
		Jamoma::SampleBundle zero(2, 64);
		
		Jamoma::UnitImpulse impulse;
		impulse.channelCount = 2;
		impulse.frameCount = 64;

		Jamoma::Delay my_delay;
		my_delay.sampleRate = 44100;
		my_delay.size = 100;
		
		Jamoma::SampleBundle out_samples1 = my_delay( impulse() );
		Jamoma::SampleBundle out_samples2 = my_delay( zero );
		Jamoma::SampleBundle out_samples3 = my_delay( zero );
		
		int badSampleCount = 0;
		
		// first 64 samples should all be zero
		for (auto& channel : out_samples1) {
			for (auto& sample : channel) {
				if (sample != 0.0) {
					badSampleCount++;
					std::cout << "bad sample " << " is " << sample << std::endl;
				}
			}
		}
		
		// 100 samples later should be 1.0
		// note: this is not the 100th sample, but 100 samples after sample 0 -- a delay of 0 is sample 0
		for (auto& channel : out_samples2) {
			int i = 0;
			for (auto& sample : channel) {
				if (i == 100-64) {
					if (sample != 1.0) {
						badSampleCount++;
						std::cout << "sample " << i << " is " << sample << std::endl;
					}
				}
				else if (sample != 0.0) {
					badSampleCount++;
					std::cout << "sample " << i << " is " << sample << std::endl;
				}
				++i;
			}
		}
		
		// last 64 samples should all be zero
		for (auto& channel : out_samples3) {
			for (auto& sample : channel) {
				if (sample != 0.0) {
					badSampleCount++;
					std::cout << "bad sample " << " is " << sample << std::endl;
				}
			}
		}

		std::cout << "the output has " << badSampleCount << " bad samples" << std::endl;
		mTest->TEST_ASSERT("Bad Sample Count", badSampleCount == 0);
	}
	
};


int main(int argc, const char * argv[])
{
	Jamoma::UnitTest<DelayTest>	aUnitTestInstance;
	return aUnitTestInstance.failureCount();
}