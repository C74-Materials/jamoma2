/** @file
	@ingroup 	jamoma2
 
	@brief 		Unit test for the Dcblocker class
 
	@author		Timothy Place, Nathan Wolek
	@copyright	Copyright (c) 2005-2015 The Jamoma Group, http://jamoma.org.
	@license	This project is released under the terms of the MIT License.

 */

#include "Jamoma.h"

namespace Jamoma {

	class SampleBundleTest {
		
		UnitTest<SampleBundleTest>*	mTest;
		
	public:
		SampleBundleTest(Jamoma::UnitTest<SampleBundleTest>* test)
		: mTest(test)
		{
			testBasic();
		}

		
		void testBasic()
		{
			Jamoma::Gain			my_gain;
			Jamoma::SampleBundle	in_samples(2, 8);
			
			my_gain.gain = 0.5;
			
			in_samples.fill(1.0);
			auto out_samples = my_gain(in_samples);
			
			my_gain.gain = 0.25;
			in_samples = out_samples;
			out_samples = my_gain(in_samples);
						
			auto bar = in_samples[0][0];
			mTest->TEST_ASSERT("in_sample casting operator", mTest->compare(bar, 0.5));
			
			in_samples[0][0] = 2.0;
			auto foo = in_samples[0][0];
			mTest->TEST_ASSERT("setting and getting", mTest->compare(foo, 2.0));
		}
	};

} // namespace Jamoma


int main(int argc, const char * argv[])
{
	Jamoma::UnitTest<Jamoma::SampleBundleTest>	aUnitTestInstance;
	return aUnitTestInstance.failureCount();
}
