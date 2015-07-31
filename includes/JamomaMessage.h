/** @file
	
	@ingroup jamoma2
	
	@brief Define a message to a JamomaObject
 
	@author Timothy Place
	
	@copyright Copyright © 2015 by Jamoma authors and contributors @n
	This code is licensed under the terms of the "BSD 3-Clause License" @n
	https://github.com/jamoma/jamoma2/blob/master/LICENSE.md @n
 */

#pragma once

#include "JamomaObject.h"

namespace Jamoma {
	
	
	/** Define a message that can be sent to a Jamoma::Object. 
		Any such message may also be called directly as method of the Object in C++.
	 */
	class Message {
	public:
		
		// TODO: Do we really need this variation of the constructor?  YAGNI!
		Message(String name, Function fn)
		: mFunction(fn)
		{}
		
		
		/** Define a Message.
			@param	name		The name of the Message when called on your Object dynamically.
			@param	synopsis	A description of what the Message does.
			@param	fn			The actual function to perform. May be a function pointer, functor, or lamba.
		 */
		Message(String name, Synopsis synopsis, Function fn)
		: mFunction(fn)
		, mSynopsis(synopsis)
		{}
		
		
		/** Execute the Message directly as a method of the C++ class. */
		Error operator ()() {
			mFunction();
			return Error::none;	// TODO: what's the point of returning an error code in this case?
		}
		
	private:
		Function		mFunction;
		Synopsis		mSynopsis;
	};
	

} // namespace Jamoma
