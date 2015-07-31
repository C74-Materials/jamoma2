/** @file
	
	@ingroup jamoma2
	
	@brief Define a parameter of a JamomaObject
 
	@author Timothy Place
	
	@copyright Copyright © 2015 by Jamoma authors and contributors @n
	This code is licensed under the terms of the "BSD 3-Clause License" @n
	https://github.com/jamoma/jamoma2/blob/master/LICENSE.md @n
 */

#pragma once

#include "JamomaObject.h"

namespace Jamoma {
	
	
	/** ParameterBase allows us a unified means for referencing a Parameter of any template-specialization.
		All Parameters must be defined using a Parameter<type> specialization and this base class is for internal use only.
	 */
	class ParameterBase {
	protected:
		Object*				mOwner;			///< The owning Jamoma::Object to which this Parameter belongs.
		String				mName;			///< The name of this Parameter as it would be addressed dynamically.
		Synopsis			mSynopsis;		///< A description of what this Parameter represents.
		RangeLimit			mRangeLimit;	///< The behavior applied to values sent to this parameter if they are outside of the suggested Range.
		Function			mSetter;		///< A function to be executed after the parameter's value has been set.
		// getter
		

		ParameterBase(Object* owner, const String& name, const Synopsis& synopsis, const RangeLimit& rangeLimit, const Function& setter)
		: mOwner(owner)
		, mName(name)
		, mSynopsis(synopsis)
		, mRangeLimit(rangeLimit)
		, mSetter(setter)
		{}
		
		virtual ~ParameterBase()
		{}
		
	public:
		String& name()
		{
			return mName;
		}
		
		friend bool operator == (const ParameterBase& lhs, const char* rhs)
		{
			return lhs.mName == rhs;
		}
		
		virtual ParameterBase& operator = (const ValueBase& input) = 0;
	};
	
	
	// NOTE: It would be nice to have an intermediate class here that dealt with the type specialization
	// from which RangeLimit specialization could inherit.
	// It doesn't work, however, because the operator overloads won't be inherited and thus we have to duplicate all the code anyway.
	
	
#pragma mark -
#pragma mark Parameters that Don't Limit
	
	
	/** Defines a Parameter with no special behavior applied if the supplied values are out of range.
	 */
	template <class T, RangeLimit rangeLimit = RangeLimit::none>
	class Parameter : public ParameterBase {
		T					mValue;
		Range<T>			mRange;
		
		void set(T input)
		{
			mValue = input;
			if (mSetter)
				mSetter();
		}
		
	public:
		Parameter() = delete;		// Can't create an unitialized Parameter
		
		
		Parameter(Object* owner, String name, T initial, Function setter = nullptr)
		: ParameterBase(owner, name, "", rangeLimit, setter)
		{
			// 1. iterate args
			// 2. determine their types
			// 3. do something appropriate for their given type
			// 4. can we make this whole process constexpr ?
			
			owner->parameters[name] = this;
			set(initial);
		}
		
		
		// setter
		Parameter& operator = (T input)
		{
			set(input);
			return *this;
		}
		
		
		// setter for case when input is a generic value
		Parameter& operator = (const ValueBase& input)
		{
			set(input);
			return *this;
		}
		
		
		// assign *values* from one attribute to another
		Parameter& operator = (Parameter& input)
		{
			set(input.mValue);
			return *this;
		}
		
		
		// getter
		operator T() const
		{
			return mValue;
		}
	};
	
	
#pragma mark -
#pragma mark Parameters that Clip
	

	/** Defines a Parameter where values are limited (clipped) to the min and max of the suggested Range.
	 */
	template<class T>
	class Parameter<T, RangeLimit::clip> : public ParameterBase {
		T					mValue;
		Range<T>			mRange;
	
		void set(T input)
		{
			mValue = Limit(input, mRange.first, mRange.second);
			if (mSetter)
				mSetter();
		}
		
	public:
		Parameter() = delete;		// Can't create an unitialized Parameter

		
		Parameter(Object* owner, String name, T initial, Range<T> range, Function setter = nullptr)
		: ParameterBase(owner, name, "", RangeLimit::clip, setter)
		, mRange(range)
		{
			// see comments above regarding arg parsing
			
			owner->parameters[name] = this;
			set(initial);
		}

		
		// setter
		Parameter& operator = (T input)
		{
			set(input);
			return *this;
		}
		
		
		// setter for case when input is a generic value
		Parameter& operator = (const ValueBase& input)
		{
			set(input);
			return *this;
		}
		
		
		// assign *values* from one attribute to another
		Parameter& operator = (Parameter& input)
		{
			set(input.mValue);
			return *this;
		}
		
		
		// getter
		operator T() const
		{
			return mValue;
		}
	};


} // namespace Jamoma
