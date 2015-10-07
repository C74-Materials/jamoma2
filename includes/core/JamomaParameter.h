/** @file
	
	@ingroup 	jamoma2
	
	@brief 		Define a parameter of a JamomaObject
 
	@author		Timothy Place
	@copyright	Copyright (c) 2005-2015 The Jamoma Group, http://jamoma.org.
	@license	This project is released under the terms of the MIT License.
 */

#pragma once

#include "JamomaObject.h"

namespace Jamoma {
	
	
	/** ParameterBase allows us a unified means for referencing a Parameter of any template-specialization.
		All Parameters must be defined using a Parameter<type> specialization and this base class is for internal use only.
	 */
	class ParameterBase {
	protected:
		Object*					mOwner;			///< The owning Jamoma::Object to which this Parameter belongs.
		String					mName;			///< The name of this Parameter as it would be addressed dynamically.
		Synopsis				mSynopsis;		///< A description of what this Parameter represents.
		RangeLimit				mRangeLimit;	///< The behavior applied to values sent to this parameter if they are outside of the suggested Range.
		Function				mSetter;		///< A function to be executed after the parameter's value has been set.
		std::vector<Observer*>	mObservers;		///< Objects receiving notifications when this parameter has been set.

		// TODO: the above raw pointer will lead to dangling references ?!?!?!?!

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
		
		virtual ParameterBase& operator = (const VarBase& input) = 0;
		
		
		void addObserver(Observer& anObserver)
		{
			mObservers.push_back(&anObserver);
		}
		
		
		void removeObserver(Observer& anObserver)
		{
			// documentation of the below: https://en.wikipedia.org/wiki/Erase%E2%80%93remove_idiom
			mObservers.erase(std::remove(mObservers.begin(), mObservers.end(), &anObserver), mObservers.end());
		}
		
	};
	
	
	// NOTE: It would be nice to have an intermediate class here that dealt with the type specialization
	// from which RangeLimit specialization could inherit.
	// It doesn't work, however, because the operator overloads won't be inherited and thus we have to duplicate all the code anyway.
	
	
#pragma mark -
#pragma mark Parameters that Fold or Wrap
	
	
	/** Defines a Parameter where values are limited (clipped) to the min and max of the suggested Range.
	 */
	template <class T, class Dataspace_Class = Dataspace::None<T, Dataspace::NoneUnit::nothing>, RangeLimit rangeLimit = RangeLimit::none>
	class Parameter : public ParameterBase {
		
		T					mValue;
		Range<T>			mRange;
		Dataspace_Class		mDataspace;
		
		
		// naked values have no dataspace unit specified, so set them directly
		void set(T input)
		{
			switch (mRangeLimit) {
				case RangeLimit::wrap:
					mValue = Wrap(input, mRange.first, mRange.second);
					break;
				case RangeLimit::fold:
					mValue = Fold(input, mRange.first, mRange.second);
					break;
				case RangeLimit::clip:
					mValue = Limit(input, mRange.first, mRange.second);
					break;
				default:
					mValue = input;
					break;
			}
			if (mSetter)
				mSetter();
			for (auto& observer : mObservers)
				(*observer)();
		}
		
		
		// set values using a dataspace conversion
		void set(T input, Unit unit)
		{
			set(mDataspace(input, (uint32_t)unit));
		}
		
		
	public:
		Parameter() = delete;		// Can't create an unitialized Parameter
		
		
		Parameter(Object* owner, String name, T initial, Range<T> range, Function setter = nullptr)
		: ParameterBase(owner, name, "", rangeLimit, setter)
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
		
		
		// setter w/ unit
		Parameter& operator = (const std::pair<T, Unit> input)
		{
			set(input.first, input.second);
			return *this;
		}
		
		
		// setter for case when input is a generic value
		// TODO: if a value has 2 members then do we use the last one as a unit? perhaps it needs some metadata so that we know?
		Parameter& operator = (const VarBase& input)
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
#pragma mark Parameters that Don't Limit
	
	
	/** Defines a Parameter with no special behavior applied if the supplied values are out of range.
	 */
	template <class T, class Dataspace_Class>
	class Parameter<T, Dataspace_Class, RangeLimit::none> : public ParameterBase {
		T					mValue;
		Range<T>			mRange;
		Dataspace_Class		mDataspace;
		
		
		// naked values have no dataspace unit specified, so set them directly
		void set(T input)
		{
			mValue = input;
			if (mSetter)
				mSetter();
			for (auto& observer : mObservers)
				(*observer)();
		}
		

		// set values using a dataspace conversion
		void set(const T& input, const Unit& unit)
		{
			set(mDataspace(input, (uint32_t&)unit));
		}
		
		
	public:
		Parameter() = delete;		// Can't create an unitialized Parameter
		
		
		Parameter(Object* owner, String name, T initial,  Function setter = nullptr)
		: ParameterBase(owner, name, "", RangeLimit::none, setter)
		{
			// 1. iterate args
			// 2. determine their types
			// 3. do something appropriate for their given type
			// 4. can we make this whole process constexpr ?
			
			owner->parameters[name] = this;
			set(initial);
		}
	
		
		Parameter(Object* owner, String name, std::pair<T, Unit>initial,  Function setter = nullptr)
		: ParameterBase(owner, name, "", RangeLimit::none, setter)
		{
			// 1. iterate args
			// 2. determine their types
			// 3. do something appropriate for their given type
			// 4. can we make this whole process constexpr ?
			
			owner->parameters[name] = this;
			set(initial.first, initial.second);
		}

		
		// setter
		Parameter& operator = (T input)
		{
			set(input);
			return *this;
		}
		
		
		// setter w/ unit
		Parameter& operator = (const std::pair<T, Unit>& input)
		{
			set(input.first, input.second);
			return *this;
		}
		
		
		// setter for case when input is a generic value
		// TODO: if a value has 2 members then do we use the last one as a unit? perhaps it needs some metadata so that we know?
		Parameter& operator = (const VarBase& input)
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
	template<class T, class Dataspace_Class>
	class Parameter<T, Dataspace_Class, RangeLimit::clip> : public ParameterBase {
		T					mValue;
		Range<T>			mRange;
		Dataspace_Class		mDataspace;


		// naked values have no dataspace unit specified, so set them directly
		void set(T input)
		{
			mValue = Limit(input, mRange.first, mRange.second);
			if (mSetter)
				mSetter();
			for (auto& observer : mObservers)
				(*observer)();
		}

		
		// set values using a dataspace conversion
		void set(T input, Unit unit)
		{
			set(mDataspace(input, (uint32_t)unit));
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
		
		
		// setter w/ unit
		Parameter& operator = (const std::pair<T, Unit> input)
		{
			set(input.first, input.second);
			return *this;
		}
		
		
		// setter for case when input is a generic value
		// TODO: if a value has 2 members then do we use the last one as a unit? perhaps it needs some metadata so that we know?
		Parameter& operator = (const VarBase& input)
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
