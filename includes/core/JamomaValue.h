/** @file
	
	@ingroup jamoma2
	
	@brief Represents an arbitrary type of value.
	
	@details Defines the primary interface used to pass values to and from methods in Jamoma. 
	
	@author Timothy Place, Nathan Wolek
	
	@copyright Copyright © 2015 by Jamoma authors and contributors @n
	This code is licensed under the terms of the "BSD 3-Clause License" @n
	https://github.com/jamoma/jamoma2/blob/master/LICENSE.md @n
 */

#pragma once

//#include "eggs/variant.hpp"

namespace Jamoma {
	
//#define VARIANT_BASED_VALUE
#ifdef VARIANT_BASED_VALUE
	
	
	
// how to make a template class...
//	template <>
//	class variant<>
//	{

	/*

		This implementation is uses Egg.Variant ( https://github.com/eggs-cpp/variant ) which internally is a discriminated union.
		The downside to this implementation is that all types that are going to be used by the union must be declared in advance.
		The upside is that all instances of the Value use C++ value semantics and are created/destroyed on the stack.
		Specific to the Eggs.Variant implementation is that much of this is a constexpr and calls can be executed at compile time.
	 
		@see http://stackoverflow.com/questions/18856824/ad-hoc-polymorphism-and-heterogeneous-containers-with-value-semantics
		@see http://boost.2283326.n4.nabble.com/Another-variant-type-was-peer-review-queue-tardiness-was-Cleaning-out-the-Boost-review-queue-Review--td4674046.html
		@see http://talesofcpp.fusionfenix.com/post-20/eggs.variant---part-ii-the-constexpr-experience
	 */
	
	using namespace eggs;


	using Value		= eggs::variant<int, double, Function>;
	using Vector	= std::vector<Value>;
	using Var		= eggs::variant<Value, Vector>;
	
	
	// Actually, a vector won't work in a variant -- maybe a shared_ptr to a vector would?  sigh...
	// One of the explicit limitations of a variant is that it doesn't support types with heap-allocated memory.
	
	
	
	

	//	template<>
//	class Value : public eggs::variant<int, double, Function> {
	

	
	template <typename ...Ts>
	class Variant : public eggs::variant<Ts...> {
	public:
		// TODO: is there a downside to accessing values like this?
		// It seems strange it would not have been implemented in eggs already?
		template<typename U>
		constexpr operator U() const {
//			U temp;
//			variants::get<U>(temp);
//			return temp;
			return variants::get<U>(*this);
		}
		
		
		constexpr Variant() noexcept
		{}
		
/*
		template<typename U>
		Variant(U arg)
		: eggs::variant<Ts...>(arg)
		{}
*/
		
		
		/**	Assignment operator.
			@param	value	The value from which to assign.
		 Type conversions will be done prior to assignment.
		 */
//		template<typename U>
//		Variant& operator = (const U& rhs) {
//			return eggs::variant<Ts...>::operator=(rhs);
//		}

		
		template<typename U>
		constexpr Variant& operator = (U rhs) {
			return eggs::variant<Ts...>::operator=(rhs);
		}

	
	};
	
	
	
	// this is how we typically pass values around generically in Jamoma
	// vector's storage is on the heap, so it's fine to include here without ballooning the mem footprint
//	using Var = eggs::variant<Value, Vector>;
/*	template<>
	class Var : public eggs::variant<Value, Vector> {
		
		// TODO: is there a downside to accessing values like this?
		// It seems strange it would not have been implemented in eggs already?
		template<typename U>
		operator U() const {
			U temp;
			variants::get<U>(temp);
			return temp;
		}
		
	};
*/	// TODO: the above is duplicated -- so maybe it should be one class...  a Jamoma::Variant
	
	
#else
	
	
	
	template <class T> class Value;
	
	/** The Value class provides a way to represent an arbitrary type of value.
		A limited amount of memory is allocated on the stack for storage of the value.
		Specializations of more extravagant types my alloc memory on the heap and store it here if needed.
	 
		This implementation uses a generic ValueBase class with templated specialization for specific data types and conversions.
		It is important that the Value<> specializations do *not* have any additional member variables, as they will be stripped when copying/moving through ValueBase instances
	 */
	class ValueBase {
		template <class T> friend class Value;
		
		
		/**	The type of the data represented by the Value.
		 */
		const std::type_info* mType = {&typeid(nullptr)};
		
		
		/**	The data, or a pointer to the data if a specialization uses it in that way.
		 */
		char mData[8] = {0,0,0,0,0,0,0,0};
		
		
	public:
		/**	Create a Value with an initial value.
			@param initial	The initial value for the Value. The type of the Value will be taken from the initial value's type.
		 */
		template <typename U>
		ValueBase(const U initial)
		: mType { &typeid(U) }
		{
			Value<U>* specialization = (Value<U>*)(this);
			*specialization = initial;
		}
		
		
		/** Default constructor with no specialization. 
			Intended for internal use only.
		 */
		ValueBase()
		{}
		
		
		/**	Destroy a Value.
		 */
		virtual ~ValueBase()
		{}
		
		
		/**	Assign a value to an Value.
			TODO: If the specialization changes then we may have resources to release!
			@param	value	The value to be assigned.
			@discussion		The actual implementation of the assignment is performed by the specialized Value, not by the ValueBase itself.
		 */
		template<typename U>
		ValueBase& operator = (const U value)
		{
			Value<U>* specialization = dynamic_cast<Value<U>*>(this);
			*specialization = value;
			return *this;
		}
		
		
		/**	Fetch the Value.
			The requested type *must* match the type of the Value.
			@return			The Value
			@discussion		The actual implementation of the assignment is performed by the specialized Value, not by the ValueBase itself.
		 */
		template<typename U>
		operator U() const
		{
			assert(typeid(U) == *mType);
			const Value<U>* specialization = (const Value<U>*)this;
			return *specialization;
		}
	};
	
	
	/**	Value specialization for most basic types of data.
		This specialization assumes that the data is a stored directly on the stack and the type is not some sort of pointer or reference.
		TODO: can we use some type deduction to perform some compile-time assertions?
	 */
	template <class T>
	class Value : public ValueBase {
	public:
		/** Create an uninitialized Value.
			Initializing the value at construction is strongly preferred.
		 */
		Value()
		{
			mType = &typeid(T);
		}
		
		/** Create an initialized Value.
			@param	initial	The value with which to initialize.
		 */
		Value(T initial)
		{
			*this = initial;
		}
		

		Value(const Value&) = default;	// using the default copy constructor

		
		/**	Assignment operator from another Value.
			@param	value	The value from which to assign.
		 */
		Value& operator = (const ValueBase& value)
		{
			mType = &typeid(T);
			mData = value.mData;
			return *this;
		}
		
		
		/**	Assignment operator.
			@param	value	The value from which to assign. 
							Type conversions will be done prior to assignment.
		 */
		Value& operator = (const T value)
		{
			mType = &typeid(T);
			memcpy(mData, *(char**)&value, sizeof(char)*8);
			return *this;
		}
		
		
		/**	Casting operator to fetch the value.
			@return	The value.
		 */
		operator T() const
		{
			T value;
			memcpy(*(char**)&value, mData, sizeof(char)*8);
			return value;
		}
	};
	
	
	/**	Value specialization for the Function type
	 */
	template <>
	class Value<Function> : public ValueBase {
	public:
		/** Create an uninitialized Function Value.
			Initializing the value at construction is preferred.
		 */
		Value()
		{
			mType = &typeid(Function);
			memcpy(mData, *(char**) new Function, sizeof(char)*8);
		}
		
		
		/** Create a Value initialized with a Function.
			@param	initial		The Function with which to initialize.
		 */
		Value(Function initial)
		{
			memcpy(mData, *(char**) new Function, sizeof(char)*8);
			*this = initial;
		}
		
		
		/** Copy constructor */
		Value(const Value& other)
		{
			memcpy(mData, *(char**) new Function, sizeof(char)*8);
			*this = other;
		}

		
		/** Destructor */
		virtual ~Value()
		{
			delete (Function*)mData;
		}
		
		
		/**	Assignment operator.
			@param	value	The Function to be assigned.
		 */
		Value& operator = (const Function value)
		{
			mType = &typeid(Function);
			*(Function*)mData = value;
			return *this;
		}
		
		
		/**	Casting operator to fetch the Function.
			@return	The Function.
		 */
		operator Function() const
		{
			return *(Function*)mData;
		}
	};
	
	using Var = ValueBase;
	
#endif

} // namespace Jamoma
