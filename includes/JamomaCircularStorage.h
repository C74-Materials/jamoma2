//
//  JamomaCircularStorage.h
//  JamomaWad
//
//  Created by Timothy Place and Nathan Wolek on 7/21/15.
//  Copyright (c) 2015 jamoma. All rights reserved.
//

#pragma once

#include <thread>
#include <cassert>

namespace Jamoma {

	
	/**	A generic circular buffer designed specifically for access from a single thread.
	 */
	template <class T>
	class CircularStorage {
		
		std::vector<T>		mItems;										///< storage for the circular buffer's data
		std::size_t			mIndex = { 0 };								///< location of the record head
		std::thread::id		mThread = { std::this_thread::get_id() };	///< used to ensure we don't access unsafely from multiple threads
	
	public:
		/** Constructor specifies a fixed size for the container.
			If you want a different size, create a new container and dispose of the one you don't want.
		 */
		CircularStorage(std::size_t itemCount)
		: mItems(itemCount)
		{}
		
		
		/** Write a block of things into the container. 
			@param	newInput	A block of things to add. May not be larger than the size of the buffer.
		 */
		void write(const std::vector<T>& newInput)
		{
			assert(std::this_thread::get_id() == mThread);
			assert(newInput.size() <= mItems.size());
			
			std::size_t count = newInput.size();
			std::size_t roomRemaining = mItems.size() - mIndex;
			bool		wrap = false;
			
			if (count > roomRemaining) {
				count = roomRemaining;
				wrap = true;
			}
			
			std::copy_n(newInput.begin(), count, mItems.begin()+mIndex);
			mIndex += count;
			
			if (wrap) {
				std::size_t offset = count;

				count = newInput.size() - offset;
				std::copy_n(newInput.begin()+offset, count, mItems.begin());
				mIndex = count;
			}
		}
		
		
		/** Read a block of things out from the container.
			These will be the N most recent items added to the history.
			@param	output	A place to write the block of things from the buffer. 
							The size of this buffer will determine the number of things to request. 
							May not be larger than the size of the buffer.
		 */
		void read(std::vector<T>& output)
		{
			assert(std::this_thread::get_id() == mThread);
			assert(mItems.size() >= output.size());
			
			long	count = output.size();
			long	start = mIndex - count;
			bool	wrap = false;
			
			if (start<0) {
				count = -start;
				start = mItems.size() + start;
				wrap = true;
			}
			
			std::copy_n(mItems.begin()+start, count, output.begin());
			
			if (wrap) {
				std::size_t offset = count;
				
				count = output.size() - offset;
				std::copy_n(mItems.begin(), count, output.begin()+offset);
			}
		}
	};
	
	
} // namespace Jamoma
