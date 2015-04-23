/*
 * buffer_manger.h
 *
 *  Created on: 07-Oct-2013
 *      Author: chandan
 */

#ifndef BUFFER_MANGER_H_
#define BUFFER_MANGER_H_

#include <stdlib.h>
#include "tbb/atomic.h"
#include "tbb/concurrent_queue.h"

struct buffer_manager{

	//ref count of buffers
	tbb::atomic<int> *ref_count;

	//atomic buffer ptrs
	tbb::atomic<double *> *data_tile_buffers;

	//queue to hold all free buffers whoes ref count is zero
	tbb::concurrent_queue<double *> *free_buffers;
};

typedef struct buffer_manager BufferManager;


BufferManager *init_buffer_manager( int size);

void *buffer_alloc_atomic(int index, int size, BufferManager *buffer_manag);

void increment_ref_count(BufferManager *buf_mang, int index);

void decrement_ref_count(BufferManager *buf_mang, int index);

void free_buffer_mang(BufferManager *buf_mang);

#endif /* BUFFER_MANGER_H_ */
