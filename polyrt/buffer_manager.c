/*
 * buffer_manger.c
 *
 *  Created on: 07-Oct-2013
 *      Author: chandan
 */


#include "buffer_manager.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>

BufferManager *init_buffer_manager(int size){

	assert(size > 0);

	BufferManager *buf_mang = (BufferManager *)malloc(sizeof(BufferManager));

	buf_mang->data_tile_buffers = (tbb::atomic<double*> *)malloc(size * sizeof(tbb::atomic<double*>));
	buf_mang->ref_count = (tbb::atomic<int> *)malloc(size * sizeof(tbb::atomic<int>));

	int i;
	for(i=0; i<size;i++)
		buf_mang->ref_count[i] = 0;

	for(i=0; i<size;i++)
		buf_mang->data_tile_buffers[i] = NULL;

	buf_mang->free_buffers = new tbb::concurrent_queue<double *>();

	return buf_mang;
}

void *buffer_alloc_atomic(int index, int size, BufferManager *buffer_manag){

	assert(buffer_manag != NULL);
	assert(index >= 0);
	assert(size > 0);


	if(buffer_manag->data_tile_buffers[index] != NULL)
		return buffer_manag->data_tile_buffers[index];

	double *ptr;
	bool found = false;

    if(!buffer_manag->free_buffers->empty())
        found = buffer_manag->free_buffers->try_pop(ptr);

   if(!found)
		ptr = (double *)malloc(size);

	assert(ptr != NULL);

	double *old_val = buffer_manag->data_tile_buffers[index].compare_and_swap(ptr, NULL);

	if(old_val != NULL){
		buffer_manag->free_buffers->push(ptr);
	}

	assert(buffer_manag->data_tile_buffers[index] != NULL);

	return buffer_manag->data_tile_buffers[index];
}

void increment_ref_count(BufferManager *buf_mang, int index){

	assert(buf_mang != NULL);

	buf_mang->ref_count[index]++;

	return;
}

void decrement_ref_count(BufferManager *buf_mang, int index){

	assert(buf_mang != NULL);

    buf_mang->ref_count[index]--;

	if(buf_mang->ref_count[index] == 0){

		double *ptr = buf_mang->data_tile_buffers[index];
		assert(ptr != NULL);
		buf_mang->free_buffers->push(ptr);
		buf_mang->data_tile_buffers[index].compare_and_swap(NULL, ptr);
	}

	return;
}

void free_buffer_mang(BufferManager *buf_mang){

	assert(buf_mang != NULL);

	free(buf_mang->data_tile_buffers);
	free(buf_mang->ref_count);

	double *ptr;
	while(buf_mang->free_buffers->try_pop(ptr)){
		free(ptr);
	}

	delete buf_mang->free_buffers;

	free(buf_mang);

	return;
}

