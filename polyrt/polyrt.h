#ifndef __POLYRT_HEADER_H
#define __POLYRT_HEADER_H
#include <stdlib.h>



// represents whether a distributed loop nest should be distributed block-cyclically or not
// can be used only by dynamic scheduling for now
// at the most 32 distributed loop nests for now
extern int __is_block_cyclic[32];
// represents whether a distributed loop nest should be multi-partitioned or not
// can be used only by dynamic scheduling for now
// at the most 32 distributed loop nests for now
extern int __is_multi_partitioned[32];

void polyrt_init(int nprocs);
void polyrt_init_grid_size();
void polyrt_loop_dist(int lb, int ub, int nprocs, int my_rank, int *my_start, int *my_end);
void polyrt_loop_dist_multi_part(int lb, int ub, int nprocs, int my_rank, int cur_phase, int *my_start, int *my_end);
void polyrt_multi_dim_loop_dist(int lb, int ub, int nprocs, int my_rank, int nploops, int cploop, int *my_start, int *my_end);
int polyrt_one_dim_pi(int iter, int __lb, int __ub, int nprocs);
int polyrt_one_dim_pi_multi_part(int iter, int __lb, int __ub, int nprocs, int phase);
int polyrt_two_dim_pi(int loop_id, int iter1, int lb1, int ub1, int iter2, int lb2, int ub2, int nprocs);
int polyrt_three_dim_pi(int iter1, int lb1, int ub1, int iter2, int lb2, int ub2, int iter3, int lb3, int ub3, int nprocs);
int polyrt_one_dim_pi_block_cyclic(int iter, int lb, int ub, int nprocs, unsigned int block_size);
int polyrt_two_dim_pi_block_cyclic(int iter1, int lb1, int ub1, int iter2, int lb2, int ub2, int nprocs, unsigned int block_size);
int polyrt_three_dim_pi_block_cyclic(int iter1, int lb1, int ub1, int iter2, int lb2, int ub2, int iter3, int lb3, int ub3, int nprocs, unsigned int block_size);
int polyrt_3d_flatten_dim_0(int my_rank, int nprocs);
int polyrt_3d_flatten_dim_1(int my_rank, int nprocs);
int polyrt_3d_flatten_dim_2(int my_rank, int nprocs);
int polyrt_3d_flatten_dims_except_0(int my_rank, int nprocs);
int polyrt_3d_flatten_dims_except_1(int my_rank, int nprocs);
int polyrt_3d_flatten_dims_except_2(int my_rank, int nprocs);
int polyrt_3d_project_dim_0(int sender_plane_color, int receiver_color, int nprocs);
int polyrt_3d_project_dim_1(int sender_plane_color, int receiver_color, int nprocs);
int polyrt_3d_project_dim_2(int sender_plane_color, int receiver_color, int nprocs);
void polyrt_3d_broadcast_dim_0(int receiver_color, int nprocs, int *receiver_list);
void polyrt_3d_broadcast_dim_1(int receiver_color, int nprocs, int *receiver_list);
void polyrt_3d_broadcast_dim_2(int receiver_color, int nprocs, int *receiver_list);
int polyrt_3d_is_broadcast_dim_0(int my_rank_color, int nprocs, int *receiver_list);
int polyrt_3d_is_broadcast_dim_1(int my_rank_color, int nprocs, int *receiver_list);
int polyrt_3d_is_broadcast_dim_2(int my_rank_color, int nprocs, int *receiver_list);
void polyrt_finalize();
void clear_sender_receiver_lists(int nprocs);
void add_proc_to_sender_list(int p, int my_rank);
void add_proc_to_receiver_list(int p, int my_rank);
void print_sender_list(int my_rank, int nprocs);
void print_receiver_list(int my_rank, int nprocs);
int need_to_send(int nprocs);
void *polyrt_max_alloc(void *buf, size_t size, size_t *curr_size);


struct buffer_list{
	struct buffer_list *next;

	void *buffer;

	int *ref_count;

};

typedef struct buffer_list BufferList;


void* buffer_alloc(int size, BufferList **list_ptr, int *ref_count);

void init_buffer_list(BufferList *buf_list, int size);

void buffer_free(BufferList **list_ptr);
#endif

