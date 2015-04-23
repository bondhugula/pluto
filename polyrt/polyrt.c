#include "polyrt.h"
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>

int debug;
#define IF_DEBUG(foo) {if (debug) { foo; } }

#define min(x,y) ((x) < (y)? (x): (y))

/* Is proc i supposed to receive anything? */
int *sender_list, *receiver_list;

int nprocs_g = -1;

// represents whether a distributed loop nest should be distributed block-cyclically or not
// can be used only by dynamic scheduling for now
// at the most 32 distributed loop nests for now
int __is_block_cyclic[32];
// represents whether a distributed loop nest should be multi-partitioned or not
// can be used only by dynamic scheduling for now
// at the most 32 distributed loop nests for now
int __is_multi_partitioned[32];

#define MAX_PROCS 128
#define MAX_GRID_DIMS 3
int grid_size[MAX_GRID_DIMS+1][MAX_PROCS+1][MAX_GRID_DIMS];

void read_grid_size(){

    FILE *file = fopen("grid.sizes", "r");

    if (!file)    return;
    int nploops, nprocs, ncloop, dim_size, ret;

    while(!feof(file)){
    	ret = fscanf(file,"n_parallel_loops=%d nprocs=%d parallel_loop_dim=%d grid_size=%d\n",&nploops, &nprocs, &ncloop, &dim_size);

    	if(ret == 4 && nprocs<=MAX_PROCS && nploops<=3 && ncloop<nploops)
			grid_size[nploops][nprocs][ncloop] = dim_size;
    }

    fclose(file);
    return;
}

void polyrt_init_grid_size(){

	int i, j, k;
	for(i=0;i<MAX_GRID_DIMS;i++)
		for(j=0;j<MAX_PROCS;j++)
			for(k=0;k<MAX_GRID_DIMS;k++)
				grid_size[i][j][k] = -1;

	//set grid size for 1d proc grid
	for(i=0;i<MAX_PROCS;i++)
		grid_size[1][i][0] = i;

	//2d proc grid, num_procs=1, row_dim = 1, col_dim = 1
	grid_size[2][1][0] = 1;
	grid_size[2][1][1] = 1;

	grid_size[2][2][0] = 2;
	grid_size[2][2][1] = 1;

	grid_size[2][3][0] = 3;
	grid_size[2][3][1] = 1;

	grid_size[2][4][0] = 2;
	grid_size[2][4][1] = 2;

	grid_size[2][5][0] = 5;
	grid_size[2][5][1] = 1;

	grid_size[2][6][0] = 3;
	grid_size[2][6][1] = 2;

	grid_size[2][7][0] = 7;
	grid_size[2][7][1] = 1;

	grid_size[2][8][0] = 4;
	grid_size[2][8][1] = 2;

	grid_size[2][9][0] = 3;
	grid_size[2][9][1] = 3;

	grid_size[2][10][0] = 5;
	grid_size[2][10][1] = 2;

	grid_size[2][11][0] = 11;
	grid_size[2][11][1] = 1;

	grid_size[2][12][0] = 4;
	grid_size[2][12][1] = 3;

	grid_size[2][13][0] = 13;
	grid_size[2][13][1] = 1;

	grid_size[2][14][0] = 7;
	grid_size[2][14][1] = 2;

	grid_size[2][15][0] = 5;
	grid_size[2][15][1] = 3;

	grid_size[2][16][0] = 8;
	grid_size[2][16][1] = 2;

	grid_size[2][17][0] = 17;
	grid_size[2][17][1] = 1;

	grid_size[2][18][0] = 6;
	grid_size[2][18][1] = 3;

	grid_size[2][19][0] = 19;
	grid_size[2][19][1] = 1;

	grid_size[2][20][0] = 5;
	grid_size[2][20][1] = 4;

	grid_size[2][21][0] = 7;
	grid_size[2][21][1] = 3;

	grid_size[2][22][0] = 11;
	grid_size[2][22][1] = 2;

	grid_size[2][23][0] = 23;
	grid_size[2][23][1] = 1;

	grid_size[2][24][0] = 6;
	grid_size[2][24][1] = 4;

	grid_size[2][25][0] = 5;
	grid_size[2][25][1] = 5;

	grid_size[2][32][0] = 8;
	grid_size[2][32][1] = 4;

	grid_size[2][36][0] = 6;
	grid_size[2][36][1] = 6;

	grid_size[2][64][0] = 8;
	grid_size[2][64][1] = 8;

	grid_size[2][128][0] = 16;
	grid_size[2][128][1] = 8;

	//3d grid, nprocs = 1, d1, d2, d3
	grid_size[3][1][0] = 1;
	grid_size[3][1][1] = 1;
	grid_size[3][1][2] = 1;

	grid_size[3][2][0] = 2;
	grid_size[3][2][1] = 1;
	grid_size[3][2][2] = 1;

	grid_size[3][3][0] = 3;
	grid_size[3][3][1] = 1;
	grid_size[3][3][2] = 1;

	grid_size[3][4][0] = 2;
	grid_size[3][4][1] = 2;
	grid_size[3][4][2] = 1;

	grid_size[3][5][0] = 5;
	grid_size[3][5][1] = 1;
	grid_size[3][5][2] = 1;

	grid_size[3][6][0] = 3;
	grid_size[3][6][1] = 2;
	grid_size[3][6][2] = 1;

	grid_size[3][7][0] = 7;
	grid_size[3][7][1] = 1;
	grid_size[3][7][2] = 1;

	grid_size[3][8][0] = 2;
	grid_size[3][8][1] = 2;
	grid_size[3][8][2] = 2;

	grid_size[3][9][0] = 3;
	grid_size[3][9][1] = 3;
	grid_size[3][9][2] = 1;

	grid_size[3][10][0] = 5;
	grid_size[3][10][1] = 2;
	grid_size[3][10][2] = 1;

	grid_size[3][11][0] = 11;
	grid_size[3][11][1] = 1;
	grid_size[3][11][2] = 1;

	grid_size[3][12][0] = 3;
	grid_size[3][12][1] = 2;
	grid_size[3][12][2] = 2;

	grid_size[3][13][0] = 13;
	grid_size[3][13][1] = 1;
	grid_size[3][13][2] = 1;

	grid_size[3][14][0] = 7;
	grid_size[3][14][1] = 2;
	grid_size[3][14][2] = 1;

	grid_size[3][15][0] = 5;
	grid_size[3][15][1] = 3;
	grid_size[3][15][2] = 1;

	grid_size[3][16][0] = 4;
	grid_size[3][16][1] = 2;
	grid_size[3][16][2] = 2;

	grid_size[3][17][0] = 17;
	grid_size[3][17][1] = 1;
	grid_size[3][17][2] = 1;

	grid_size[3][18][0] = 3;
	grid_size[3][18][1] = 3;
	grid_size[3][18][2] = 2;

	grid_size[3][19][0] = 19;
	grid_size[3][19][1] = 1;
	grid_size[3][19][2] = 1;

	grid_size[3][20][0] = 5;
	grid_size[3][20][1] = 2;
	grid_size[3][20][2] = 2;

	grid_size[3][21][0] = 7;
	grid_size[3][21][1] = 3;
	grid_size[3][21][2] = 1;

	grid_size[3][22][0] = 11;
	grid_size[3][22][1] = 2;
	grid_size[3][22][2] = 1;

	grid_size[3][23][0] = 23;
	grid_size[3][23][1] = 1;
	grid_size[3][23][2] = 1;

	grid_size[3][24][0] = 4;
	grid_size[3][24][1] = 3;
	grid_size[3][24][2] = 2;

	grid_size[3][25][0] = 5;
	grid_size[3][25][1] = 5;
	grid_size[3][25][2] = 1;

	grid_size[3][27][0] = 3;
	grid_size[3][27][1] = 3;
	grid_size[3][27][2] = 3;

	grid_size[3][32][0] = 4;
	grid_size[3][32][1] = 4;
	grid_size[3][32][2] = 2;

	grid_size[3][64][0] = 4;
	grid_size[3][64][1] = 4;
	grid_size[3][64][2] = 4;

	grid_size[3][128][0] = 8;
	grid_size[3][128][1] = 4;
	grid_size[3][128][2] = 4;

//	TODO: Fill all remaining values

	read_grid_size();
}

void polyrt_init(int nprocs)
{
    debug = 0;

    nprocs_g = nprocs;

    sender_list = (int *)malloc(nprocs*sizeof(int));
    receiver_list = (int *)malloc(nprocs*sizeof(int));

    polyrt_init_grid_size();
}


/* Performs a load-balanced block distribution of the loop. Sets my_start and
 * my_end to the lower and upper bound for a process with id 'my_rank' */
void polyrt_loop_dist(int lb, int ub, int nprocs, int my_rank, int *my_start, int *my_end)
{
    long n = ub - lb + 1;

    if (my_rank < n%nprocs)  {
        *my_start =  lb + (n/nprocs)*my_rank + my_rank;
        /* procs with id < n%nprocs get an extra iteration to distribute
         * the remainder */
        *my_end = *my_start + (n/nprocs)-1 + 1;
    }else{
        *my_start =  lb + (n/nprocs)*my_rank + n%nprocs;
        *my_end = *my_start + (n/nprocs) - 1;
    }

    IF_DEBUG(fprintf(stderr, "Proc %d:  lb: %d; ub: %d\n", my_rank, lb, ub));
    IF_DEBUG(fprintf(stderr, "Proc %d executing %d iterations from %d to %d\n", 
                my_rank, *my_end-*my_start+1, *my_start, *my_end));
}


void polyrt_loop_dist_multi_part(int lb, int ub, int nprocs, int my_rank, int cur_phase, int *my_start, int *my_end){

	polyrt_loop_dist(lb, ub, nprocs, (my_rank + cur_phase)%nprocs, my_start, my_end);

}


void polyrt_multi_dim_loop_dist(int lb, int ub, int nprocs, int my_rank, int nploops, int cploop, int *my_start, int *my_end)
{

	assert(nprocs<=MAX_PROCS);
	assert(nploops<=3);
	assert(cploop<nploops);

	int curr_grid_dim = grid_size[nploops][nprocs][cploop];
	if(curr_grid_dim <= 0){
        printf("No of dimensions in grid: %d\n", nploops);
        printf("No of processors: %d\n", nprocs);
		printf("ERROR: Couldn't determine grid size. Please specify grid size in grid.sizes file\n");
        exit(EXIT_FAILURE);
	}
	int i = 0, p = my_rank;

	for(i=0;i<cploop;i++){
		p /=  grid_size[nploops][nprocs][i];
	}

	p %= curr_grid_dim;

	polyrt_loop_dist(lb, ub, curr_grid_dim, p, my_start, my_end);

	return;
}

int polyrt_one_dim_pi(int iter, int __lb, int __ub, int nprocs){
#ifndef __DONT_USE_INV_BLOCK_DIST_FUNCTION

    /* Check: iter out of the range of the distributed set of iterations */
    if ( iter < __lb || iter > __ub ) {
        printf("iter %d lb %d ub %d nprocs %d\n", iter, __lb, __ub, nprocs);
        assert(0);
        return -1;
    }

    /* Inverse function */
    int size = __ub - __lb + 1;
    int local_index = iter - __lb;
    int limit = (size / nprocs + 1) * (size % nprocs);

    int rank;
    if ( local_index < limit ) {
        rank = local_index / (size/nprocs + 1);
    } else {
        rank = (local_index - limit) / (size/nprocs) + (size % nprocs);
    }

#ifdef DEBUG
    /* Check */
    int b,e;
    polyrt_loop_dist(  __lb, __ub, nprocs, rank, &b, &e );
    if ( iter < b || iter > e ) {
        printf("CTRL ERROR: rank=%d, b=%d, e=%d\n", rank, b, e );
    }
#endif

    return rank;

#else

    int __p;
    long compute_start, compute_end;
    long __n = __ub - __lb + 1;

    for (__p=0; __p<nprocs; __p++)    {
		if (__p < __n%nprocs)  {
			compute_start =  __lb + (__n/nprocs)*__p + __p;
			/* procs with id < __n%nprocs get an extra iteration to distributed
			 * the remainder */
			compute_end = compute_start + (__n/nprocs)-1 + 1;
		}else{
			compute_start =  __lb + (__n/nprocs)*__p + __n%nprocs;
			compute_end = compute_start + (__n/nprocs) - 1;
		}
		if (iter >= compute_start && iter <= compute_end) return __p;
    }
	//assert(0);
	return -1;

#endif
}

int polyrt_one_dim_pi_multi_part(int iter, int __lb, int __ub, int nprocs, int phase){

    return (polyrt_one_dim_pi(iter, __lb, __ub, nprocs) + phase ) % nprocs;
}

int polyrt_two_dim_pi(int loop_id, int iter1, int lb1, int ub1, int iter2, int lb2, int ub2, int nprocs){

    int p1, p2, pdim1, pdim2;

    assert(nprocs<=MAX_PROCS);

    if (__is_multi_partitioned[loop_id]) {
        int phase = polyrt_one_dim_pi(iter1, lb1, ub1, nprocs);
        return (polyrt_one_dim_pi(iter2, lb2, ub2, nprocs) + phase) % nprocs;
    }

    pdim1 = grid_size[2][nprocs][0];
	if(pdim1<=0){
        printf("No of dimensions in grid: 2\n");
        printf("No of processors: %d\n", nprocs);
		printf("ERROR: Couldn't determine grid size. Please specify grid size in grid.sizes file\n");
        exit(EXIT_FAILURE);
	}
    p1 = polyrt_one_dim_pi(iter1, lb1, ub1, pdim1);

    pdim2 = grid_size[2][nprocs][1];
	if(pdim2<=0){
        printf("No of dimensions in grid: 2\n");
        printf("No of processors: %d\n", nprocs);
		printf("ERROR: Couldn't determine grid size. Please specify grid size in grid.sizes file\n");
        exit(EXIT_FAILURE);
	}
    p2 = polyrt_one_dim_pi(iter2, lb2, ub2, pdim2);

    if(p1<0 || p2 < 0)
		return -1;

    return p1 + p2 * pdim1;
}

int polyrt_three_dim_pi(int iter1, int lb1, int ub1, int iter2, int lb2, int ub2, int iter3, int lb3, int ub3, int nprocs){

    int p1, p2, p3,  pdim1, pdim2, pdim3;

    assert(nprocs<=MAX_PROCS);

    pdim1 = grid_size[3][nprocs][0];
	if(pdim1<=0){
        printf("No of dimensions in grid: 3\n");
        printf("No of processors: %d\n", nprocs);
		printf("ERROR: Couldn't determine grid size. Please specify grid size in grid.sizes file\n");
        exit(EXIT_FAILURE);
	}
    p1 = polyrt_one_dim_pi(iter1, lb1, ub1, pdim1);

    pdim2 = grid_size[3][nprocs][1];
	if(pdim2<=0){
        printf("No of dimensions in grid: 3\n");
        printf("No of processors: %d\n", nprocs);
		printf("ERROR: Couldn't determine grid size. Please specify grid size in grid.sizes file\n");
        exit(EXIT_FAILURE);
	}
    p2 = polyrt_one_dim_pi(iter2, lb2, ub2, pdim2);

    pdim3 = grid_size[3][nprocs][2];
	if(pdim3<=0){
        printf("No of dimensions in grid: 3\n");
        printf("No of processors: %d\n", nprocs);
		printf("ERROR: Couldn't determine grid size. Please specify grid size in grid.sizes file\n");
        exit(EXIT_FAILURE);
	}
    p3 = polyrt_one_dim_pi(iter3, lb3, ub3, pdim3);

    if(p1<0 || p2 < 0 || p3 < 0)
		return -1;

    return p1 + p2 * pdim1 + p3 * pdim1 * pdim2;
}

int polyrt_one_dim_pi_block_cyclic(int iter, int lb, int ub, int nprocs, unsigned int block_size){
    /* Check: iter out of the range of the distributed set of iterations */
    if ( iter < lb || iter > ub ) {
        printf("iter %d lb %d ub %d nprocs %d block_size %u\n", iter, lb, ub, nprocs, block_size);
        assert(0);
        return -1;
    }

    int local_index = iter - lb;
    int cycle_size = nprocs * block_size;
    int index_in_cycle = local_index % cycle_size;
    int rank = index_in_cycle/block_size;

    return rank;
}

int polyrt_two_dim_pi_block_cyclic(int iter1, int lb1, int ub1, int iter2, int lb2, int ub2, int nprocs, unsigned int block_size){

    int p1, p2, pdim1, pdim2;

    assert(nprocs<=MAX_PROCS);

    pdim1 = grid_size[2][nprocs][0];
    if(pdim1<=0){
        printf("No of dimensions in grid: 2\n");
        printf("No of processors: %d\n", nprocs);
        printf("ERROR: Couldn't determine grid size. Please specify grid size in grid.sizes file\n");
        exit(EXIT_FAILURE);
    }
    p1 = polyrt_one_dim_pi_block_cyclic(iter1, lb1, ub1, pdim1, block_size);

    pdim2 = grid_size[2][nprocs][1];
    if(pdim2<=0){
        printf("No of dimensions in grid: 2\n");
        printf("No of processors: %d\n", nprocs);
        printf("ERROR: Couldn't determine grid size. Please specify grid size in grid.sizes file\n");
        exit(EXIT_FAILURE);
    }
    p2 = polyrt_one_dim_pi_block_cyclic(iter2, lb2, ub2, pdim2, block_size);

    if(p1<0 || p2 < 0)
        return -1;

    return p1 + p2 * pdim1;
}

int polyrt_three_dim_pi_block_cyclic(int iter1, int lb1, int ub1, int iter2, int lb2, int ub2, int iter3, int lb3, int ub3, int nprocs, unsigned int block_size){

    int p1, p2, p3,  pdim1, pdim2, pdim3;

    assert(nprocs<=MAX_PROCS);

    pdim1 = grid_size[3][nprocs][0];
	if(pdim1<=0){
        printf("No of dimensions in grid: 3\n");
        printf("No of processors: %d\n", nprocs);
		printf("ERROR: Couldn't determine grid size. Please specify grid size in grid.sizes file\n");
        exit(EXIT_FAILURE);
	}
    p1 = polyrt_one_dim_pi_block_cyclic(iter1, lb1, ub1, pdim1, block_size);

    pdim2 = grid_size[3][nprocs][1];
	if(pdim2<=0){
        printf("No of dimensions in grid: 3\n");
        printf("No of processors: %d\n", nprocs);
		printf("ERROR: Couldn't determine grid size. Please specify grid size in grid.sizes file\n");
        exit(EXIT_FAILURE);
	}
    p2 = polyrt_one_dim_pi_block_cyclic(iter2, lb2, ub2, pdim2, block_size);

    pdim3 = grid_size[3][nprocs][2];
	if(pdim3<=0){
        printf("No of dimensions in grid: 3\n");
        printf("No of processors: %d\n", nprocs);
		printf("ERROR: Couldn't determine grid size. Please specify grid size in grid.sizes file\n");
        exit(EXIT_FAILURE);
	}
    p3 = polyrt_one_dim_pi_block_cyclic(iter3, lb3, ub3, pdim3, block_size);

    if(p1<0 || p2 < 0 || p3 < 0)
		return -1;

    return p1 + p2 * pdim1 + p3 * pdim1 * pdim2;
}

int polyrt_3d_flatten_dim_0(int my_rank, int nprocs) {
    int pdim1 = grid_size[3][nprocs][0];
    return (int)(my_rank / pdim1);
}

int polyrt_3d_flatten_dim_1(int my_rank, int nprocs) {
    int pdim1 = grid_size[3][nprocs][0];
    int pdim2 = grid_size[3][nprocs][1];
    return (my_rank % pdim1) + (pdim1 * (int)(my_rank / (pdim1*pdim2)));
}

int polyrt_3d_flatten_dim_2(int my_rank, int nprocs) {
    int pdim1 = grid_size[3][nprocs][0];
    int pdim2 = grid_size[3][nprocs][1];
    return my_rank % (pdim1*pdim2);
}

int polyrt_3d_flatten_dims_except_0(int my_rank, int nprocs) {
    int pdim1 = grid_size[3][nprocs][0];
    return my_rank % pdim1;
}

int polyrt_3d_flatten_dims_except_1(int my_rank, int nprocs) {
    int pdim1 = grid_size[3][nprocs][0];
    int pdim2 = grid_size[3][nprocs][1];
    return ((int)(my_rank / pdim1)) % pdim2;
}

int polyrt_3d_flatten_dims_except_2(int my_rank, int nprocs) {
    int pdim1 = grid_size[3][nprocs][0];
    int pdim2 = grid_size[3][nprocs][1];
    return (int)(my_rank / (pdim1 * pdim2));
}

int polyrt_3d_project_dim_0(int sender_plane_color, int receiver_color, int nprocs) {
    int pdim1 = grid_size[3][nprocs][0];
    return (int)(receiver_color * pdim1) + sender_plane_color;
}

int polyrt_3d_project_dim_1(int sender_plane_color, int receiver_color, int nprocs) {
    int pdim1 = grid_size[3][nprocs][0];
    int pdim2 = grid_size[3][nprocs][1];
    return (int)((receiver_color/pdim1) * (pdim1*pdim2)) + (receiver_color%pdim1) + (sender_plane_color*pdim1);
}

int polyrt_3d_project_dim_2(int sender_plane_color, int receiver_color, int nprocs) {
    int pdim1 = grid_size[3][nprocs][0];
    int pdim2 = grid_size[3][nprocs][1];
    return (int)(receiver_color + (sender_plane_color*pdim1*pdim2));
}

void polyrt_3d_broadcast_dim_0(int receiver_color, int nprocs, int *receiver_list) {
    int pdim1 = grid_size[3][nprocs][0];
    int i, rank;
    for (i=0; i<pdim1; i++) {
        rank = (receiver_color * pdim1) + i;
        receiver_list[rank] = 1;
    }
}

void polyrt_3d_broadcast_dim_1(int receiver_color, int nprocs, int *receiver_list) {
    int pdim1 = grid_size[3][nprocs][0];
    int pdim2 = grid_size[3][nprocs][1];
    int i, rank;
    for (i=0; i<pdim2; i++) {
        rank = ((receiver_color/pdim1) * (pdim1*pdim2)) + (receiver_color%pdim1) + (i*pdim1);
        receiver_list[rank] = 1;
    }
}

void polyrt_3d_broadcast_dim_2(int receiver_color, int nprocs, int *receiver_list) {
    int pdim1 = grid_size[3][nprocs][0];
    int pdim2 = grid_size[3][nprocs][1];
    int pdim3 = grid_size[3][nprocs][2];
    int i, rank;
    for (i=0; i<pdim3; i++) {
        rank = receiver_color + (i*pdim1*pdim2);
        receiver_list[rank] = 1;
    }
}

int polyrt_3d_is_broadcast_dim_0(int my_rank_color, int nprocs, int *receiver_list) {
    int pdim1 = grid_size[3][nprocs][0];
    int i, rank;
    int broadcast = 0;
    for (i=0; i<pdim1; i++) {
        rank = (my_rank_color * pdim1) + i;
        if (receiver_list[rank] == 1) {
            broadcast = 1;
            break;
        }
    }
    if (broadcast) {
        for (i=0; i<pdim1; i++) {
            rank = (my_rank_color * pdim1) + i;
            receiver_list[rank] = 0;
        }
    }
    return broadcast;
}

int polyrt_3d_is_broadcast_dim_1(int my_rank_color, int nprocs, int *receiver_list) {
    int pdim1 = grid_size[3][nprocs][0];
    int pdim2 = grid_size[3][nprocs][1];
    int i, rank;
    int broadcast = 0;
    for (i=0; i<pdim2; i++) {
        rank = ((my_rank_color/pdim1) * (pdim1*pdim2)) + (my_rank_color%pdim1) + (i*pdim1);
        if (receiver_list[rank] == 1) {
            broadcast = 1;
            break;
        }
    }
    if (broadcast) {
        for (i=0; i<pdim2; i++) {
            rank = ((my_rank_color/pdim1) * (pdim1*pdim2)) + (my_rank_color%pdim1) + (i*pdim1);
            receiver_list[rank] = 0;
        }
    }
    return broadcast;
}

int polyrt_3d_is_broadcast_dim_2(int my_rank_color, int nprocs, int *receiver_list) {
    int pdim1 = grid_size[3][nprocs][0];
    int pdim2 = grid_size[3][nprocs][1];
    int pdim3 = grid_size[3][nprocs][2];
    int i, rank;
    int broadcast = 0;
    for (i=0; i<pdim3; i++) {
        rank = my_rank_color + (i*pdim1*pdim2);
        if (receiver_list[rank] == 1) {
            broadcast = 1;
            break;
        }
    }
    if (broadcast) {
        for (i=0; i<pdim3; i++) {
            rank = my_rank_color + (i*pdim1*pdim2);
            receiver_list[rank] = 0;
        }
    }
    return broadcast;
}

void polyrt_finalize()
{
    free(sender_list);
    free(receiver_list);
}

void clear_sender_receiver_lists(int nprocs)
{
    int i;
    for (i=0; i<nprocs; i++) {
        sender_list[i] = 0;
        receiver_list[i] = 0;
    }
}       

void add_proc_to_sender_list(int p, int my_rank)
{
    assert(p >= 0  && p <= nprocs_g-1); 
    if (p != my_rank)   {
        sender_list[p] = 1;
    }
}


void add_proc_to_receiver_list(int p, int my_rank)
{
    // fprintf(stderr, "%d\n", p );
    assert(p >= 0 && p <= nprocs_g-1); 
    if (p != my_rank)   {
        receiver_list[p] = 1;
    }
}


void print_sender_list(int my_rank, int nprocs)
{   
    int i;
    fprintf(stderr, "procs: ");
    for (i=0; i<nprocs; i++) {
        if (sender_list[i]) {                    
            fprintf(stderr, "%d ", i);                    
        }
    }
    fprintf(stderr, "\n");
}       


void print_receiver_list(int my_rank, int nprocs)
{   
    int i;
    fprintf(stderr, "procs: ");
    for (i=0; i<nprocs; i++) {
        if (receiver_list[i]) {                    
            fprintf(stderr, "%d ", i);                    
        }
    }
    fprintf(stderr, "\n");
}       


/* Does the calling have anything in its receiver list? */
int need_to_send(int nprocs)
{
    int sum = 0, i;

    /* reciver_list[my_rank] would already be zero from the way receiver_list
     * is created from sigma */
    for (i=0; i<nprocs; i++) {
        sum += receiver_list[i];
    }

    return (sum==0)? 0: 1;
}


/* Reallocates if needed; update prev_size; buf should be initialized
 * to NULL or should point to dynamically allocated storage  */
void *polyrt_max_alloc(void *buf, size_t size, size_t *curr_size)
{
    if (size > *curr_size) {
        buf = realloc(buf, size);
        *curr_size = size;
        if (buf == NULL) {
            printf("[polyrt] Memory allocation failure: asked for %lu bytes\n", size);
            assert(0);
        }
    }
    return buf;
}

/* buffer management code for data tile allocation
 */
void* buffer_alloc(int size, BufferList **list_ptr, int *ref_count){
	assert(list_ptr != NULL);
	BufferList *list = *list_ptr;

	if(list == NULL){
		list = (BufferList*)malloc(sizeof(BufferList));
		list->next = NULL;
		list->buffer = malloc(size);
        assert(list->buffer != NULL);
		list->ref_count = ref_count;

		return list->buffer;
	}

	/* Check to see if there is any buffer with ref_count as zero
	 * if so return that buffer else alloc a new buffer
	 */

	BufferList *prev, *curr;
	curr = list;
	while(curr != NULL){
		if(curr->ref_count == 0){
			curr->ref_count = ref_count;
			return curr->buffer;
		}

		prev = curr;
		curr = curr->next;
	}

	assert(prev != NULL);
	BufferList *new_buff = (BufferList*)malloc(sizeof(BufferList));
	prev->next = new_buff;
	new_buff->next = NULL;
	new_buff->buffer = malloc(size);
	new_buff->ref_count = ref_count;

	return new_buff->buffer;

}

void buffer_free(BufferList **list_ptr){
	assert(list_ptr != NULL);
	BufferList *list = *list_ptr;

	if(list == NULL) return;

	BufferList *next;
	while(list != NULL){
		free(list->buffer);
		next = list->next;
		free(list);
		list = next;
	}

	return;
}
