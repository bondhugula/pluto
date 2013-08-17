#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

int debug;
#define IF_DEBUG(foo) {if (debug) { foo; } }

#define min(x,y) ((x) < (y)? (x): (y))

/* Is proc i supposed to receive anything? */
int *sender_list, *receiver_list;

long *compute_start, *compute_end;

int nprocs_g = -1;

void polyrt_init(int nprocs)
{
    debug = 0;

    nprocs_g = nprocs;

    compute_start = malloc(nprocs*sizeof(long));
    compute_end = malloc(nprocs*sizeof(long));

    sender_list = malloc(nprocs*sizeof(int));
    receiver_list = malloc(nprocs*sizeof(int));
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


void polyrt_finalize()
{
    free(compute_start);
    free(compute_end);

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
