void task_computation(int i, int j, int k, int l);
void generate_dag(void);
void dag_execute();
void init_dag(int nVertices);
void dag_add_vertex(int i, int j, int k, int l, double vertexWt, int *vertexId);
void dag_add_edge(int i1, int j1, int k1, int l1, int i2, int j2, int k2, int l2, double reuseVolume);
void update_number_vertices(void);
void free_dag(void);
void print_dag(void);

#ifdef __STATIC_SCHEDULE__
void init_schedule();
void free_schedule(void);
void schedule_dag_using_mcp(void);
#endif

#ifndef __USE_DAG_SCHEDULER_AS_LIB__
#include "scheduler.c"
#endif

