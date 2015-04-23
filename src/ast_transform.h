#include "pluto.h"
#include "cloog/cloog.h"

void pluto_mark_parallel(struct clast_stmt *root, const PlutoProg *prog, CloogOptions *options);
void pluto_mark_parallel_dynschedule(struct clast_stmt *root, const PlutoProg *prog, CloogOptions *options);
void pluto_mark_vector(struct clast_stmt *root, const PlutoProg *prog, CloogOptions *options);
