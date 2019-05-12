#include "cloog/cloog.h"
#include "pluto.h"

void pluto_mark_parallel(struct clast_stmt *root, const PlutoProg *prog,
                         CloogOptions *options);
void pluto_mark_vector(struct clast_stmt *root, const PlutoProg *prog,
                       CloogOptions *options);
