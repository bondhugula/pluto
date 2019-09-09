#ifndef AST_TRANSFORM_H
#define AST_TRANSFORM_H

#include "cloog/cloog.h"

typedef struct plutoProg PlutoProg;

void pluto_mark_parallel(struct clast_stmt *root, const PlutoProg *prog,
                         CloogOptions *options);
void pluto_mark_vector(struct clast_stmt *root, const PlutoProg *prog,
                       CloogOptions *options);
#endif // AST_TRANSFORM_H
