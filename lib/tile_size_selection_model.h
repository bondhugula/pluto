#include "pluto.h"
#ifndef TILESIZESELECTIONMODEL_H
#define TILESIZESELECTIONMODEL_H

#if defined(__cplusplus)
extern "C" {
#endif

/// C wrapper around tile size selection model. Tile sizes for SCALAR dimensions
/// in the input band are set to 42, rest are inferred from the tile size
/// selection model.
void find_tile_sizes(Band *band, PlutoProg *prog, int *tile_sizes);

#if defined(__cplusplus)
}
#endif

#endif // TILESIZESELECTIONMODEL_H
