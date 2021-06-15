#include "tile_size_selection_model.h"
#include "math_support.h"
#include "pluto/matrix.h"
#include "pluto/pluto.h"
#include "post_transform.h"
#include "program.h"
#include <algorithm>
#include <cassert>
#include <math.h>
#include <set>
#include <sys/time.h>
#include <vector>

class TileSizeSelectionModel {
private:
  std::vector<float> dimReuse;
  std::vector<bool> vectorizableLoops;
  unsigned vectorDim;
  unsigned cache_size_per_core;
  unsigned available_cache_size_per_core;
  unsigned num_data_elements_in_cache;
  std::vector<double> reuse_polynomial;
  void add_tile_size_constraints_for_diamond_tiling();
  unsigned compute_tile_footprint_coefficient_for_access(PlutoAccess *access,
                                                         Band *band,
                                                         PlutoProg *prog,
                                                         double &coefficient);
  unsigned tile_size_for_vector_dimension;
  unsigned unroll_jam_factor;
  unsigned tile_size_for_parallel_dimension;
  unsigned par_loop_depth;
  bool has_pipeline_parallelism;

public:
  void construct_expression_for_tile_volume(Band *band, PlutoProg *prog);
  TileSizeSelectionModel(const std::vector<float> &_dimReuse,
                         const std::vector<bool> &_vectorizableLoops,
                         unsigned cache_size, unsigned data_element_size,
                         unsigned ufactor, unsigned _par_loop_depth,
                         bool _has_pipeline_parallelism)
      : dimReuse(_dimReuse), vectorizableLoops(_vectorizableLoops),
        cache_size_per_core(cache_size), unroll_jam_factor(ufactor),
        par_loop_depth(_par_loop_depth),
        has_pipeline_parallelism(_has_pipeline_parallelism) {
    num_data_elements_in_cache = cache_size_per_core / data_element_size;
    tile_size_for_vector_dimension = 512;
    tile_size_for_parallel_dimension = 32;
    vectorDim = _vectorizableLoops.size();
    for (unsigned i = 0; i < _vectorizableLoops.size(); i++) {
      if (!_vectorizableLoops[i])
        continue;
      if (vectorDim == _vectorizableLoops.size()) {
        vectorDim = i;
        continue;
      }

      if (dimReuse[i] >= dimReuse[vectorDim])
        vectorDim = i;
    }
  }
  unsigned solve_reuse_expression();
  std::vector<unsigned> infer_tile_sizes(unsigned root);
};

/// TODO: have a single central implementation
static double rtclock() {
  struct timeval Tp;
  int stat = gettimeofday(&Tp, NULL);
  if (stat != 0)
    printf("Error return from gettimeofday: %d", stat);
  return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}

/// Solves the quadratic equation whose coefficients are given in the input
/// vector. The coefficient of i^th degree is present at the i^th location of
/// the vector.
static double solve_quadratic_equation(const std::vector<double> &equation) {
  assert(equation.size() == 3);
  double a = equation[2];
  double b = equation[1];
  double c = equation[0];

  double discriminant = (b * b) - 4 * a * c;
  assert(discriminant >= 0);
  double tmp = sqrt(discriminant);
  double root1 = (-b + tmp) / (2 * a);
  double root2 = (-b - tmp) / (2 * a);
  return (root1 > root2) ? root1 : root2;
}

/// Returns the approximation of the input number 'num' to the nearest multiple
/// of the input 'factor'.
static unsigned approximate_to_nearest_multiple(unsigned num, unsigned factor) {
  assert(factor > 0);
  if (factor == 1)
    return num;
  unsigned r = num % factor;
  if (r == 0)
    return num;
  // round below
  if (2 * r < factor)
    return num - r;
  // Round up
  return num - r + factor;
}

/// Solves the reuse polynomial and returns the maximum positive root. The
/// current implementation gives precise results when the degree of the reuse
/// polynomial is less than 3. When degree > 3, the reuse polynomaial is
/// approximated to the highest power. This can be very imprecise, but we did
/// not see any cases where the degree was larger than 2. Integrate with an
/// external library to solve higher degee polynomials.
unsigned TileSizeSelectionModel::solve_reuse_expression() {

  // Just a constant. return anything
  if (reuse_polynomial.size() < 2)
    return 0;
  assert(reuse_polynomial.size() >= 2 and reuse_polynomial[0] < 0.0f);

  // If quadratic, find the exact roots. Else approximate.
  if (reuse_polynomial.size() == 3) {
    double max_root = solve_quadratic_equation(reuse_polynomial);
    assert(max_root >= 0 && "Negative root for reuse polynomial");
    return (unsigned)max_root;
  }
  double base = (-reuse_polynomial[0]) / reuse_polynomial.back();
  double exponent = (double)1 / (double)(reuse_polynomial.size() - 1);
  double root = pow(base, exponent);
  return (unsigned)root;
}

/// Rounds up the tile size to the nearest multiple of unroll factor.
static void round_tile_sizes_wrt_unroll_jam_factor(unsigned &tile_size,
                                                   unsigned ufactor) {
  // TODO: Check if the loop is suitable for unroll jam ?
  if (tile_size < ufactor)
    tile_size = ufactor;
  else
    tile_size = approximate_to_nearest_multiple(tile_size, ufactor);
}

/// Infers the tile sizes of each dimension, based on the solution to the reuse
/// polynomial.
std::vector<unsigned> TileSizeSelectionModel::infer_tile_sizes(unsigned root) {
  std::vector<unsigned> tile_sizes(dimReuse.size());
  float min_dim_reuse = *std::max_element(dimReuse.begin(), dimReuse.end());

  for (unsigned i = 0; i < dimReuse.size(); i++)
    if (dimReuse[i] > 0.0f && !vectorDim && dimReuse[i] < min_dim_reuse) {
      min_dim_reuse = dimReuse[i];
    }

  unsigned min_reuse_dim =
      find(dimReuse.begin(), dimReuse.end(), min_dim_reuse) - dimReuse.begin();
  assert(min_reuse_dim != dimReuse.size());
  // This can happen only when the vectorizable dimension also has a dimensional
  // reuse the same min_dim_reuse.
  if (min_reuse_dim == vectorDim || min_reuse_dim == par_loop_depth) {
    for (unsigned i = 0; i < dimReuse.size(); i++) {
      if (i == vectorDim || i == par_loop_depth)
        continue;
      if (dimReuse[i] == min_dim_reuse) {
        min_reuse_dim = i;
        break;
      }
    }
  }

  tile_sizes[min_reuse_dim] = dimReuse[min_reuse_dim] * root;
  round_tile_sizes_wrt_unroll_jam_factor(tile_sizes[min_reuse_dim],
                                         unroll_jam_factor);

  unsigned min_tile_size = tile_sizes[min_reuse_dim];

  // Compute tile size for the dimension with the smallest dimensional reuse.
  for (unsigned i = 0; i < dimReuse.size(); i++) {
    if (has_pipeline_parallelism && (i == 0 || i == 1)) {
      tile_sizes[i] = 32;
      continue;
    }
    if (i == vectorDim) {
      tile_sizes[i] = tile_size_for_vector_dimension;
      continue;
    }
    if (i == par_loop_depth) {
      tile_sizes[i] = tile_size_for_parallel_dimension;
      continue;
    }
    if (i == min_reuse_dim)
      continue;
    float dim_reuse_ratio = dimReuse[i] / min_dim_reuse;
    tile_sizes[i] = dim_reuse_ratio * min_tile_size;
    round_tile_sizes_wrt_unroll_jam_factor(tile_sizes[i], unroll_jam_factor);
  }
  return tile_sizes;
}

/// Returns distinct accesses in the band as a set.
std::set<PlutoAccess *> get_distinct_accesses_in_band(Band *band,
                                                      PlutoProg *prog) {
  unsigned nstmts = band->loop->nstmts;
  std::set<PlutoAccess *> distinctAccesses;
  Stmt **stmts = band->loop->stmts;
  unsigned num_unique_accesses;
  PlutoAccess **accesses =
      get_unique_accesses_in_stmts(stmts, nstmts, prog, &num_unique_accesses);

  for (unsigned i = 0; i < num_unique_accesses; i++)
    distinctAccesses.insert(accesses[i]);
  free(accesses);

  return distinctAccesses;
}

/// Returns a new matrix whose rows correspond to linearly independent rows of
/// the input matrix.
PlutoMatrix *get_linearly_independent_rows_in_mat(PlutoMatrix *mat) {
  PlutoMatrix *newMat = pluto_matrix_dup(mat);
  if (pluto_matrix_get_rank(mat) == mat->nrows)
    return newMat;
  unsigned num_independent_rows = 0;
  for (unsigned i = 0; i < mat->nrows; i++) {
    newMat->nrows = num_independent_rows + 1;
    if (pluto_matrix_get_rank(newMat) == num_independent_rows + 1) {
      num_independent_rows++;
      continue;
    }
    pluto_matrix_remove_row(newMat, num_independent_rows);
  }
  assert(newMat->nrows == num_independent_rows);
  return newMat;
}

/// The routines computes the tile footprint of the input access in the band. It
/// returns the corresponding coefficient and the degree of the variable in the
/// reuse polynomial.
unsigned TileSizeSelectionModel::compute_tile_footprint_coefficient_for_access(
    PlutoAccess *acc, Band *band, PlutoProg *prog, double &coefficient) {
  unsigned num_tiled_levels = 0;
  unsigned depth = band->loop->depth + num_tiled_levels * band->width;

  PlutoMatrix *accMat = acc->mat;

  // Scalars
  if (accMat->nrows == 0)
    return 0;

  PlutoMatrix *newAccMat = get_linearly_independent_rows_in_mat(accMat);
  if (newAccMat->nrows == 0) {
    pluto_matrix_free(newAccMat);
    return 0;
  }

  unsigned nloops;
  Ploop **loops = pluto_get_loops_under(band->loop->stmts, band->loop->nstmts,
                                        depth, prog, &nloops);

  PlutoContext *context = prog->context;
  unsigned degree = 0;
  double gamma = 0.0f;
  for (unsigned i = 0; i < newAccMat->nrows; i++) {
    double gamma_i = 0.0f;
    bool has_vectorizable_loop = false;
    bool has_non_vector_dim_in_access = false;
    bool has_parallel_loop = false;
    unsigned num_pipe_parallel_access_dims = 0;
    for (unsigned j = 0; j < nloops; j++) {
      unsigned loop_depth = depth + loops[j]->depth;
      IF_DEBUG2(printf("loop_depth: %d, loop[%d]->depth: %d\n", loop_depth, j,
                       loops[j]->depth););
      if (newAccMat->val[i][loops[j]->depth] >= 1) {

        if ((loops[j]->depth) - depth == vectorDim) {
          has_vectorizable_loop = true;
          IF_DEBUG2(printf("Access has a vectorizable dimension at loop %d\n",
                           loops[j]->depth););
          continue;
        }
        if ((loops[j]->depth) - depth == par_loop_depth) {
          has_parallel_loop = true;
          IF_DEBUG2(printf("Access has a parallel dimension at loop %d\n",
                           loops[j]->depth););
          continue;
        }

        if (has_pipeline_parallelism &&
            (loops[j]->depth == depth || loops[j]->depth == depth + 1)) {
          IF_DEBUG2(printf("Access along pipe parallel dimension %d\n",
                           loops[j]->depth););
          num_pipe_parallel_access_dims++;
          continue;
        }

        has_non_vector_dim_in_access = true;
        IF_DEBUG2(
            printf("loop_depth : %d, loops[%d]->depth: %d, dim_reuse: %0.5f, "
                   "depth: %d, DimReuse size: %ld\n",
                   loop_depth, j, loops[j]->depth, dimReuse[loops[j]->depth],
                   depth, dimReuse.size()););
        gamma_i = gamma_i + dimReuse[loops[j]->depth - depth];
      }
    }
    IF_DEBUG(printf("gamma_i for access row  %d: %0.6f\n", i, gamma_i););
    if (gamma_i > 0.0f) {
      if (gamma == 0.0f)
        gamma = gamma_i;
      else
        gamma = gamma * gamma_i;
    }

    if (has_vectorizable_loop)
      gamma = (gamma == 0.0f) ? tile_size_for_vector_dimension
                              : gamma * tile_size_for_vector_dimension;

    // outer parallelism
    if (has_parallel_loop) {
      gamma = gamma * tile_size_for_parallel_dimension;
    }
    // Pipeline parallelism
    if (num_pipe_parallel_access_dims > 0)
      gamma = gamma * num_pipe_parallel_access_dims *
              tile_size_for_parallel_dimension;
    if (has_non_vector_dim_in_access)
      degree++;
  }
  pluto_matrix_free(newAccMat);
  pluto_loops_free(loops, nloops);
  coefficient = gamma;
  return degree;
}

/// Constructs a single variable n-degree polynomial for tile volume based on
/// dimensional reuse. The degree of the polynomial is bounded by the number of
/// loops in the permutable band.
void TileSizeSelectionModel::construct_expression_for_tile_volume(
    Band *band, PlutoProg *prog) {

  PlutoContext *context = prog->context;
  PlutoOptions *options = context->options;

  int num_tiled_levels = 0;
  int depth = band->loop->depth + num_tiled_levels * band->width;
  unsigned nstmts = band->loop->nstmts;
  Stmt **stmts = band->loop->stmts;
  unsigned num_loop_dimensions = 0;

  for (unsigned i = 0; i < band->width; i++) {
    unsigned j = 0;
    for (; j < nstmts; j++) {
      if (pluto_is_hyperplane_loop(stmts[j], depth + i))
        break;
    }
    if (j < nstmts)
      num_loop_dimensions++;
  }

  std::set<PlutoAccess *> distinctAccesses =
      get_distinct_accesses_in_band(band, prog);
  unsigned max_degree = 0;

  // Compute the tile footprint access-wise. Tile footprint is computed
  // according to dimensional reuse.
  for (auto acc : distinctAccesses) {
    double coefficient = 0.0f;
    if (pluto_matrix_get_rank(acc->mat) == 0) {
      pluto_access_free(acc);
      continue;
    }

    if (options->moredebug) {
      printf("Computing tile footprint coefficient for access\n");
      pluto_matrix_print(stdout, acc->mat);
    }

    unsigned degree = compute_tile_footprint_coefficient_for_access(
        acc, band, prog, coefficient);

    IF_DEBUG2(printf("Coefficient: %0.5f, degree :%d\n", coefficient, degree););

    if (degree == 0) {
      if (reuse_polynomial.size() == 0)
        reuse_polynomial.resize(1);
      reuse_polynomial[0] += coefficient;
      pluto_access_free(acc);
      continue;
    }

    if (degree > max_degree) {
      max_degree = degree;
      reuse_polynomial.resize(max_degree + 1);
      reuse_polynomial[max_degree] = 0.0f;
    }

    reuse_polynomial[degree] += coefficient;
    pluto_access_free(acc);
  }

  // f(\tau) - Cache_size = 0;
  reuse_polynomial[0] -= num_data_elements_in_cache;

  if (options->debug || options->moredebug) {
    printf("[tile-size-selection] Reuse Polynomial \n");
    for (unsigned i = 0; i < reuse_polynomial.size(); i++) {
      printf("\t Coefficient of degree %d: %0.5f\n", i, reuse_polynomial[i]);
    }
  }
}

/// For the input band, the routine returns a boolean vector whose values
/// corresponding to the dimensions that will be moved to the innermost level
/// are set to true.
static std::vector<bool> get_vectorizable_dimensions(Band *band,
                                                     PlutoProg *prog) {
  PlutoContext *context = prog->context;
  PlutoOptions *options = context->options;
  std::vector<bool> vectorizableDims(band->width, false);

  // TODO:Handle cases when intra tile optimizations are disabled
  if (!options->intratileopt)
    return vectorizableDims;

  // TODO: Refactor this code with pluto_intra_tile_optimize_band
  int num_tiled_levels = 0;
  unsigned nstmt_bands;
  Band **per_stmt_bands = get_per_stmt_band(band, &nstmt_bands);

  unsigned num_fused_bands;
  Band **ibands =
      fuse_per_stmt_bands(per_stmt_bands, nstmt_bands, num_tiled_levels,
                          &num_fused_bands, prog->context);

  if (num_fused_bands != nstmt_bands) {
    pluto_bands_free(per_stmt_bands, nstmt_bands);
  }
  // TODO:Remove the following debug code / update debug message.
  if (options->debug) {
    printf("Bands for intra tile optimization \n");
    pluto_bands_print(ibands, num_fused_bands);
  }
  for (unsigned i = 0; i < num_fused_bands; i++) {
    Band *band = ibands[i];
    int depth = band->loop->depth + num_tiled_levels * band->width;
    IF_DEBUG(printf("Getting loop at depth  %d\n", depth););

    unsigned nloops;
    Ploop **loops = pluto_get_loops_under(band->loop->stmts, band->loop->nstmts,
                                          depth, prog, &nloops);
    Ploop *best_loop = get_best_vectorizable_loop(loops, nloops, prog);
    if (options->debug) {
      printf("loop at depth %d considered for innermost tile_size\n",
             best_loop->depth);
      pluto_loop_print(best_loop);
    }

    // Do not select outermost parallel loop
    if (best_loop->depth - band->loop->depth == 0)
      continue;
    // Do not assign vector tile sizes to concurrent start dimensions.
    if (prog->is_diamond_tiled && best_loop->depth - band->loop->depth <= 1)
      continue;
    vectorizableDims[best_loop->depth - depth] = true;
    pluto_loops_free(loops, nloops);
  }
  pluto_bands_free(ibands, num_fused_bands);
  return vectorizableDims;
}

/// Returns true if dimensional reuse for all dimensions is zero.
static bool has_zero_dim_reuse(const std::vector<float> &dimReuse) {
  float max = *std::max_element(dimReuse.begin(), dimReuse.end());
  if (max == 0.0f)
    return true;
  return false;
}

/// Returns dimensional reuse for each dimension as a vector of floats for the
/// input band.
static std::vector<float> get_dimensional_reuse(Band *band, PlutoProg *prog) {
  std::vector<float> dimReuse(band->width, 0.0f);
  unsigned depth = band->loop->depth;
  unsigned nloops;
  PlutoContext *context = prog->context;
  Ploop **loops = pluto_get_loops_under(band->loop->stmts, band->loop->nstmts,
                                        depth, prog, &nloops);

  for (unsigned i = 0; i < nloops; i++) {
    dimReuse[loops[i]->depth - depth] = get_dimensional_reuse(loops[i]);
    // TODO: Use the following if we want to only consider distinct accesses.
    // To be experimented if this gives better results.
    // dimReuse[loops[i]->depth - depth] = get_num_invariant_accesses_in_stmts(
    //     loops[i]->stmts, loops[i]->nstmts, loops[i]->depth, prog);
    IF_DEBUG(printf("Dimensional reuse of dimension %d: %0.6f\n",
                    loops[i]->depth - depth,
                    dimReuse[loops[i]->depth - depth]););
  }
  pluto_loops_free(loops, nloops);

  return dimReuse;
}

/// Assign tile sizes for bands that have zero dimensional reuse for all
/// dimensions.
static std::vector<unsigned> get_tile_sizes_for_band_with_zero_dim_reuse(
    Band *band, const std::vector<float> &dimReuse,
    const std::vector<bool> &vectorizableDims, unsigned par_loop_depth) {
  std::vector<unsigned> tile_sizes(dimReuse.size());
  for (unsigned i = 0; i < dimReuse.size(); i++) {
    if (vectorizableDims[i]) {
      tile_sizes[i] = 512;
      continue;
    }
    if (i == par_loop_depth) {
      tile_sizes[i] = 32;
      continue;
    }
    tile_sizes[i] = 32;
  }
  return tile_sizes;
}

/// Returns true if the input band has wavefront parallelism, return false
/// otherwise.
// TODO: Check if this is sufficient for wavefront parallelism.
static bool pluto_band_has_wavefront_parallelism(PlutoProg *prog, Band *band) {
  if (pluto_loop_is_parallel(prog, band->loop))
    return false;

  /* Number of dimensions which are loops for all statements in this band. */
  int nloops = 1;
  for (unsigned depth = band->loop->depth + 1;
       depth < band->loop->depth + band->width; depth++) {
    unsigned j;
    for (j = 0; j < band->loop->nstmts; j++) {
      if (pluto_is_hyperplane_scalar(band->loop->stmts[j], depth))
        break;
    }
    if (j == band->loop->nstmts) {
      /* All of them are loops */
      nloops++;
    }
  }

  if (nloops <= 1) {
    /* Band doesn't have at least two dimensions for which all
     * statements have loops at those dimensions */
    return false;
  }
  return true;
}

/// Returns a vector of tile sizes. The tile sizes are proportional to the
/// dimensional reuse of each loop in the input band.
static std::vector<unsigned> find_tile_sizes(Band *band, PlutoProg *prog) {
  // Trivial band.
  if (band->width == 1)
    return std::vector<unsigned>(1, 32);

  PlutoContext *context = prog->context;
  PlutoOptions *options = prog->context->options;

  std::vector<bool> vectorizableDims = get_vectorizable_dimensions(band, prog);

  unsigned nloops;
  Ploop **loops = pluto_get_loops_under(band->loop->stmts, band->loop->nstmts,
                                        band->loop->depth, prog, &nloops);
  unsigned int num_tiled_levels = 0;
  unsigned depth = band->loop->depth + num_tiled_levels * band->width;
  unsigned par_loop_depth = vectorizableDims.size();
  for (unsigned i = 0; i < nloops; i++) {
    bool parallel = true;
    for (unsigned j = 0; j < loops[i]->nstmts; j++) {
      if (!pluto_loop_is_parallel_for_stmt(prog, loops[i],
                                           loops[i]->stmts[j])) {
        parallel = false;
        break;
      }
    }

    if (parallel) {
      par_loop_depth = loops[i]->depth - depth;
      IF_DEBUG(printf("Parallel loop\n"););
      IF_DEBUG(pluto_loop_print(loops[i]););
    }
  }
  pluto_loops_free(loops, nloops);

  IF_DEBUG(printf("Parallel dimension %d\n", par_loop_depth););
  if (vectorizableDims[par_loop_depth])
    vectorizableDims[par_loop_depth] = false;

  bool is_band_wavefront_parallel =
      pluto_band_has_wavefront_parallelism(prog, band);
  if (is_band_wavefront_parallel) {
    vectorizableDims[0] = false;
    vectorizableDims[1] = false;
    IF_DEBUG(printf("Band has wavefront parallelism\n"););
  }
  if (is_band_wavefront_parallel && band->width == 2) {
    std::vector<unsigned> tile_sizes(2, 32);
    return tile_sizes;
  }

  for (unsigned i = 0; i < vectorizableDims.size(); i++) {
    if (vectorizableDims[i])
      IF_DEBUG(
          printf(
              "Dimension %d is vectorizable for some statement in the band\n",
              i););
  }

  std::vector<float> dimReuse = get_dimensional_reuse(band, prog);
  // Early bailout when dimReuse of all dimensions is zero.
  if (has_zero_dim_reuse(dimReuse)) {
    auto tile_sizes = get_tile_sizes_for_band_with_zero_dim_reuse(
        band, dimReuse, vectorizableDims, par_loop_depth);
    return tile_sizes;
  }

  // Set the dimensional reuse of diamond tiling hyperplanes to be the same.
  if (prog->is_diamond_tiled)
    dimReuse[0] = dimReuse[1];

  IF_DEBUG(printf("After normalization\n"););
  for (unsigned i = 0; i < dimReuse.size(); i++) {
    IF_DEBUG(
        printf("Dimensional reuse of dimension %d: %0.5f \n", i, dimReuse[i]););
  }
  TileSizeSelectionModel *tss = new TileSizeSelectionModel(
      dimReuse, vectorizableDims, options->cache_size,
      options->data_element_size, options->ufactor, par_loop_depth,
      is_band_wavefront_parallel);
  tss->construct_expression_for_tile_volume(band, prog);

  unsigned root = tss->solve_reuse_expression();
  IF_DEBUG(printf("Solution of the reuse polynomial: %d\n", root););
  auto tile_sizes = tss->infer_tile_sizes(root);

  delete tss;
  IF_DEBUG(printf("Tile sizes\n"););
  for (auto i : tile_sizes)
    IF_DEBUG(printf("%d\t", i););
  IF_DEBUG(printf("\n"););
  return tile_sizes;
}

void find_tile_sizes(Band *band, PlutoProg *prog, int *tile_sizes) {
  unsigned firstLoop = band->loop->depth;
  unsigned num_tile_dims = band->width;
  unsigned nstmts = band->loop->nstmts;
  Stmt **stmts = band->loop->stmts;

  double t_start = rtclock();
  std::vector<unsigned> tileSizeVec = find_tile_sizes(band, prog);
  double t_end = rtclock();
  prog->tss_time += t_end - t_start;
  assert(tileSizeVec.size() == num_tile_dims);

  for (unsigned i = 0; i < num_tile_dims; i++) {
    unsigned j = 0;
    for (; j < nstmts; j++) {
      if (pluto_is_hyperplane_loop(stmts[j], firstLoop + i))
        break;
    }
    if (j < nstmts)
      tile_sizes[i] = tileSizeVec[i];
    else
      tile_sizes[i] = 42;
  }

  PlutoOptions *options = prog->context->options;
  if (options->debug) {
    printf("Tile sizes for band\n");
    pluto_band_print(band);
    for (unsigned i = 0; i < num_tile_dims; i++) {
      printf("TileSize for dimension %d: %d\n", i, tile_sizes[i]);
    }
  }

  return;
}
