#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

/* Forward declarations */
SEXP has_resolution_depth_c();
SEXP calculate_resolution_depth_c(SEXP read_names, SEXP chrom1_vec,
                                  SEXP pos1_vec, SEXP chrom2_vec,
                                  SEXP pos2_vec, SEXP bin_size_r);
SEXP process_pairs_file_c(SEXP pairs_file_r, SEXP bin_size_r);
SEXP calculate_coverage_from_positions_c(SEXP chrom1_vec, SEXP pos1_vec,
                                         SEXP chrom2_vec, SEXP pos2_vec,
                                         SEXP bin_size_r, SEXP min_contacts_r);
SEXP read_pairs_cache_c(SEXP pairs_file_r);
SEXP calculate_coverage_from_cache_c(SEXP cache_ptr, SEXP bin_size_r,
                                     SEXP min_contacts_r);
SEXP process_pairs_from_cache_c(SEXP cache_ptr, SEXP bin_size_r);
SEXP read_pairs_chrom_c(SEXP r_filename, SEXP r_chrom, SEXP r_inter_chrom);
SEXP read_pairs_all_intra_c(SEXP r_filename);
SEXP read_pairs_all_c(SEXP r_filename);

/* Method registration table */
static const R_CallMethodDef callMethods[] = {
    {"has_resolution_depth_c", (DL_FUNC)&has_resolution_depth_c, 0},
    {"calculate_resolution_depth_c", (DL_FUNC)&calculate_resolution_depth_c, 6},
    {"process_pairs_file_c", (DL_FUNC)&process_pairs_file_c, 2},
    {"calculate_coverage_from_positions_c", (DL_FUNC)&calculate_coverage_from_positions_c, 6},
    {"read_pairs_cache_c", (DL_FUNC)&read_pairs_cache_c, 1},
    {"calculate_coverage_from_cache_c", (DL_FUNC)&calculate_coverage_from_cache_c, 3},
    {"process_pairs_from_cache_c", (DL_FUNC)&process_pairs_from_cache_c, 2},
    {"read_pairs_chrom_c", (DL_FUNC)&read_pairs_chrom_c, 3},
    {"read_pairs_all_intra_c", (DL_FUNC)&read_pairs_all_intra_c, 1},
    {"read_pairs_all_c", (DL_FUNC)&read_pairs_all_c, 1},
    {NULL, NULL, 0}};

/* Registration routine */
void R_init_gghic(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
