#include <R.h>
#include <Rinternals.h>
#include "bin_table.h"

/* Check if C functions are available */
SEXP has_resolution_depth_c()
{
    return ScalarLogical(TRUE);
}

/* Main C function for calculating resolution depth */
SEXP calculate_resolution_depth_c(
    SEXP read_names,
    SEXP chrom1_vec,
    SEXP pos1_vec,
    SEXP chrom2_vec,
    SEXP pos2_vec,
    SEXP bin_size_r)
{

    int n = LENGTH(read_names);
    int bin_size = asInteger(bin_size_r);

    /* Create hash table for bins */
    BinTable *table = create_bin_table(1024);

    /* Process each pair */
    for (int i = 0; i < n; i++)
    {
        int pos1 = INTEGER(pos1_vec)[i];
        int pos2 = INTEGER(pos2_vec)[i];
        const char *chrom1 = CHAR(STRING_ELT(chrom1_vec, i));
        const char *chrom2 = CHAR(STRING_ELT(chrom2_vec, i));

        /* Calculate bins (ceiling division) */
        int bin1 = (pos1 + bin_size - 1) / bin_size;
        int bin2 = (pos2 + bin_size - 1) / bin_size;

        /* Add to table */
        bin_table_add(table, chrom1, bin1);
        bin_table_add(table, chrom2, bin2);
    }

    /* Get all entries as array for R export */
    BinEntry *entries_array;
    int entries_count;
    bin_table_get_all_entries(table, &entries_array, &entries_count);

    /* Create output data frame */
    SEXP chrom_out, bin_out, count_out, df, names;

    PROTECT(chrom_out = allocVector(STRSXP, entries_count));
    PROTECT(bin_out = allocVector(INTSXP, entries_count));
    PROTECT(count_out = allocVector(INTSXP, entries_count));

    for (int i = 0; i < entries_count; i++)
    {
        SET_STRING_ELT(chrom_out, i, mkChar(entries_array[i].chrom));
        INTEGER(bin_out)
        [i] = entries_array[i].bin;
        INTEGER(count_out)
        [i] = entries_array[i].count;
    }

    /* Create data frame */
    PROTECT(df = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(df, 0, chrom_out);
    SET_VECTOR_ELT(df, 1, bin_out);
    SET_VECTOR_ELT(df, 2, count_out);

    /* Set names */
    PROTECT(names = allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, mkChar("chrom"));
    SET_STRING_ELT(names, 1, mkChar("bin"));
    SET_STRING_ELT(names, 2, mkChar("count"));
    setAttrib(df, R_NamesSymbol, names);

    /* Set row names */
    SEXP row_names;
    PROTECT(row_names = allocVector(INTSXP, 2));
    INTEGER(row_names)
    [0] = NA_INTEGER;
    INTEGER(row_names)
    [1] = -entries_count;
    setAttrib(df, R_RowNamesSymbol, row_names);

    /* Set class */
    SEXP class;
    PROTECT(class = allocVector(STRSXP, 1));
    SET_STRING_ELT(class, 0, mkChar("data.frame"));
    setAttrib(df, R_ClassSymbol, class);

    R_Free(entries_array); /* Free the temporary array */
    free_bin_table(table);

    UNPROTECT(7);
    return df;
}
