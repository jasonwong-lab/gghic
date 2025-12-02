#include <R.h>
#include <Rinternals.h>
#include <R_ext/Memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "bin_table.h"

#define BUFFER_SIZE 65536
#define INITIAL_CAPACITY 100000

/* Line buffer for file reading */
typedef struct
{
    char *buffer;
    size_t size;
    size_t capacity;
} LineBuffer;

LineBuffer *create_line_buffer(size_t initial_capacity)
{
    LineBuffer *buf = (LineBuffer *)R_alloc(1, sizeof(LineBuffer));
    buf->buffer = (char *)R_alloc(initial_capacity, sizeof(char));
    buf->size = 0;
    buf->capacity = initial_capacity;
    return buf;
}

void free_line_buffer(LineBuffer *buf)
{
    /* R_alloc memory is automatically freed, so this is now a no-op */
    /* Keep function for API compatibility */
}

void expand_line_buffer(LineBuffer *buf, size_t needed)
{
    while (buf->capacity < needed)
    {
        buf->capacity *= 2;
    }
    char *new_buffer = (char *)R_alloc(buf->capacity, sizeof(char));
    memcpy(new_buffer, buf->buffer, buf->size);
    buf->buffer = new_buffer;
}

/* Read line from file (handles both gzip and plain text) */
int read_line(gzFile gz, FILE *f, LineBuffer *buf)
{
    buf->size = 0;

    if (gz != NULL)
    {
        /* Reading from gzip file */
        char c;
        while ((c = gzgetc(gz)) != -1)
        {
            if (buf->size + 1 >= buf->capacity)
            {
                expand_line_buffer(buf, buf->capacity * 2);
            }

            buf->buffer[buf->size++] = c;

            if (c == '\n')
            {
                buf->buffer[buf->size - 1] = '\0'; /* Replace \n with \0 */
                return 1;
            }
        }
    }
    else if (f != NULL)
    {
        /* Reading from plain text file */
        char c;
        while ((c = fgetc(f)) != EOF)
        {
            if (buf->size + 1 >= buf->capacity)
            {
                expand_line_buffer(buf, buf->capacity * 2);
            }

            buf->buffer[buf->size++] = c;

            if (c == '\n')
            {
                buf->buffer[buf->size - 1] = '\0';
                return 1;
            }
        }
    }

    return 0; /* EOF */
}

/* Parse a single line of pairs data */
int parse_pairs_line(const char *line, int bin_size, int *pos1, int *pos2,
                     char *chrom1, char *chrom2, int max_chrom_len)
{
    char read_name[256];

    /* Parse first 5 required fields: read_name chrom1 pos1 chrom2 pos2
       Any additional fields after pos2 are ignored */
    int parsed = sscanf(line, "%s %s %d %s %d",
                        read_name, chrom1, pos1, chrom2, pos2);

    if (parsed >= 5)
    {
        return 1; /* Successfully parsed required fields */
    }

    return 0; /* Invalid line - not enough fields */
}

static int is_gzip(const char *filename)
{
    size_t len = strlen(filename);
    return (len > 3 && strcmp(filename + len - 3, ".gz") == 0);
}

/* ===== PairsData structure for hypergraph functions ===== */

typedef struct
{
    char **read_names;
    char **chrom1;
    int *pos1;
    char **chrom2;
    int *pos2;
    size_t size;
    size_t capacity;
} PairsData;

static void init_pairs_data(PairsData *data)
{
    data->capacity = INITIAL_CAPACITY;
    data->size = 0;
    data->read_names = (char **)R_Calloc(data->capacity, char *);
    data->chrom1 = (char **)R_Calloc(data->capacity, char *);
    data->pos1 = (int *)R_Calloc(data->capacity, int);
    data->chrom2 = (char **)R_Calloc(data->capacity, char *);
    data->pos2 = (int *)R_Calloc(data->capacity, int);

    if (!data->read_names || !data->chrom1 || !data->pos1 ||
        !data->chrom2 || !data->pos2)
    {
        error("Memory allocation failed");
    }
}

static void expand_pairs_data(PairsData *data)
{
    data->capacity *= 2;
    data->read_names = (char **)R_Realloc(data->read_names, data->capacity, char *);
    data->chrom1 = (char **)R_Realloc(data->chrom1, data->capacity, char *);
    data->pos1 = (int *)R_Realloc(data->pos1, data->capacity, int);
    data->chrom2 = (char **)R_Realloc(data->chrom2, data->capacity, char *);
    data->pos2 = (int *)R_Realloc(data->pos2, data->capacity, int);

    if (!data->read_names || !data->chrom1 || !data->pos1 ||
        !data->chrom2 || !data->pos2)
    {
        error("Memory reallocation failed");
    }
}

static void add_pair(PairsData *data, const char *read_name,
                     const char *chr1, int p1, const char *chr2, int p2)
{
    if (data->size >= data->capacity)
    {
        expand_pairs_data(data);
    }

    /* Allocate and copy strings using R_alloc */
    size_t len_read = strlen(read_name) + 1;
    data->read_names[data->size] = (char *)R_alloc(len_read, sizeof(char));
    strcpy(data->read_names[data->size], read_name);

    size_t len_chr1 = strlen(chr1) + 1;
    data->chrom1[data->size] = (char *)R_alloc(len_chr1, sizeof(char));
    strcpy(data->chrom1[data->size], chr1);

    data->pos1[data->size] = p1;

    size_t len_chr2 = strlen(chr2) + 1;
    data->chrom2[data->size] = (char *)R_alloc(len_chr2, sizeof(char));
    strcpy(data->chrom2[data->size], chr2);

    data->pos2[data->size] = p2;
    data->size++;
}

static void free_pairs_data(PairsData *data)
{
    /* Strings allocated with R_alloc are auto-freed */
    /* Free the arrays allocated with Calloc/Realloc */
    R_Free(data->read_names);
    R_Free(data->chrom1);
    R_Free(data->pos1);
    R_Free(data->chrom2);
    R_Free(data->pos2);
}

/* Filter modes for reading pairs */
typedef enum
{
    FILTER_NONE = 0,            /* Keep all contacts */
    FILTER_INTRA_ONLY = 1,      /* Keep only intra-chromosomal */
    FILTER_CHROM = 2,           /* Keep all pairs involving specified chromosome(s) */
    FILTER_CHROM_INTRA_ONLY = 3 /* Keep only intra-chromosomal pairs for specified chromosome(s) */
} FilterMode;

/* Check if a chromosome is in the target chromosome array */
static int chrom_in_list(const char *chrom, const char **chrom_list, int n_chroms)
{
    if (chrom_list == NULL || n_chroms == 0)
        return 0;

    for (int i = 0; i < n_chroms; i++)
    {
        if (strcmp(chrom, chrom_list[i]) == 0)
            return 1;
    }
    return 0;
}

/* Check if a pair passes the filter */
static int passes_filter(const char *chrom1, const char *chrom2,
                         FilterMode mode, const char **target_chroms, int n_chroms)
{
    switch (mode)
    {
    case FILTER_NONE:
        return 1;
    case FILTER_INTRA_ONLY:
        return strcmp(chrom1, chrom2) == 0;
    case FILTER_CHROM:
        /* For single chromosome: keep if either is in list (for inter-chrom with that chrom) */
        /* For multiple chromosomes: keep if BOTH are in list (for inter/intra between those chroms) */
        if (n_chroms == 1)
        {
            return chrom_in_list(chrom1, target_chroms, n_chroms) ||
                   chrom_in_list(chrom2, target_chroms, n_chroms);
        }
        else
        {
            return chrom_in_list(chrom1, target_chroms, n_chroms) &&
                   chrom_in_list(chrom2, target_chroms, n_chroms);
        }
    case FILTER_CHROM_INTRA_ONLY:
        /* Keep only if both chroms are in target list AND it's intra-chromosomal */
        return strcmp(chrom1, chrom2) == 0 &&
               chrom_in_list(chrom1, target_chroms, n_chroms);
    default:
        return 0;
    }
}

/* Unified function to read pairs with different filters */
static void read_pairs_filtered(const char *filename, PairsData *data,
                                FilterMode mode, const char **target_chroms, int n_chroms)
{
    char buffer[BUFFER_SIZE];
    char read_name[256], strand1[2], chrom1[64], strand2[2], chrom2[64];
    int pos1, frag1, pos2, frag2, mapq1, mapq2;
    size_t line_count = 0;

    if (is_gzip(filename))
    {
        gzFile fp = gzopen(filename, "rb");
        if (!fp)
        {
            error("Cannot open gzipped file: %s", filename);
        }

        while (gzgets(fp, buffer, BUFFER_SIZE) != NULL)
        {
            line_count++;

            if (buffer[0] == '\n' || buffer[0] == '\0')
                continue;

            int parsed = sscanf(buffer,
                                "%255s\t%1s\t%63s\t%d\t%d\t%1s\t%63s\t%d\t%d\t%d\t%d",
                                read_name, strand1, chrom1, &pos1, &frag1,
                                strand2, chrom2, &pos2, &frag2, &mapq1, &mapq2);

            if (parsed == 11)
            {
                if (passes_filter(chrom1, chrom2, mode, target_chroms, n_chroms))
                {
                    add_pair(data, read_name, chrom1, pos1, chrom2, pos2);
                }
            }
            else
            {
                parsed = sscanf(buffer, "%255s\t%63s\t%d\t%63s\t%d",
                                read_name, chrom1, &pos1, chrom2, &pos2);

                if (parsed >= 5 && passes_filter(chrom1, chrom2, mode, target_chroms, n_chroms))
                {
                    add_pair(data, read_name, chrom1, pos1, chrom2, pos2);
                }
            }

            if (line_count % 1000000 == 0)
            {
                Rprintf("\rProcessed %zu million lines, kept %zu contacts",
                        line_count / 1000000, data->size);
            }
        }

        gzclose(fp);
    }
    else
    {
        FILE *fp = fopen(filename, "r");
        if (!fp)
        {
            error("Cannot open file: %s", filename);
        }

        while (fgets(buffer, BUFFER_SIZE, fp) != NULL)
        {
            line_count++;

            if (buffer[0] == '\n' || buffer[0] == '\0')
                continue;

            int parsed = sscanf(buffer,
                                "%255s\t%1s\t%63s\t%d\t%d\t%1s\t%63s\t%d\t%d\t%d\t%d",
                                read_name, strand1, chrom1, &pos1, &frag1,
                                strand2, chrom2, &pos2, &frag2, &mapq1, &mapq2);

            if (parsed == 11)
            {
                if (passes_filter(chrom1, chrom2, mode, target_chroms, n_chroms))
                {
                    add_pair(data, read_name, chrom1, pos1, chrom2, pos2);
                }
            }
            else
            {
                parsed = sscanf(buffer, "%255s\t%63s\t%d\t%63s\t%d",
                                read_name, chrom1, &pos1, chrom2, &pos2);

                if (parsed >= 5 && passes_filter(chrom1, chrom2, mode, target_chroms, n_chroms))
                {
                    add_pair(data, read_name, chrom1, pos1, chrom2, pos2);
                }
            }

            if (line_count % 1000000 == 0)
            {
                Rprintf("Processed %zu million lines, kept %zu contacts\n",
                        line_count / 1000000, data->size);
            }
        }

        fclose(fp);
    }

    Rprintf("\nTotal: %zu lines processed, %zu contacts retained\n",
            line_count, data->size);
}

/* Convert PairsData to R data frame */
static SEXP pairs_data_to_df(PairsData *data)
{
    SEXP result = PROTECT(allocVector(VECSXP, 5));
    SEXP read_names_r = PROTECT(allocVector(STRSXP, data->size));
    SEXP chrom1_r = PROTECT(allocVector(STRSXP, data->size));
    SEXP pos1_r = PROTECT(allocVector(INTSXP, data->size));
    SEXP chrom2_r = PROTECT(allocVector(STRSXP, data->size));
    SEXP pos2_r = PROTECT(allocVector(INTSXP, data->size));

    for (size_t i = 0; i < data->size; i++)
    {
        SET_STRING_ELT(read_names_r, i, mkChar(data->read_names[i]));
        SET_STRING_ELT(chrom1_r, i, mkChar(data->chrom1[i]));
        INTEGER(pos1_r)
        [i] = data->pos1[i];
        SET_STRING_ELT(chrom2_r, i, mkChar(data->chrom2[i]));
        INTEGER(pos2_r)
        [i] = data->pos2[i];
    }

    SET_VECTOR_ELT(result, 0, read_names_r);
    SET_VECTOR_ELT(result, 1, chrom1_r);
    SET_VECTOR_ELT(result, 2, pos1_r);
    SET_VECTOR_ELT(result, 3, chrom2_r);
    SET_VECTOR_ELT(result, 4, pos2_r);

    SEXP colnames = PROTECT(allocVector(STRSXP, 5));
    SET_STRING_ELT(colnames, 0, mkChar("read_name"));
    SET_STRING_ELT(colnames, 1, mkChar("chrom1"));
    SET_STRING_ELT(colnames, 2, mkChar("pos1"));
    SET_STRING_ELT(colnames, 3, mkChar("chrom2"));
    SET_STRING_ELT(colnames, 4, mkChar("pos2"));
    setAttrib(result, R_NamesSymbol, colnames);

    SEXP rownames = PROTECT(allocVector(INTSXP, 2));
    INTEGER(rownames)
    [0] = NA_INTEGER;
    INTEGER(rownames)
    [1] = -(int)data->size;
    setAttrib(result, R_RowNamesSymbol, rownames);

    SEXP df_class = PROTECT(allocVector(STRSXP, 1));
    SET_STRING_ELT(df_class, 0, mkChar("data.frame"));
    setAttrib(result, R_ClassSymbol, df_class);

    UNPROTECT(9);
    return result;
}

/* ===== Exported R functions ===== */

/* Read pairs filtered by chromosome */
SEXP read_pairs_chrom_c(SEXP r_filename, SEXP r_chrom, SEXP r_inter_chrom)
{
    const char *filename = CHAR(STRING_ELT(r_filename, 0));
    int inter_chrom = asLogical(r_inter_chrom);

    /* Handle single or multiple chromosomes */
    int n_chroms = length(r_chrom);
    const char **target_chroms = (const char **)R_alloc(n_chroms, sizeof(char *));

    for (int i = 0; i < n_chroms; i++)
    {
        target_chroms[i] = CHAR(STRING_ELT(r_chrom, i));
    }

    PairsData data;
    init_pairs_data(&data);

    Rprintf("Reading pairs from: %s\n", filename);
    if (n_chroms == 1)
    {
        Rprintf("Filtering for chromosome: %s\n", target_chroms[0]);
    }
    else
    {
        Rprintf("Filtering for %d chromosomes\n", n_chroms);
    }

    FilterMode mode = inter_chrom ? FILTER_CHROM : FILTER_CHROM_INTRA_ONLY;
    const char *mode_str = inter_chrom ? "(intra- and inter-chromosomal)" : "(intra-chromosomal only)";
    Rprintf("Mode: %s\n", mode_str);

    read_pairs_filtered(filename, &data, mode, target_chroms, n_chroms);

    if (data.size == 0)
    {
        /* R_alloc memory is auto-freed */
        free_pairs_data(&data);
        if (n_chroms == 1)
        {
            error("No contacts found for specified chromosome");
        }
        else
        {
            error("No contacts found for specified chromosomes");
        }
    }

    /* R_alloc memory is auto-freed */
    SEXP result = pairs_data_to_df(&data);
    free_pairs_data(&data);
    return result;
}

/* Read all intra-chromosomal contacts */
SEXP read_pairs_all_intra_c(SEXP r_filename)
{
    const char *filename = CHAR(STRING_ELT(r_filename, 0));

    PairsData data;
    init_pairs_data(&data);

    Rprintf("Reading all intra-chromosomal pairs from: %s\n", filename);

    read_pairs_filtered(filename, &data, FILTER_INTRA_ONLY, NULL, 0);

    if (data.size == 0)
    {
        free_pairs_data(&data);
        error("No intra-chromosomal contacts found");
    }

    SEXP result = pairs_data_to_df(&data);
    free_pairs_data(&data);
    return result;
}

/* Read all contacts (both intra- and inter-chromosomal) */
SEXP read_pairs_all_c(SEXP r_filename)
{
    const char *filename = CHAR(STRING_ELT(r_filename, 0));

    PairsData data;
    init_pairs_data(&data);

    Rprintf("Reading all pairs (intra- and inter-chromosomal) from: %s\n", filename);

    read_pairs_filtered(filename, &data, FILTER_NONE, NULL, 0);

    if (data.size == 0)
    {
        free_pairs_data(&data);
        error("No contacts found");
    }

    SEXP result = pairs_data_to_df(&data);
    free_pairs_data(&data);
    return result;
}

/* Process pairs file in streaming fashion */
SEXP process_pairs_file_c(SEXP pairs_file_r, SEXP bin_size_r)
{
    const char *pairs_file = CHAR(STRING_ELT(pairs_file_r, 0));
    int bin_size = asInteger(bin_size_r);

    gzFile gz = NULL;
    FILE *f = NULL;

    if (is_gzip(pairs_file))
    {
        gz = gzopen(pairs_file, "rb");
        if (!gz)
        {
            error("Cannot open gzip file: %s", pairs_file);
        }
    }
    else
    {
        f = fopen(pairs_file, "r");
        if (!f)
        {
            error("Cannot open file: %s", pairs_file);
        }
    }

    BinTable *table = create_bin_table(1024);
    LineBuffer *buf = create_line_buffer(4096);

    int line_count = 0;
    int chunk_count = 0;

    while (read_line(gz, f, buf))
    {
        if (buf->size == 0)
            continue;

        char chrom1[256], chrom2[256];
        int pos1, pos2;

        if (parse_pairs_line(buf->buffer, bin_size, &pos1, &pos2,
                             chrom1, chrom2, 255))
        {
            int bin1 = (pos1 + bin_size - 1) / bin_size;
            int bin2 = (pos2 + bin_size - 1) / bin_size;

            bin_table_add(table, chrom1, bin1);
            bin_table_add(table, chrom2, bin2);

            line_count++;
            if (line_count % 100000 == 0)
            {
                chunk_count++;
                Rprintf("\rProcessed %d lines (%d chunks)...", line_count, chunk_count);
                R_FlushConsole();
            }
        }

        R_CheckUserInterrupt();
    }

    if (gz)
        gzclose(gz);
    if (f)
        fclose(f);

    BinEntry *entries_array;
    int entries_count;
    bin_table_get_all_entries(table, &entries_array, &entries_count);

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

    PROTECT(df = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(df, 0, chrom_out);
    SET_VECTOR_ELT(df, 1, bin_out);
    SET_VECTOR_ELT(df, 2, count_out);

    PROTECT(names = allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, mkChar("chrom"));
    SET_STRING_ELT(names, 1, mkChar("bin"));
    SET_STRING_ELT(names, 2, mkChar("count"));
    setAttrib(df, R_NamesSymbol, names);

    SEXP row_names;
    PROTECT(row_names = allocVector(INTSXP, 2));
    INTEGER(row_names)
    [0] = NA_INTEGER;
    INTEGER(row_names)
    [1] = -entries_count;
    setAttrib(df, R_RowNamesSymbol, row_names);

    SEXP class;
    PROTECT(class = allocVector(STRSXP, 1));
    SET_STRING_ELT(class, 0, mkChar("data.frame"));
    setAttrib(df, R_ClassSymbol, class);

    free_line_buffer(buf);
    R_Free(entries_array);
    free_bin_table(table);

    UNPROTECT(7);
    return df;
}

/* Calculate coverage from cached positions (vector input) */
SEXP calculate_coverage_from_positions_c(
    SEXP chrom1_vec, SEXP pos1_vec, SEXP chrom2_vec, SEXP pos2_vec,
    SEXP bin_size_r, SEXP min_contacts_r)
{
    int n = LENGTH(chrom1_vec);
    int bin_size = asInteger(bin_size_r);
    int min_contacts = asInteger(min_contacts_r);

    BinTable *table = create_bin_table(1024);

    for (int i = 0; i < n; i++)
    {
        int pos1 = INTEGER(pos1_vec)[i];
        int pos2 = INTEGER(pos2_vec)[i];
        const char *chrom1 = CHAR(STRING_ELT(chrom1_vec, i));
        const char *chrom2 = CHAR(STRING_ELT(chrom2_vec, i));

        int bin1 = (pos1 + bin_size - 1) / bin_size;
        int bin2 = (pos2 + bin_size - 1) / bin_size;

        bin_table_add(table, chrom1, bin1);
        bin_table_add(table, chrom2, bin2);
    }

    int total_bins = table->size;
    int bins_above_threshold = 0;

    for (int i = 0; i < table->capacity; i++)
    {
        BinEntry *entry = table->buckets[i];
        while (entry)
        {
            if (entry->count >= min_contacts)
            {
                bins_above_threshold++;
            }
            entry = entry->next;
        }
    }

    double coverage = (total_bins > 0) ? ((double)bins_above_threshold / total_bins) : 0.0;

    free_bin_table(table);

    return ScalarReal(coverage);
}

/* PairsCache structure for external pointer */
typedef struct
{
    char **chrom1;
    int *pos1;
    char **chrom2;
    int *pos2;
    int size;
} PairsCache;

static void finalize_pairs_cache(SEXP ext_ptr)
{
    if (TYPEOF(ext_ptr) != EXTPTRSXP)
        return;

    PairsCache *cache = (PairsCache *)R_ExternalPtrAddr(ext_ptr);
    if (cache == NULL)
        return;

    for (int i = 0; i < cache->size; i++)
    {
        if (cache->chrom1[i])
            R_Free(cache->chrom1[i]);
        if (cache->chrom2[i])
            R_Free(cache->chrom2[i]);
    }
    R_Free(cache->chrom1);
    R_Free(cache->pos1);
    R_Free(cache->chrom2);
    R_Free(cache->pos2);
    R_Free(cache);

    R_ClearExternalPtr(ext_ptr);
}

/* Read pairs file and cache in C memory */
SEXP read_pairs_cache_c(SEXP pairs_file_r)
{
    const char *pairs_file = CHAR(STRING_ELT(pairs_file_r, 0));

    gzFile gz = NULL;
    FILE *f = NULL;

    if (is_gzip(pairs_file))
    {
        gz = gzopen(pairs_file, "rb");
        if (!gz)
        {
            error("Cannot open gzip file: %s", pairs_file);
        }
    }
    else
    {
        f = fopen(pairs_file, "r");
        if (!f)
        {
            error("Cannot open file: %s", pairs_file);
        }
    }

    PairsCache *cache = (PairsCache *)R_Calloc(1, PairsCache);
    int capacity = 100000;
    cache->size = 0;
    cache->chrom1 = (char **)R_Calloc(capacity, char *);
    cache->pos1 = (int *)R_Calloc(capacity, int);
    cache->chrom2 = (char **)R_Calloc(capacity, char *);
    cache->pos2 = (int *)R_Calloc(capacity, int);

    LineBuffer *buf = create_line_buffer(4096);
    int line_count = 0;

    while (read_line(gz, f, buf))
    {
        if (buf->size == 0)
            continue;

        char chrom1[256], chrom2[256];
        int pos1, pos2;

        if (parse_pairs_line(buf->buffer, 1, &pos1, &pos2, chrom1, chrom2, 255))
        {
            if (cache->size >= capacity)
            {
                capacity *= 2;
                cache->chrom1 = (char **)R_Realloc(cache->chrom1, capacity, char *);
                cache->pos1 = (int *)R_Realloc(cache->pos1, capacity, int);
                cache->chrom2 = (char **)R_Realloc(cache->chrom2, capacity, char *);
                cache->pos2 = (int *)R_Realloc(cache->pos2, capacity, int);
            }

            size_t len1 = strlen(chrom1) + 1;
            cache->chrom1[cache->size] = (char *)R_Calloc(len1, char);
            strcpy(cache->chrom1[cache->size], chrom1);
            size_t len2 = strlen(chrom2) + 1;
            cache->chrom2[cache->size] = (char *)R_Calloc(len2, char);
            strcpy(cache->chrom2[cache->size], chrom2);
            cache->pos1[cache->size] = pos1;
            cache->pos2[cache->size] = pos2;

            cache->size++;
            line_count++;

            if (line_count % 100000 == 0)
            {
                Rprintf("\rReading positions: %d lines...", line_count);
                R_FlushConsole();
            }
        }

        R_CheckUserInterrupt();
    }

    Rprintf("\nRead %d positions into C memory\n", line_count);

    if (gz)
        gzclose(gz);
    if (f)
        fclose(f);
    free_line_buffer(buf);

    SEXP ext_ptr = PROTECT(R_MakeExternalPtr(cache, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(ext_ptr, finalize_pairs_cache, TRUE);

    UNPROTECT(1);
    return ext_ptr;
}

/* Calculate coverage from cached C memory */
SEXP calculate_coverage_from_cache_c(
    SEXP cache_ptr, SEXP bin_size_r, SEXP min_contacts_r)
{
    if (TYPEOF(cache_ptr) != EXTPTRSXP)
    {
        error("cache_ptr must be an external pointer");
    }

    PairsCache *cache = (PairsCache *)R_ExternalPtrAddr(cache_ptr);
    if (cache == NULL)
    {
        error("Invalid external pointer (NULL)");
    }

    int bin_size = asInteger(bin_size_r);
    int min_contacts = asInteger(min_contacts_r);

    BinTable *table = create_bin_table(16384);

    for (int i = 0; i < cache->size; i++)
    {
        int bin1 = (cache->pos1[i] + bin_size - 1) / bin_size;
        int bin2 = (cache->pos2[i] + bin_size - 1) / bin_size;

        bin_table_add(table, cache->chrom1[i], bin1);
        bin_table_add(table, cache->chrom2[i], bin2);
    }

    int total_bins = table->size;
    int bins_above_threshold = 0;

    for (int i = 0; i < table->capacity; i++)
    {
        BinEntry *entry = table->buckets[i];
        while (entry)
        {
            if (entry->count >= min_contacts)
            {
                bins_above_threshold++;
            }
            entry = entry->next;
        }
    }

    double coverage = (total_bins > 0) ? ((double)bins_above_threshold / total_bins) : 0.0;

    free_bin_table(table);

    return ScalarReal(coverage);
}

/* Process pairs from cache and return depth table */
SEXP process_pairs_from_cache_c(SEXP cache_ptr, SEXP bin_size_r)
{
    if (TYPEOF(cache_ptr) != EXTPTRSXP)
    {
        error("cache_ptr must be an external pointer");
    }

    PairsCache *cache = (PairsCache *)R_ExternalPtrAddr(cache_ptr);
    if (cache == NULL)
    {
        error("Invalid external pointer (NULL)");
    }

    int bin_size = asInteger(bin_size_r);
    BinTable *table = create_bin_table(16384);

    /* Process cached positions */
    for (int i = 0; i < cache->size; i++)
    {
        int bin1 = (cache->pos1[i] + bin_size - 1) / bin_size;
        int bin2 = (cache->pos2[i] + bin_size - 1) / bin_size;

        bin_table_add(table, cache->chrom1[i], bin1);
        bin_table_add(table, cache->chrom2[i], bin2);
    }

    /* Convert to data frame */
    BinEntry *entries_array;
    int entries_count;
    bin_table_get_all_entries(table, &entries_array, &entries_count);

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

    PROTECT(df = allocVector(VECSXP, 3));
    SET_VECTOR_ELT(df, 0, chrom_out);
    SET_VECTOR_ELT(df, 1, bin_out);
    SET_VECTOR_ELT(df, 2, count_out);

    PROTECT(names = allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, mkChar("chrom"));
    SET_STRING_ELT(names, 1, mkChar("bin"));
    SET_STRING_ELT(names, 2, mkChar("count"));
    setAttrib(df, R_NamesSymbol, names);

    SEXP row_names;
    PROTECT(row_names = allocVector(INTSXP, 2));
    INTEGER(row_names)
    [0] = NA_INTEGER;
    INTEGER(row_names)
    [1] = -entries_count;
    setAttrib(df, R_RowNamesSymbol, row_names);

    SEXP class;
    PROTECT(class = allocVector(STRSXP, 1));
    SET_STRING_ELT(class, 0, mkChar("data.frame"));
    setAttrib(df, R_ClassSymbol, class);

    R_Free(entries_array);
    free_bin_table(table);

    UNPROTECT(7);
    return df;
}
