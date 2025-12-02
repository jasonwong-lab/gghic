#ifndef BIN_TABLE_H
#define BIN_TABLE_H

#include <R.h>
#include <R_ext/Memory.h>
#include <stdlib.h>
#include <string.h>

/* Hash table for unique bins using proper hashing */
typedef struct BinEntry
{
    char *chrom;
    int bin;
    int count;
    struct BinEntry *next; /* For chaining collisions */
} BinEntry;

typedef struct
{
    BinEntry **buckets; /* Array of bucket pointers */
    int capacity;       /* Number of buckets */
    int size;           /* Total entries */
} BinTable;

/* Forward declarations */
static inline BinTable *create_bin_table(int initial_capacity);
static inline void free_bin_table(BinTable *table);
static inline unsigned int hash_chrom_bin(const char *chrom, int bin, int capacity);
static inline BinEntry *bin_table_lookup(BinTable *table, const char *chrom, int bin);
static inline void bin_table_add(BinTable *table, const char *chrom, int bin);
static inline void bin_table_get_all_entries(BinTable *table, BinEntry **out_entries, int *out_size);

/* Hash function for chrom+bin combination */
static inline unsigned int hash_chrom_bin(const char *chrom, int bin, int capacity)
{
    unsigned int hash = 5381;

    /* Hash the chromosome string */
    const char *str = chrom;
    while (*str)
    {
        hash = ((hash << 5) + hash) + (*str++); /* hash * 33 + c */
    }

    /* Incorporate bin number */
    hash = ((hash << 5) + hash) + bin;

    return hash % capacity;
}

/* Implementations */
static inline BinTable *create_bin_table(int initial_capacity)
{
    BinTable *table = (BinTable *)R_Calloc(1, BinTable);
    table->buckets = (BinEntry **)R_Calloc(initial_capacity, BinEntry *);
    table->capacity = initial_capacity;
    table->size = 0;
    return table;
}

static inline void free_bin_table(BinTable *table)
{
    if (table)
    {
        /* Free all entries in all buckets */
        for (int i = 0; i < table->capacity; i++)
        {
            BinEntry *entry = table->buckets[i];
            while (entry)
            {
                BinEntry *next = entry->next;
                R_Free(entry->chrom);
                R_Free(entry);
                entry = next;
            }
        }
        R_Free(table->buckets);
        R_Free(table);
    }
}

static inline BinEntry *bin_table_lookup(BinTable *table, const char *chrom, int bin)
{
    unsigned int hash = hash_chrom_bin(chrom, bin, table->capacity);
    BinEntry *entry = table->buckets[hash];

    while (entry)
    {
        if (entry->bin == bin && strcmp(entry->chrom, chrom) == 0)
        {
            return entry;
        }
        entry = entry->next;
    }

    return NULL;
}

static inline void bin_table_add(BinTable *table, const char *chrom, int bin)
{
    BinEntry *existing = bin_table_lookup(table, chrom, bin);

    if (existing)
    {
        existing->count++;
    }
    else
    {
        /* Create new entry */
        BinEntry *new_entry = (BinEntry *)R_Calloc(1, BinEntry);
        size_t chrom_len = strlen(chrom) + 1;
        new_entry->chrom = (char *)R_Calloc(chrom_len, char);
        strcpy(new_entry->chrom, chrom);
        new_entry->bin = bin;
        new_entry->count = 1;

        /* Insert at head of bucket */
        unsigned int hash = hash_chrom_bin(chrom, bin, table->capacity);
        new_entry->next = table->buckets[hash];
        table->buckets[hash] = new_entry;

        table->size++;
    }
}

/* Get all entries as an array (for R export) */
static inline void bin_table_get_all_entries(BinTable *table, BinEntry **out_entries, int *out_size)
{
    *out_entries = (BinEntry *)R_Calloc(table->size, BinEntry);
    *out_size = 0;

    for (int i = 0; i < table->capacity; i++)
    {
        BinEntry *entry = table->buckets[i];
        while (entry)
        {
            (*out_entries)[*out_size].chrom = entry->chrom; /* Share pointer */
            (*out_entries)[*out_size].bin = entry->bin;
            (*out_entries)[*out_size].count = entry->count;
            (*out_entries)[*out_size].next = NULL;
            (*out_size)++;
            entry = entry->next;
        }
    }
}

#endif /* BIN_TABLE_H */
