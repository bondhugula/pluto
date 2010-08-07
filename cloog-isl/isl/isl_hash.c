#include <stdlib.h>
#include "isl_hash.h"
#include "isl_ctx.h"

uint32_t isl_hash_string(uint32_t hash, const char *s)
{
	for (; *s; s++)
		isl_hash_byte(hash, *s);
	return hash;
}

static unsigned int round_up(unsigned int v)
{
	int old_v = v;

	while (v) {
		old_v = v;
		v ^= v & -v;
	}
	return old_v << 1;
}

int isl_hash_table_init(struct isl_ctx *ctx, struct isl_hash_table *table,
			int min_size)
{
	size_t size;

	if (!table)
		return -1;

	if (min_size < 2)
		min_size = 2;
	table->bits = ffs(round_up(4 * (min_size + 1) / 3 - 1)) - 1;
	table->n = 0;

	size = 2 << table->bits;
	table->entries = isl_calloc_array(ctx, struct isl_hash_table_entry,
					  size);
	if (!table->entries)
		return -1;

	return 0;
}

struct isl_hash_table *isl_hash_table_alloc(struct isl_ctx *ctx, int min_size)
{
	struct isl_hash_table *table = NULL;

	table = isl_alloc_type(ctx, struct isl_hash_table);
	if (isl_hash_table_init(ctx, table, min_size))
		goto error;
	return table;
error:
	isl_hash_table_free(ctx, table);
	return NULL;
}

void isl_hash_table_clear(struct isl_hash_table *table)
{
	if (!table)
		return;
	free(table->entries);
}

void isl_hash_table_free(struct isl_ctx *ctx, struct isl_hash_table *table)
{
	if (!table)
		return;
	isl_hash_table_clear(table);
	free(table);
}

struct isl_hash_table_entry *isl_hash_table_find(struct isl_ctx *ctx,
				struct isl_hash_table *table,
				uint32_t key_hash,
				int (*eq)(const void *entry, const void *val),
				const void *val, int reserve)
{
	size_t size;
	uint32_t h, key_bits;

	key_bits = isl_hash_bits(key_hash, table->bits);
	size = 2 << table->bits;
	for (h = key_bits; table->entries[h].data; h = (h+1) % size)
		if (table->entries[h].hash == key_hash &&
		    eq(table->entries[h].data, val))
			return &table->entries[h];

	if (!reserve)
		return NULL;

	assert(4 * table->n < 3 * size); /* XXX */
	table->n++;
	table->entries[h].hash = key_hash;

	return &table->entries[h];
}

void isl_hash_table_remove(struct isl_ctx *ctx,
				struct isl_hash_table *table,
				struct isl_hash_table_entry *entry)
{
	int h, h2, last_h;
	size_t size;

	if (!table || !entry)
		return;

	size = 2 << table->bits;
	h = entry - table->entries;
	isl_assert(ctx, h >= 0 && h < size, return);

	for (h2 = h+1; table->entries[h2 % size].data; h2++) {
		uint32_t bits = isl_hash_bits(table->entries[h2 % size].hash,
						table->bits);
		uint32_t offset = (size + bits - (h+1)) % size;
		if (offset <= h2 - (h+1))
			continue;
		*entry = table->entries[h2 % size];
		h = h2;
		entry = &table->entries[h % size];
	}

	entry->hash = 0;
	entry->data = NULL;
	table->n--;
}
