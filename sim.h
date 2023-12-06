#pragma once
#include "evar.h"
#include <stdlib.h>

#define GENE_LEN 24

typedef struct {
	u8 c[GENE_LEN];
} gene;

typedef struct {
	gene code;
	u8 chem[4];
	u8 state; // type(2) | dir(2) | junk(3) | lock(1)
	u8 energy;
} cell;

typedef struct {
	cell *c;
	u64 *exec;
	u32 x, y;
	u32 exec_size, exec_alloc;
} grid;

extern double mutation_rate;
extern double death_rate;

void kill_cell(grid *g, u32 x, u32 y);

int cell_run(grid *g, u32 x, u32 y);
