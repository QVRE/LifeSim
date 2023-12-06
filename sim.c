#include "sim.h"

double mutation_rate = (1./2048);
double death_rate = (1./1500);

int check_boundary(grid *g, u32 x, u32 y) {
	return x >= g->x || y >= g->y;
}

int check_blocked(grid *g, u32 x, u32 y) {
	if (check_boundary(g,x,y)) return 1;
	cell *c = &g->c[y*g->x+x];
	return (c->state >> 6 > 0);
}

int check_alive(grid *g, u32 x, u32 y) {
	if (check_boundary(g,x,y)) return 0;
	cell *c = &g->c[y*g->x+x];
	return c->state >> 7;
}

void kill_cell(grid *g, u32 x, u32 y) {
	if (!check_alive(g,x,y)) return; //if not alive, return
	cell *c = &g->c[y*g->x+x];
	c->state &= 0b00111110; //remove last two bits and lock
	c->state |= 0b01000000; //set type to dead
	u64 e = ((u64)y << 32) | x;
	for (u32 i = 0; i < g->exec_size; i++)
		if (g->exec[i] == e) {
			if (i < g->exec_size - 1) g->exec[i] = g->exec[g->exec_size - 1];
			g->exec_size--;
			return;
		}
}

void cell_rep_exec(grid *g, u32 x, u32 y, u32 new_x, u32 new_y) {
	u64 e = ((u64)y << 32) | x;
	for (u32 i = 0; i < g->exec_size; i++)
		if (g->exec[i] == e) {
			g->exec[i] = ((u64)new_y << 32) | new_x;
			return;
		}
}

void cell_copy_code(cell *c, cell *v) {
	for (u32 j = 0; j < GENE_LEN; j++) {
		v->code.c[j] = 0;
		for (u32 k = 0; k < 8; k++) {
			u8 m = (((double)rand() / RAND_MAX) < mutation_rate);
			v->code.c[j] |= (((c->code.c[j] >> k) & 1) ^ m) << k;
		}
	}
}

int cell_run(grid *g, u32 x, u32 y) {
	//return if cell is dead, non-existent, outside boundary or has no energy
	if (!check_alive(g,x,y)) return 1;
	cell *c = &g->c[y*g->x+x];
	u32 i = 0, run = 0, ran = 0;
	u32 type = c->state >> 6, lock = c->state & 1;
	u32 dir = (c->state >> 4) & 3;
	int dir_x = (int)(dir != 2) - (int)(dir != 0);
	int dir_y = (int)(dir != 3) - (int)(dir != 1);
	//process photosynthesis
	if (type == 2) {
		int product = c->chem[0] >> 5;
		if ((c->chem[0] & 31) > 8) product++;
		product -= c->chem[1] >> 6;
		if (product < 0) product = 0;
		if (255 - c->energy > product) c->energy += product;
	}
	//run genetic code
	do {
		u8 ins = c->code.c[i++];
		if (ins >> 7) { //action
			if (run) { //if 'promoter' has been triggered
				ran = 1;
				ins &= 0x7F;
				if (ins < 8) { //new cell
					if (check_blocked(g, x+dir_x, y+dir_y)) continue;
					if (g->exec_size >= g->exec_alloc) continue;
					//energy given to new cell
					u32 given = 1 << (2 * ((ins >> 1) & 3) + 1);
					u32 cost = given + 2;
					if (c->energy < cost) continue;
					c->energy -= cost;
					cell *new_c = &g->c[(y+dir_y)*g->x + x+dir_x];
					for (u32 j = 0; j < GENE_LEN; j++) new_c->code.c[j] = 0;
					new_c->energy = given;
					new_c->state = (1 << 7) | ((ins & 1) << 6) | (dir << 4);
					g->exec[g->exec_size++] = ((u64)(y+dir_y) << 32) | (x+dir_x);
				} else if (ins < 16) { //write code
					if (check_boundary(g, x+dir_x, y+dir_y)) continue;
					if (!check_blocked(g, x+dir_x, y+dir_y)) continue;
					cell *v = &g->c[(y+dir_y)*g->x + x+dir_x];
					if (v->state & 1) continue;
					cell_copy_code(c,v);
				} else if (ins < 24) { //currently NOP
					continue;
				} else if (ins < 32) { //turn
					dir += (ins & 1) ? 1 : -1;
					dir = dir & 3;
					dir_x = (int)(dir != 2) - (int)(dir != 0);
					dir_y = (int)(dir != 3) - (int)(dir != 1);
					c->state &= 0b11001111;
					c->state |= dir << 4;
				} else if (ins < 40) { //lock self
					lock = 1;
					c->state |= 1;
				} else if (ins < 48) { //unlock self
					lock = 0;
					c->state &= 0b11111110; 
				} else if (ins < 56) { //lock other cell
					if (!check_alive(g, x+dir_x, y+dir_y)) continue;
					g->c[(y+dir_y)*g->x + x+dir_x].state |= 1;
				} else if (ins < 64) { //push
					if (check_boundary(g, x+dir_x, y+dir_y)) continue;
					if (!check_blocked(g, x+dir_x, y+dir_y)) continue;
					if (check_boundary(g, x+dir_x*2, y+dir_y*2)) continue;
					if (check_blocked(g, x+dir_x*2, y+dir_y*2)) continue;
					cell nul = {(gene){{0}}, {0}, 0, 0};
					g->c[(y+dir_y*2)*g->x + x+dir_x*2] = g->c[(y+dir_y)*g->x + x+dir_x];
					g->c[(y+dir_y)*g->x + x+dir_x] = nul;
					cell_rep_exec(g, x+dir_x, y+dir_y, x+dir_x*2, y+dir_y*2);
				} else if (ins < 72) { //pull
					if (check_boundary(g, x+dir_x, y+dir_y)) continue;
					if (check_blocked(g, x+dir_x, y+dir_y)) continue;
					if (check_boundary(g, x+dir_x*2, y+dir_y*2)) continue;
					if (!check_blocked(g, x+dir_x*2, y+dir_y*2)) continue;
					cell nul = {(gene){{0}}, {0}, 0, 0};
					g->c[(y+dir_y)*g->x + x+dir_x] = g->c[(y+dir_y*2)*g->x + x+dir_x*2];
					g->c[(y+dir_y*2)*g->x + x+dir_x*2] = nul;
					cell_rep_exec(g, x+dir_x*2, y+dir_y*2, x+dir_x, y+dir_y);
				} else if (ins < 80) { //give energy
					if (!check_alive(g, x+dir_x, y+dir_y)) continue;
					cell *v = &g->c[(y+dir_y)*g->x + x+dir_x];
					u8 a = 1 << (ins & 3);
					if (c->energy < a || 255 - v->energy < a) continue;
					c->energy -= a;
					v->energy += a;
				} else if (ins < 88) { //take energy
					if (!check_alive(g, x+dir_x, y+dir_y)) continue;
					cell *v = &g->c[(y+dir_y)*g->x + x+dir_x];
					if (v->state & 1) continue;
					u8 a = 1 << (ins & 3);
					if (v->energy < a || 255 - c->energy < a) continue;
					v->energy -= a;
					c->energy += a;
				} else if (ins < 96) { //make chem
					u8 chem = ins & 3;
					if (c->energy == 0) continue;
					c->energy--;
					if (c->chem[chem] < 255) c->chem[chem]++;
				} else if (ins < 104) { //give chem
					if (!check_alive(g, x+dir_x, y+dir_y)) continue;
					cell *v = &g->c[(y+dir_y)*g->x + x+dir_x];
					if (v->state & 1) continue;
					u8 chem = ins & 3;
					if (c->chem[chem] == 0 || v->chem[chem] == 255) continue;
					c->chem[chem]--;
					v->chem[chem]++;
				} else if (ins < 112) { //digest chem
					u8 chem = ins & 3;
					if (c->chem[chem] == 0) continue;
					c->chem[chem]--;
					if (c->energy < 255) c->energy++;
				} else if (ins < 120) { //copy code
					if (check_boundary(g, x+dir_x, y+dir_y)) continue;
					if (!check_blocked(g, x+dir_x, y+dir_y)) continue;
					if (check_boundary(g, x-dir_x, y-dir_y)) continue;
					if (!check_blocked(g, x-dir_x, y-dir_y)) continue;
					cell *v1 = &g->c[(y+dir_y)*g->x + x+dir_x];
					if (v1->state & 1) continue;
					cell *v2 = &g->c[(y-dir_y)*g->x + x-dir_x];
					if (v2->state & 1) continue;
					cell_copy_code(v1,v2);
				} else if (ins < 128) { //digest cell
					if (check_boundary(g, x+dir_x, y+dir_y)) continue;
					cell *v = &g->c[(y+dir_y)*g->x + x+dir_x];
					if (!c->chem[1] || (v->state >> 6) == 3) continue;
					kill_cell(g, x+dir_x, y+dir_y);
					u32 gain = (u32)v->energy / 2 + 1;
					gain += (u32)v->chem[0] / 2;
					gain += (u32)v->chem[1] / 2;
					gain += (u32)v->chem[2] / 2;
					gain += (u32)v->chem[3] / 2;
					c->chem[1]--;
					if (gain > 255 - c->energy) {
						c->energy = 255;
					} else c->energy += gain;
					//delete dead cell
					v->state = 0;
					v->energy = 0;
					v->chem[0] = 0, v->chem[1] = 0, v->chem[2] = 0, v->chem[3] = 0;
					for (u32 j = 0; j < GENE_LEN; j++) v->code.c[j] = 0;
				}
			}
		} else { //promoter
			if (run ^ (ins >> 6)) continue;
			u8 r = !(ins >> 6);
			u8 ptype = (ins >> 3) & 7;
			u8 parg = ins & 7;
			switch (ptype) {
				case 0: //time
					parg = parg * 2 + 1;
					if (parg <= 1) {
						run = r;
						break;
					}
					if (((double)rand() / RAND_MAX) < (1. / (double)parg)) run = r;
					break;
				case 1: //energy
					u8 p = c->energy;
					goto promoter_presence;
				case 2: //chem 0
					p = c->chem[0];
					goto promoter_presence;
				case 3: //chem 1
					p = c->chem[1];
					goto promoter_presence;
				case 4: //chem 2
					p = c->chem[2];
					goto promoter_presence;
				case 5: //chem 3
					p = c->chem[3];
					goto promoter_presence;
				case 6: //cell perception
					if (check_boundary(g, x+dir_x, y+dir_y)) continue;
					cell *v = &g->c[(y+dir_y)*g->x + x+dir_x];
					u8 vtype = v->state >> 6, vlock = v->state & 1;
					switch (parg) {
						case 0: if (vlock) run = r; break;
						case 1: if (!vlock) run = r; break;
						case 2: if (vtype == 0) run = r; break;
						case 3: if (vtype == 1) run = r; break;
						case 4: if (vtype == 2) run = r; break;
						case 5: if (vtype == 3) run = r; break;
						case 6: if (type == 2) run = r; break;
						case 7: if (type == 3) run = r; break;
					}
					break;
promoter_presence:
				switch (parg) {
					case 0: if (p) run = r; break;
					case 1: if (!p) run = r; break;
					case 2: if (p < 32) run = r; break;
					case 3: if (p >= 16) run = r; break;
					case 4: if (p >= 64) run = r; break;
					case 5: if (p >= 192) run = r; break;
					case 6: if (p >= 16 && p < 64) run = r; break;
					case 7: if (p >= 64 && p < 192) run = r; break;
				}
				break;
			}
		}
	} while(i < GENE_LEN);
	if (c->energy == 0 || ((double)rand() / RAND_MAX) < death_rate) {
		kill_cell(g,x,y);
		return 1;
	} else c->energy -= ran;
	return 0;
}
