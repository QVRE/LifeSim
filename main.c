#include "Elafri.h"
#include "gui.h"
#include "pkg.h"
#include "sim.h"
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define RES_X 600
#define RES_Y 600
#define FPS 60

int main() {
	ElafriInit("Sim", RES_X, RES_Y);

	gr Gr = GrBuffer(res.x, res.y);

	//make grid
	grid g;
	g.c = calloc(RES_X*RES_Y, sizeof(cell));
	g.x = RES_X, g.y = RES_Y;
	g.exec_size = 0, g.exec_alloc = 65536;
	g.exec = malloc(g.exec_alloc * sizeof(u64));

	//variables
	u32 mode = 0, show_pop = 1, mouse_click = 1, zoom = 1, running = 1;
	int scx = 0, scy = 0;
	float dscx = 0, dscy = 0, dtime = 0;
	cell stored_cell;
	srand(time(NULL));

	//starting cell
	u32 startx = RES_X / 2, starty = RES_Y / 2;
	cell *start_cell = &g.c[starty * g.x + startx];
	start_cell->state = 2 << 6;
	start_cell->energy = 100;
	start_cell->chem[0] = 64;
	{
	u32 cii = 0; //counter
	start_cell->code.c[cii++] = 7; //time promoter
	start_cell->code.c[cii++] = 220; //make chem 0
	start_cell->code.c[cii++] = 132; //new cell
	start_cell->code.c[cii++] = 138; //copy over code
	start_cell->code.c[cii++] = 158; //turn
	start_cell->code.c[cii++] = 253; //digest cell
	start_cell->code.c[cii++] = 166; //lock self
	start_cell->code.c[cii++] = 64; //stop
	}
	//append organism to execution list
	g.exec[g.exec_size++] = ((u64)starty << 32) | startx;

	RenderLoop:
	Input(); //https://wiki.libsdl.org/SDL_Keycode
	if (kbd[SDL_SCANCODE_M] == KEY_PRESSED) mode = (mode + 1) & 3;
	if (kbd[SDL_SCANCODE_P] == KEY_PRESSED) show_pop = !show_pop;
	//camera movement code
	dscx = fmod(dscx, 1), dscy = fmod(dscy, 1), dtime *= 64;
	if (kbd[SDL_SCANCODE_W]) dscy -= dtime;
	if (kbd[SDL_SCANCODE_A]) dscx -= dtime;
	if (kbd[SDL_SCANCODE_S]) dscy += dtime;
	if (kbd[SDL_SCANCODE_D]) dscx += dtime;
	scx += (int)dscx, scy += (int)dscy;
	if (kbd[SDL_SCANCODE_Q] == KEY_PRESSED && zoom > 1) zoom /= 2, scx -= RES_X / 4 / zoom, scy -= RES_Y / 4 / zoom;
	if (kbd[SDL_SCANCODE_E] == KEY_PRESSED && zoom < 32) scx += RES_X / 4 / zoom, scy += RES_Y / 4 / zoom, zoom *= 2;
	//fix values so we don't get a segmentation fault
	if (scx < 0) scx = 0;
	if (scy < 0) scy = 0;
	if (scx + RES_X / zoom >= RES_X) scx = RES_X - RES_X / zoom;
	if (scy + RES_Y / zoom >= RES_Y) scy = RES_Y - RES_Y / zoom;
	//save / load world
	if (kbd[SDL_SCANCODE_K] == KEY_PRESSED) {
		pfile savedat[2] = {{0, g.c, g.x*g.y*sizeof(cell), 1}, {1, g.exec, g.exec_size*sizeof(u64), 1}};
		pkg p = {savedat, 2};
		WritePkg("world.bin", p, 2);
	}
	if (kbd[SDL_SCANCODE_L] == KEY_PRESSED) {
		pkg p = ReadPkg("world.bin");
		//not a full check but good enough to assume the file is fine
		if (p.filecount != 2) goto load_world_error;
		if (p.files[0].ID != 0 || p.files[1].ID != 1) goto load_world_error;
		if (p.files[0].size != g.x * g.y * sizeof(cell)) goto load_world_error;
		memcpy(g.c, p.files[0].dat, g.x * g.y * sizeof(cell));
		memcpy(g.exec, p.files[1].dat, p.files[1].size);
		g.exec_size = p.files[1].size / sizeof(u64);
		FreePkg(p);
	}
load_world_error:
	if (kbd[SDL_SCANCODE_SPACE] == KEY_PRESSED) running = !running;
	if (kbd[SDL_SCANCODE_X] == KEY_PRESSED) running = 2;
	if (mC & SDL_BUTTON_LMASK && mouse_click) {
		printf("BEGIN GENOME:\n");
		u64 pos = (scy + m.y / zoom) * g.x + scx + m.x / zoom;
		for (u32 i = 0; i < GENE_LEN; i++) {
			printf("%u\n", g.c[pos].code.c[i]);
		}
		printf("END\n");
		mouse_click = 0;
	} else mouse_click = ~mC & SDL_BUTTON_LMASK;

	//run execution list
	if (running) {
		for (u32 i = 0; i < g.exec_size; i++) {
			u64 e = g.exec[i];
			u32 ex = e & 0xFFFFFFFF;
			u32 ey = e >> 32;
			cell_run(&g, ex, ey);
		}
		running = (running != 2); //set running to 0 if it's 2 else 1
	}

	if (mode == 0) {
		for (u32 i = 0; i < RES_Y; i++)
			for (u32 j = 0; j < RES_X; j++) {
				u32 iy = scy + i / zoom, jx = scx + j / zoom;
				u64 pos = iy * g.x + jx;
				switch (g.c[pos].state >> 6) {
					case 0: break;
					case 1:
						color clr = WHITE;
						GrPixel(&Gr, j, i, clr);
						break;
					case 2:
						u8 sum1 = 0, sum2 = 0, sum3 = 0;
						sum1 = g.c[pos].code.c[0] ^ g.c[pos].code.c[1] ^ g.c[pos].code.c[2];
						sum2 = g.c[pos].code.c[3] ^ g.c[pos].code.c[4] ^ g.c[pos].code.c[5];
						sum3 = g.c[pos].code.c[6] ^ g.c[pos].code.c[7] ^ g.c[pos].code.c[8];
						clr = (color){sum1, 192 + (sum2 / 4), sum3};
						GrPixel(&Gr, j, i, clr);
						break;
					case 3:
						sum1 = 0, sum2 = 0, sum3 = 0;
						sum1 = g.c[pos].code.c[0] ^ g.c[pos].code.c[1] ^ g.c[pos].code.c[2];
						sum2 = g.c[pos].code.c[3] ^ g.c[pos].code.c[4] ^ g.c[pos].code.c[5];
						sum3 = g.c[pos].code.c[6] ^ g.c[pos].code.c[7] ^ g.c[pos].code.c[8];
						clr = (color){192 + (sum1 / 4), sum2, sum3};
						GrPixel(&Gr, j, i, clr);
						break;
				}
			}
	} else {
		for (u32 i = 0; i < RES_Y; i++)
			for (u32 j = 0; j < RES_X; j++) {
				u32 iy = scy + i / zoom, jx = scx + j / zoom;
				u64 pos = iy * g.x + jx;
				u8 val = 0;
				switch (mode) {
					case 0: break;
					case 1: val = g.c[pos].energy; break;
					case 2: val = g.c[pos].chem[0]; break;
					case 3: val = g.c[pos].chem[1]; break;
				}
				color clr = {(mode == 3) * val, (mode == 2) * val, (mode == 1) * val};
				GrPixel(&Gr, j, i, clr);
			}
	}

	char execstr[32];
	if (show_pop) {
		sprintf(execstr, "%u", g.exec_size);
		GrText(&Gr, 16, 16, strlen(execstr)*24, 24, strlen(execstr), BLUE, execstr);
	}
	DrawFill(&Gr, BLACK);

	dtime = FramerateHandler(FPS);
	goto RenderLoop;
}
