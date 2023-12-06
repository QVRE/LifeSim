libs = -lm -lSDL2
obj = sim.o Elafri.o graphics.o gui.o pkg.o

main: main.c $(obj)
	$(CC) $(libs) -O3 $(obj) main.c -o out

$(obj): %.o: %.c
	$(CC) -O3 -c $< -o $@

clean:
	rm out $(obj)
