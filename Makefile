CC=gcc -lm -O2 -Wall

DEMOS=tetra sphere dna

$(DEMOS): %: %.c renderer.o phong.h
	$(CC) $< renderer.o -o $@

renderer.o: renderer.c renderer.h constants.h
	$(CC) -c renderer.c

clean:
	rm -f renderer.o $(DEMOS)

