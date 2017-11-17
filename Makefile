CC=gcc -lm -O2 -Wall

DEMOS=tetra sphere dna

$(DEMOS): %: %.c renderer.o
	$(CC) $< renderer.o -o $@

renderer.o: renderer.c
	$(CC) -c renderer.c

clean:
	rm -f renderer.o $(DEMOS)

