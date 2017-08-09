renderer: renderer.c
	gcc -o renderer -lm -O2 -Wall renderer.c

clean:
	rm -f renderer

