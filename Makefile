make: src/calcLtcs.c 
	mkdir bin
	gcc -o bin/calcLtcs src/calcLtcs.c -pthread -Wall -O3
