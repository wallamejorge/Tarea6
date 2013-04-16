a.out: main.c functions.c
	cc -lm main.c functions.c
	./a.out 
clean:
	rm -f *.o 
	rm -f *.txt 
