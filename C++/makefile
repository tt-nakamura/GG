NTL = -lntl -lgmp -L/usr/local/lib
OBJ = GG.o GGFactoring.o ZZlib.o ZZFactoring.o mpqs.o rho.o

example: example.o $(OBJ)
	g++ example.o $(OBJ) $(NTL)
