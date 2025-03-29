NTL = -lntl -lgmp -L/usr/local/lib
OBJ = GG.o GGFactoring.o ZZlib.o ZZFactoring.o mpqs.o rho.o

example: example.o QrtRootMod.o $(OBJ)
	g++ example.o QrtRootMod.o $(OBJ) $(NTL)
fig1: fig1.o $(OBJ)
	g++ fig1.o $(OBJ) $(NTL)
