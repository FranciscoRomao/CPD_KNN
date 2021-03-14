MAIN = ballAlg.c

TARGET = ballAlg

OBJS = 

FLAGS = -O3 -fopenmp -lm
debugFLAGS = -fopenmp -lm

all: ballAlg

$(TARGET): gen_points.o distance.o
	gcc $(MAIN)  $^ -o $@ $(debugFLAGS)

gen_points.o: gen_points.c
	gcc $^ -c $(debugFLAGS)

distance.o:distance.c
	gcc $^ -c $(debugFLAGS)

# $^ - representa as dependencias do comando
# $@ - representa o ficheiro de saida do comando

.PHONY: clean #Serve para indicar que o comando "make clean" não tem nenhum ficheiro em específico

clean:
	rm -f $(wildcard *.o) $(TARGET)
