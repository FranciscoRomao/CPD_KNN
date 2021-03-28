MAIN = ballAlg.c

TARGET = ballAlg

OBJS = 

FLAGS = -O3 -fopenmp -lm
debugFLAGS = -fopenmp -g -lm -Wall

all: ballAlg

$(TARGET): gen_points.o geometry.o
	gcc $(MAIN)  $^ -o $@ $(debugFLAGS)

gen_points.o: gen_points.c
	gcc $^ -c $(debugFLAGS)

geometry.o:geometry.c
	gcc $^ -c $(debugFLAGS)

# $^ - representa as dependencias do comando
# $@ - representa o ficheiro de saida do comando

.PHONY: clean #Serve para indicar que o comando "make clean" não tem nenhum ficheiro em específico

clean:
	rm -f $(wildcard *.o) $(TARGET)
