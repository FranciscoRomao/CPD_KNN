MAIN = ballAlg.c

TARGET = ballAlg

OBJS = 

FLAGS = -O3 -fopenmp
debugFLAGS = -fopenmp

all: ballAlg

$(TARGET): gen_points.o
	gcc $(MAIN) $(debugFLAGS) $^ -o $@

gen_points.o: gen_points.c
	gcc $^ -c $(debugFLAGS)

# $^ - representa as dependencias do comando
# $@ - representa o ficheiro de saida do comando

.PHONY: clean #Serve para indicar que o comando "make clean" não tem nenhum ficheiro em específico

clean:
	rm -f $(wildcard *.o) $(TARGET)
