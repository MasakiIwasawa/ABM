#PS_PATH = ../../FDPS/src/
PS_PATH = ../../fdps-devel/src/ 
INC = -I$(PS_PATH)

#CC = time g++
CC = time mpicxx
CFLAGS = -std=c++17
CFLAGS += -Wall
#CFLAGS = -O3
#CFLAGS += -ffast-math
#CFLAGS += -funroll-loops
CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

all: main.out

main.out:main.cpp
	$(CC) $(INC) $(CFLAGS) -o $@ $^
clean:
	rm main.out

