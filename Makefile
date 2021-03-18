# Makefile for Project 3
#


# List all the drivers
all: main


# list the test program all classes (cpp files)
#  cpp files must use #include , to include *.h or *.hpp files

main: main.cpp Linear_solve.hpp linalg.cpp
	g++ -std=c++11 -o main.o main.cpp Linear_solve.hpp linalg.cpp

# List all the executables under 'all:'
clean:
	rm main.o
