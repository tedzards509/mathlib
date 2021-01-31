CXX = clang++

CXX_ARGS = -O3 -ffast-math

mandelbrot: mandelbrot.cpp
	$(CXX) $(CXX_ARGS) mandelbrot.cpp -o mandelbrot

mandelbrot_win: mandelbrot.cpp
	x86_64-w64-mingw32-g++ -O3 -ffast-math mandelbrot.cpp -o mandelbrot.exe -march=native  -static-libstdc++ -static-libgcc -static

all: mandelbrot