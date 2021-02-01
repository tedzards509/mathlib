emcc main.cpp -O3 -ffast-math -s WASM=1 -o hello.html -s NO_EXIT_RUNTIME=1
emcc samples/mandelbrot.cpp -O3 -ffast-math -s WASM=1 -o mandelbrot.html -s NO_EXIT_RUNTIME=1