
g++ -O3 -march=native -std=c++20 -mavx512f -ffast-math -fno-math-errno -mprefer-vector-width=512 block_bench.cpp -lm -o sw_bench -save-temps && ./sw_bench
