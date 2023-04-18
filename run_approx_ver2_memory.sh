g++ -O3 -std=c++11 main_approx_ver2_memory.cpp -o approx_ver2_memory;
./approx_ver2_memory 10000 0.01 none; # h-motif counting
./approx_ver2_memory 10000 0.01 ab1; # 3h-motif counitng
rm approx_ver2_memory;
