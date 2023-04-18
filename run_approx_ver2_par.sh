g++ -O3 -std=c++11 -lgomp -fopenmp main_approx_ver2_par.cpp -o approx_ver2_par;
./approx_ver2_par 10000 4 none; # h-motif counting
./approx_ver2_par 10000 4 ab1; # 3h-motif counting
rm approx_ver2_par;
