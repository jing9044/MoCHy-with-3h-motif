g++ -O3 -std=c++11 -lgomp -fopenmp main_exact_par.cpp -o exact_par;
./exact_par 4 none; # h-motif counting
./exact_par 4 ab1; # 3h-motif counitng
rm exact_par;
