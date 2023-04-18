g++ -O3 -std=c++11 main_exact.cpp -o exact;
./exact none; # h-motif counting
./exact ab1; # 3h-motif counitng
rm exact;
