g++ -O3 -std=c++11 main_approx_ver1.cpp -o approx_ver1;
./approx_ver1 1500 none; # h-motif counting
./approx_ver1 1500 ab1; # 3h-motif counting
rm approx_ver1;
