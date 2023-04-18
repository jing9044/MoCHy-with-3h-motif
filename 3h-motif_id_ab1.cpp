#include <iostream>

using namespace std;

pair <int, int> get_motif_index_ab1(int deg_a, int deg_b, int deg_c, int C_ab, int C_bc, int C_ca, int g_abc){
	int a = deg_a - (C_ab + C_ca) + g_abc; 
	int b = deg_b - (C_bc + C_ab) + g_abc; 
	int c = deg_c - (C_ca + C_bc) + g_abc; 
	int d = C_ab - g_abc; 
	int e = C_bc - g_abc; 
	int f = C_ca - g_abc; 
	int g = g_abc; 
	int motif_id = (a > 0) + ((b > 0) << 1) + ((c > 0) << 2) + ((d > 0) << 3) + ((e > 0) << 4) + ((f > 0) << 5) + ((g > 0) << 6);
	int absolute_motif_id = (a == 1) + ((b == 1) << 1) + ((c == 1) << 2) + ((d == 1) << 3) + ((e == 1) << 4) + ((f == 1) << 5) + ((g == 1) << 6);
	
	if (id_to_index[motif_id] - 2 == 1) { // h-motif 1
		if (absolute_motif_id == 67 || absolute_motif_id == 69 || absolute_motif_id == 70) // 1000011 1000101 1000110 : 1 0 2
			return {1, 0};
		else if (absolute_motif_id == 3 || absolute_motif_id == 5 || absolute_motif_id == 6) // 0000011 0000101 0000110 : 0 0 2
			return {1, 1};
		else if (absolute_motif_id == 65 || absolute_motif_id == 66 || absolute_motif_id == 68) // 1000001 1000010 1000100 : 1 0 1
			return {1, 2};
		else if (absolute_motif_id == 1 || absolute_motif_id == 2 || absolute_motif_id == 4) // 0000001 0000010 0000100 : 0 0 1
			return {1, 3};
		else if (absolute_motif_id == 64) // 1000000 : 1 0 0
			return {1, 4};
		else if (absolute_motif_id == 0) // 0000000 : 0 0 0
			return {1, 5};
	}
	else if (id_to_index[motif_id] - 2 == 2) { // h-motif 2
		if (absolute_motif_id == 71) // 1000111 : 1 0 3
			return {2, 0};
		else if (absolute_motif_id == 7) // 0000111 : 0 0 3
			return {2, 1};
		else if (absolute_motif_id == 67 || absolute_motif_id == 69 || absolute_motif_id == 70) // 1000011 1000101 1000110 : 1 0 2
			return {2, 2};
		else if (absolute_motif_id == 3 || absolute_motif_id == 5 || absolute_motif_id == 6) // 0000011 0000101 0000110 : 0 0 2
			return {2, 3};
		else if (absolute_motif_id == 65 || absolute_motif_id == 66 || absolute_motif_id == 68) // 1000001 1000010 1000100 : 1 0 1
			return {2, 4};
		else if (absolute_motif_id == 1 || absolute_motif_id == 2 || absolute_motif_id == 4) // 0000001 0000010 0000100 : 0 0 1
			return {2, 5};
		else if (absolute_motif_id == 64) // 1000000 : 1 0 0
			return {2, 6};
		else if (absolute_motif_id == 0) // 0000000 : 0 0 0
			return {2, 7};
	}
	else if (id_to_index[motif_id] - 3 == 3) { // h-motif 3
		if (absolute_motif_id == 73 || absolute_motif_id == 74 || absolute_motif_id == 76 || 
			absolute_motif_id == 81 || absolute_motif_id == 82 || absolute_motif_id == 84 || 
			absolute_motif_id == 97 || absolute_motif_id == 98 || absolute_motif_id == 100) // 1001001 1001010 1001100 1010001 1010010 1010100 1100001 1100010 1100100 : 1 1 1
			return {3, 0};
		else if (absolute_motif_id == 9 || absolute_motif_id == 10 || absolute_motif_id == 12 || 
				 absolute_motif_id == 17 || absolute_motif_id == 18 || absolute_motif_id == 20 || 
				 absolute_motif_id == 33 || absolute_motif_id == 34 || absolute_motif_id == 36) // 0001001 0001010 0001100 0010001 0010010 0010100 0100001 0100010 0100100 : 0 1 1
			return {3, 1};
		else if (absolute_motif_id == 65 || absolute_motif_id == 66 || absolute_motif_id == 68) // 1000001 1000010 1000100 : 1 0 1
			return {3, 2};
		else if (absolute_motif_id == 1 || absolute_motif_id == 2 || absolute_motif_id == 4) // 0000001 0000010 0000100 : 0 0 1
			return {3, 3};
		else if (absolute_motif_id == 72 || absolute_motif_id == 80 || absolute_motif_id == 96) // 1001000 1010000 1100000 : 1 1 0
			return {3, 4};
		else if (absolute_motif_id == 8 || absolute_motif_id == 16 || absolute_motif_id == 32) // 0001000 0010000 0100000 : 0 1 0
			return {3, 5};
		else if (absolute_motif_id == 64) // 1000000 : 1 0 0
			return {3, 6};
		else if (absolute_motif_id == 0) // 0000000 : 0 0 0
			return {3, 7};
	}
	else if (id_to_index[motif_id] - 4 == 4) { // h-motif 4
		if (absolute_motif_id == 75 || absolute_motif_id == 77 || absolute_motif_id == 78 || 
			absolute_motif_id == 83 || absolute_motif_id == 85 || absolute_motif_id == 86 || 
			absolute_motif_id == 99 || absolute_motif_id == 101 || absolute_motif_id == 102) // 1001011 1001101 1001110 1010011 1010101 1010110 1100011 1100101 1100110 : 1 1 2
			return {4, 0};
		else if (absolute_motif_id == 11 || absolute_motif_id == 13 || absolute_motif_id == 14 || 
				 absolute_motif_id == 19 || absolute_motif_id == 21 || absolute_motif_id == 22 || 
				 absolute_motif_id == 35 || absolute_motif_id == 37 || absolute_motif_id == 38) // 0001011 0001101 0001110 0010011 0010101 0010110 0100011 0100101 0100110 : 0 1 2
			return {4, 1};
		else if (absolute_motif_id == 67 || absolute_motif_id == 69 || absolute_motif_id == 70) // 1000011 1000101 1000110 : 1 0 2
			return {4, 2};
		else if (absolute_motif_id == 3 || absolute_motif_id == 5 || absolute_motif_id == 6) // 0000011 0000101 0000110 : 0 0 2
			return {4, 3};
		else if (absolute_motif_id == 73 || absolute_motif_id == 74 || absolute_motif_id == 76 || 
				 absolute_motif_id == 81 || absolute_motif_id == 82 || absolute_motif_id == 84 || 
				 absolute_motif_id == 97 || absolute_motif_id == 98 || absolute_motif_id == 100) // 1001001 1001010 1001100 1010001 1010010 1010100 1100001 1100010 1100100 : 1 1 1
			return {4, 4};
		else if (absolute_motif_id == 9 || absolute_motif_id == 10 || absolute_motif_id == 12 || 
				 absolute_motif_id == 17 || absolute_motif_id == 18 || absolute_motif_id == 20 || 
				 absolute_motif_id == 33 || absolute_motif_id == 34 || absolute_motif_id == 36) // 0001001 0001010 0001100 0010001 0010010 0010100 0100001 0100010 0100100 : 0 1 1
			return {4, 5};
		else if (absolute_motif_id == 65 || absolute_motif_id == 66 || absolute_motif_id == 68) // 1000001 1000010 1000100 : 1 0 1
			return {4, 6};
		else if (absolute_motif_id == 1 || absolute_motif_id == 2 || absolute_motif_id == 4) // 0000001 0000010 0000100 : 0 0 1
			return {4, 7};
		else if (absolute_motif_id == 72 || absolute_motif_id == 80 || absolute_motif_id == 96) // 1001000 1010000 1100000 : 1 1 0
			return {4, 8};
		else if (absolute_motif_id == 8 || absolute_motif_id == 16 || absolute_motif_id == 32) // 0001000 0010000 0100000 : 0 1 0
			return {4, 9};
		else if (absolute_motif_id == 64) // 1000000 : 1 0 0
			return {4, 10};
		else if (absolute_motif_id == 0) // 0000000 : 0 0 0
			return {4, 11};
	}
	else if (id_to_index[motif_id] - 4 == 5) { // h-motif 5 : not symmetric
		if (absolute_motif_id == 75 || absolute_motif_id == 77 || absolute_motif_id == 78 || 
			absolute_motif_id == 83 || absolute_motif_id == 85 || absolute_motif_id == 86 || 
			absolute_motif_id == 99 || absolute_motif_id == 101 || absolute_motif_id == 102) // 1001011 1001101 1001110 1010011 1010101 1010110 1100011 1100101 1100110 : 1 1 2
			return {5, 0};
		else if (absolute_motif_id == 11 || absolute_motif_id == 13 || absolute_motif_id == 14 || 
				 absolute_motif_id == 19 || absolute_motif_id == 21 || absolute_motif_id == 22 || 
				 absolute_motif_id == 35 || absolute_motif_id == 37 || absolute_motif_id == 38) // 0001011 0001101 0001110 0010011 0010101 0010110 0100011 0100101 0100110 : 0 1 2
			return {5, 1};
		else if (absolute_motif_id == 67 || absolute_motif_id == 69 || absolute_motif_id == 70) // 1000011 1000101 1000110 : 1 0 2
			return {5, 2};
		else if (absolute_motif_id == 3 || absolute_motif_id == 5 || absolute_motif_id == 6) // 0000011 0000101 0000110 : 0 0 2
			return {5, 3};
		else if (absolute_motif_id == 73 || absolute_motif_id == 74 || absolute_motif_id == 82 || 
				 absolute_motif_id == 84 || absolute_motif_id == 97 || absolute_motif_id == 100) // 1001001 1001010 1010010 1010100 1100001 1100100
			return {5, 4};
		else if (absolute_motif_id == 76 || absolute_motif_id == 81 || absolute_motif_id == 98) // 1001100 1010001 1100010
			return {5, 5};
		else if (absolute_motif_id == 9 || absolute_motif_id == 10 || absolute_motif_id == 18 || 
				 absolute_motif_id == 20 || absolute_motif_id == 33 || absolute_motif_id == 36) // 0001001 0001010 0010010 0010100 0100001 0100100
			return {5, 6};
		else if (absolute_motif_id == 12 || absolute_motif_id == 17 || absolute_motif_id == 34) // 0001100 0010001 0100010
			return {5, 7};
		else if (motif_id - absolute_motif_id == 12 || motif_id - absolute_motif_id == 17 || motif_id - absolute_motif_id == 34) // (motif_id - absolute_motif_id) 0001100 0010001 0100010
			return {5, 8};	
		else if (motif_id - absolute_motif_id == 9 || motif_id - absolute_motif_id == 10 || motif_id - absolute_motif_id == 18 || 
				 motif_id - absolute_motif_id == 20 || motif_id - absolute_motif_id == 33 || motif_id - absolute_motif_id == 36) // (motif_id - absolute_motif_id) 0001001 0001010 0010010 0010100 0100001 0100100
			return {5, 9};
		else if (motif_id - absolute_motif_id == 76 || motif_id - absolute_motif_id == 81 || motif_id - absolute_motif_id == 98) // (motif_id - absolute_motif_id) 1001100 1010001 1100010
			return {5, 10};
		else if (motif_id - absolute_motif_id == 73 || motif_id - absolute_motif_id == 74 || motif_id - absolute_motif_id == 82 || 
				 motif_id - absolute_motif_id == 84 || motif_id - absolute_motif_id == 97 || motif_id - absolute_motif_id == 100) // (motif_id - absolute_motif_id) 1001001 1001010 1010010 1010100 1100001 1100100
			return {5, 11};
		else if (absolute_motif_id == 72 || absolute_motif_id == 80 || absolute_motif_id == 96) // 1001000 1010000 1100000 : 1 1 0
			return {5, 12};
		else if (absolute_motif_id == 8 || absolute_motif_id == 16 || absolute_motif_id == 32) // 0001000 0010000 0100000 : 0 1 0
			return {5, 13};
		else if (absolute_motif_id == 64) // 1000000 : 1 0 0
			return {5, 14};
		else if (absolute_motif_id == 0) // 0000000 : 0 0 0
			return {5, 15};
	}
	else if (id_to_index[motif_id] - 4 == 6) { // h-motif 6 : not symmetric
		if (absolute_motif_id == 79 || absolute_motif_id == 87 || absolute_motif_id == 103)// 1001111 1010111 1100111 : 1 1 3
			return {6, 0};
		else if (absolute_motif_id == 15 || absolute_motif_id == 23 || absolute_motif_id == 39)// 0001111 0010111 0100111 : 0 1 3 
			return {6, 1};
		else if (absolute_motif_id == 71) // 1000111 : 1 0 3
			return {6, 2};
		else if (absolute_motif_id == 7) // 0000111 : 0 0 3
			return {6, 3};
		else if (absolute_motif_id == 75 || absolute_motif_id == 86 || absolute_motif_id == 101) // 1001011 1010110 1100101
			return {6, 4};
		else if (absolute_motif_id == 77 || absolute_motif_id == 78 || absolute_motif_id == 83 || 
				 absolute_motif_id == 85 || absolute_motif_id == 99 || absolute_motif_id == 102) // 1001101 1001110 1010011 1010101 1100011 1100110
			return {6, 5};
		else if (absolute_motif_id == 11 || absolute_motif_id == 22 || absolute_motif_id == 37) // 0001011 0010110 0100101
			return {6, 6};
		else if (absolute_motif_id == 13 || absolute_motif_id == 14 || absolute_motif_id == 19 || 
				 absolute_motif_id == 21 || absolute_motif_id == 35 || absolute_motif_id == 38) // 0001101 0001110 0010011 0010101 0100011 0100110
			return {6, 7};
		else if (motif_id - absolute_motif_id == 12 || motif_id - absolute_motif_id == 17 || motif_id - absolute_motif_id == 34) // (motif_id - absolute_motif_id) 0001100 0010001 0100010
			return {6, 8};	
		else if (motif_id - absolute_motif_id == 9 || motif_id - absolute_motif_id == 10 || motif_id - absolute_motif_id == 18 || 
				 motif_id - absolute_motif_id == 20 || motif_id - absolute_motif_id == 33 || motif_id - absolute_motif_id == 36) // (motif_id - absolute_motif_id) 0001001 0001010 0010010 0010100 0100001 0100100
			return {6, 9};
		else if (motif_id - absolute_motif_id == 76 || motif_id - absolute_motif_id == 81 || motif_id - absolute_motif_id == 98) // (motif_id - absolute_motif_id) 1001100 1010001 1100010
			return {6, 10};
		else if (motif_id - absolute_motif_id == 73 || motif_id - absolute_motif_id == 74 || motif_id - absolute_motif_id == 82 || 
				 motif_id - absolute_motif_id == 84 || motif_id - absolute_motif_id == 97 || motif_id - absolute_motif_id == 100) // (motif_id - absolute_motif_id) 1001001 1001010 1010010 1010100 1100001 1100100
			return {6, 11};
		else if (absolute_motif_id == 73 || absolute_motif_id == 74 || absolute_motif_id == 82 || 
				 absolute_motif_id == 84 || absolute_motif_id == 97 || absolute_motif_id == 100) // 1001001 1001010 1010010 1010100 1100001 1100100
			return {6, 12};
		else if (absolute_motif_id == 76 || absolute_motif_id == 81 || absolute_motif_id == 98) // 1001100 1010001 1100010
			return {6, 13};
		else if (absolute_motif_id == 9 || absolute_motif_id == 10 || absolute_motif_id == 18 || 
				 absolute_motif_id == 20 || absolute_motif_id == 33 || absolute_motif_id == 36) // 0001001 0001010 0010010 0010100 0100001 0100100
			return {6, 14};
		else if (absolute_motif_id == 12 || absolute_motif_id == 17 || absolute_motif_id == 34) // 0001100 0010001 0100010
			return {6, 15};
		else if (motif_id - absolute_motif_id == 13 || motif_id - absolute_motif_id == 14 || motif_id - absolute_motif_id == 19 || 
				 motif_id - absolute_motif_id == 21 || motif_id - absolute_motif_id == 35 || motif_id - absolute_motif_id == 38) // (motif_id - absolute_motif_id) 0001101 0001110 0010011 0010101 0100011 0100110
			return {6, 16};
		else if (motif_id - absolute_motif_id == 11 || motif_id - absolute_motif_id == 22 || motif_id - absolute_motif_id == 37) // (motif_id - absolute_motif_id) 0001011 0010110 0100101
			return {6, 17};
		else if (motif_id - absolute_motif_id == 77 || motif_id - absolute_motif_id == 78 || motif_id - absolute_motif_id == 83 || 
				 motif_id - absolute_motif_id == 85 || motif_id - absolute_motif_id == 99 || motif_id - absolute_motif_id == 102) // (motif_id - absolute_motif_id) 1001101 1001110 1010011 1010101 1100011 1100110
			return {6, 18};
		else if (motif_id - absolute_motif_id == 75 || motif_id - absolute_motif_id == 86 || motif_id - absolute_motif_id == 101) // (motif_id - absolute_motif_id) 1001011 1010110 1100101
			return {6, 19};
		else if (absolute_motif_id == 72 || absolute_motif_id == 80 || absolute_motif_id == 96) // 1001000 1010000 1100000 : 1 1 0
			return {6, 20};
		else if (absolute_motif_id == 8 || absolute_motif_id == 16 || absolute_motif_id == 32) // 0001000 0010000 0100000 : 0 1 0
			return {6, 21};
		else if (absolute_motif_id == 64) // 1000000 : 1 0 0
			return {6, 22};
		else if (absolute_motif_id == 0) // 0000000 : 0 0 0
			return {6, 23};	
	}
	else if (id_to_index[motif_id] - 4 == 7) { // h-motif 7
		if (absolute_motif_id == 88 || absolute_motif_id == 104 || absolute_motif_id == 112) // 1011000 1101000 1110000 : 1 2 0
			return {7, 0};
		else if (absolute_motif_id == 24 || absolute_motif_id == 40 || absolute_motif_id == 48) // 0011000 0101000 0110000 : 0 2 0
			return {7, 1};
		else if (absolute_motif_id == 72 || absolute_motif_id == 80 || absolute_motif_id == 96) // 1001000 1010000 1100000 : 1 1 0
			return {7, 2};
		else if (absolute_motif_id == 8 || absolute_motif_id == 16 || absolute_motif_id == 32) // 0001000 0010000 0100000 : 0 1 0
			return {7, 3};
		else if (absolute_motif_id == 64) // 1000000 : 1 0 0
			return {7, 4};
		else if (absolute_motif_id == 0) // 0000000 : 0 0 0
			return {7, 5};
	}
	else if (id_to_index[motif_id] - 4 == 8) { // h-motif 8
		if (absolute_motif_id == 89 || absolute_motif_id == 105 || absolute_motif_id == 113 || 
			absolute_motif_id == 90 || absolute_motif_id == 106 || absolute_motif_id == 114 || 
			absolute_motif_id == 92 || absolute_motif_id == 108 || absolute_motif_id == 116) // 1011001 1101001 1110001 1011010 1101010 1110010 1011100 1101100 1110100 : 1 2 1
			return {8, 0};
		else if (absolute_motif_id == 25 || absolute_motif_id == 41 || absolute_motif_id == 49 || 
				 absolute_motif_id == 26 || absolute_motif_id == 42 || absolute_motif_id == 50 || 
				 absolute_motif_id == 28 || absolute_motif_id == 44 || absolute_motif_id == 52) // 0011001 0101001 0110001 0011010 0101010 0110010 0011100 0101100 0110100 : 0 2 1
			return {8, 1};
		else if (absolute_motif_id == 73 || absolute_motif_id == 74 || absolute_motif_id == 76 || 
				 absolute_motif_id == 81 || absolute_motif_id == 82 || absolute_motif_id == 84 || 
				 absolute_motif_id == 97 || absolute_motif_id == 98 || absolute_motif_id == 100) // 1001001 1001010 1001100 1010001 1010010 1010100 1100001 1100010 1100100 : 1 1 1
			return {8, 2};
		else if (absolute_motif_id == 9 || absolute_motif_id == 10 || absolute_motif_id == 12 || 
				 absolute_motif_id == 17 || absolute_motif_id == 18 || absolute_motif_id == 20 || 
				 absolute_motif_id == 33 || absolute_motif_id == 34 || absolute_motif_id == 36) // 0001001 0001010 0001100 0010001 0010010 0010100 0100001 0100010 0100100 : 0 1 1
			return {8, 3};
		else if (absolute_motif_id == 65 || absolute_motif_id == 66 || absolute_motif_id == 68) // 1000001 1000010 1000100 : 1 0 1
			return {8, 4};
		else if (absolute_motif_id == 1 || absolute_motif_id == 2 || absolute_motif_id == 4) // 0000001 0000010 0000100 : 0 0 1
			return {8, 5};
		else if (absolute_motif_id == 88 || absolute_motif_id == 104 || absolute_motif_id == 112) // 1011000 1101000 1110000 : 1 2 0
			return {8, 6};
		else if (absolute_motif_id == 24 || absolute_motif_id == 40 || absolute_motif_id == 48) // 0011000 0101000 0110000 : 0 2 0
			return {8, 7};
		else if (absolute_motif_id == 72 || absolute_motif_id == 80 || absolute_motif_id == 96) // 1001000 1010000 1100000 : 1 1 0
			return {8, 8};
		else if (absolute_motif_id == 8 || absolute_motif_id == 16 || absolute_motif_id == 32) // 0001000 0010000 0100000 : 0 1 0
			return {8, 9};
		else if (absolute_motif_id == 64) // 1000000 : 1 0 0
			return {8, 10};
		else if (absolute_motif_id == 0) // 0000000 : 0 0 0
			return {8, 11};
	}
	else if (id_to_index[motif_id] - 4 == 9) { // h-motif 9 : not symmetric
		if (absolute_motif_id == 89 || absolute_motif_id == 105 || absolute_motif_id == 113 || 
			absolute_motif_id == 90 || absolute_motif_id == 106 || absolute_motif_id == 114 || 
			absolute_motif_id == 92 || absolute_motif_id == 108 || absolute_motif_id == 116) // 1011001 1101001 1110001 1011010 1101010 1110010 1011100 1101100 1110100 : 1 2 1
			return {9, 0};
		else if (absolute_motif_id == 25 || absolute_motif_id == 41 || absolute_motif_id == 49 || 
				 absolute_motif_id == 26 || absolute_motif_id == 42 || absolute_motif_id == 50 || 
				 absolute_motif_id == 28 || absolute_motif_id == 44 || absolute_motif_id == 52) // 0011001 0101001 0110001 0011010 0101010 0110010 0011100 0101100 0110100 : 0 2 1
			return {9, 1};
		else if (absolute_motif_id == 73 || absolute_motif_id == 74 || absolute_motif_id == 82 || 
				 absolute_motif_id == 84 || absolute_motif_id == 97 || absolute_motif_id == 100) // 1001001 1001010 1010010 1010100 1100001 1100100
			return {9, 2};
		else if (absolute_motif_id == 76 || absolute_motif_id == 81 || absolute_motif_id == 98) // 1001100 1010001 1100010
			return {9, 3};
		else if (absolute_motif_id == 9 || absolute_motif_id == 10 || absolute_motif_id == 18 || 
				 absolute_motif_id == 20 || absolute_motif_id == 33 || absolute_motif_id == 36) // 0001001 0001010 0010010 0010100 0100001 0100100
			return {9, 4};
		else if (absolute_motif_id == 12 || absolute_motif_id == 17 || absolute_motif_id == 34) // 0001100 0010001 0100010
			return {9, 5};
		else if (absolute_motif_id == 65 || absolute_motif_id == 66 || absolute_motif_id == 68) // 1000001 1000010 1000100 : 1 0 1
			return {9, 6};
		else if (absolute_motif_id == 1 || absolute_motif_id == 2 || absolute_motif_id == 4) // 0000001 0000010 0000100 : 0 0 1
			return {9, 7};
		else if (absolute_motif_id == 88 || absolute_motif_id == 104 || absolute_motif_id == 112) // 1011000 1101000 1110000 : 1 2 0
			return {9, 8};
		else if (absolute_motif_id == 24 || absolute_motif_id == 40 || absolute_motif_id == 48) // 0011000 0101000 0110000 : 0 2 0
			return {9, 9};
		else if (motif_id - absolute_motif_id == 12 || motif_id - absolute_motif_id == 17 || motif_id - absolute_motif_id == 34) // (motif_id - absolute_motif_id) 0001100 0010001 0100010
			return {9, 10};
		else if (motif_id - absolute_motif_id == 9 || motif_id - absolute_motif_id == 10 || motif_id - absolute_motif_id == 18 || 
				 motif_id - absolute_motif_id == 20 || motif_id - absolute_motif_id == 33 || motif_id - absolute_motif_id == 36) // (motif_id - absolute_motif_id) 0001001 0001010 0010010 0010100 0100001 0100100
			return {9, 11};
		else if (motif_id - absolute_motif_id == 76 || motif_id - absolute_motif_id == 81 || motif_id - absolute_motif_id == 98) // (motif_id - absolute_motif_id) 1001100 1010001 1100010
			return {9, 12};
		else if (motif_id - absolute_motif_id == 73 || motif_id - absolute_motif_id == 74 || motif_id - absolute_motif_id == 82 || 
				 motif_id - absolute_motif_id == 84 || motif_id - absolute_motif_id == 97 || motif_id - absolute_motif_id == 100) // (motif_id - absolute_motif_id) 1001001 1001010 1010010 1010100 1100001 1100100
			return {9, 13};
		else if (absolute_motif_id == 64) // 1000000 : 1 0 0
			return {9, 14};
		else if (absolute_motif_id == 0) // 0000000 : 0 0 0
			return {9, 15};
	}
	else if (id_to_index[motif_id] - 4 == 10) { // h-motif 10 : not symmetric
		if (absolute_motif_id == 91 || absolute_motif_id == 93 || absolute_motif_id == 94 || 
			absolute_motif_id == 107 || absolute_motif_id == 109 || absolute_motif_id == 110 || 
			absolute_motif_id == 115 || absolute_motif_id == 117 || absolute_motif_id == 118) // 1011011 1011101 1011110 1101011 1101101 1101110 1110011 1110101 1110110 : 1 2 2
			return {10, 0};
		else if (absolute_motif_id == 27 || absolute_motif_id == 29 || absolute_motif_id == 30 || 
				 absolute_motif_id == 43 || absolute_motif_id == 45 || absolute_motif_id == 46 || 
				 absolute_motif_id == 51 || absolute_motif_id == 53 || absolute_motif_id == 54) // 0011011 0011101 0011110 0101011 0101101 0101110 0110011 0110101 0110110 : 0 2 2
			return {10, 1};
		else if (absolute_motif_id == 75 || absolute_motif_id == 86 || absolute_motif_id == 101) // 1001011 1010110 1100101
			return {10, 2};
		else if (absolute_motif_id == 77 || absolute_motif_id == 78 || absolute_motif_id == 83 || 
				 absolute_motif_id == 85 || absolute_motif_id == 99 || absolute_motif_id == 102) // 1001101 1001110 1010011 1010101 1100011 1100110
			return {10, 3};
		else if (absolute_motif_id == 11 || absolute_motif_id == 22 || absolute_motif_id == 37) // 0001011 0010110 0100101
			return {10, 4};
		else if (absolute_motif_id == 13 || absolute_motif_id == 14 || absolute_motif_id == 19 || 
				 absolute_motif_id == 21 || absolute_motif_id == 35 || absolute_motif_id == 38) // 0001101 0001110 0010011 0010101 0100011 0100110
			return {10, 5};
		else if (absolute_motif_id == 67 || absolute_motif_id == 69 || absolute_motif_id == 70) // 1000011 1000101 1000110 : 1 0 2
			return {10, 6};
		else if (absolute_motif_id == 3 || absolute_motif_id == 5 || absolute_motif_id == 6) // 0000011 0000101 0000110 : 0 0 2
			return {10, 7};
		else if (absolute_motif_id == 90 || absolute_motif_id == 105 || absolute_motif_id == 116) // 1011010 1101001 1110100
			return {10, 8};
		else if (absolute_motif_id == 89 || absolute_motif_id == 92 || absolute_motif_id == 106 || 
				 absolute_motif_id == 108 || absolute_motif_id == 113 || absolute_motif_id == 114) // 1011001 1011100 1101010 1101100 1110001 1110010  
			return {10, 9};
		else if (absolute_motif_id == 26 || absolute_motif_id == 41 || absolute_motif_id == 52) // 0011010 0101001 0110100
			return {10, 10};
		else if (absolute_motif_id == 25 || absolute_motif_id == 28 || absolute_motif_id == 42 || 
				 absolute_motif_id == 44 || absolute_motif_id == 49 || absolute_motif_id == 50) // 0011001 0011100 0101010 0101100 0110001 0110010  
			return {10, 11};
		else if (absolute_motif_id == 73) { // 1001001
			if (motif_id - absolute_motif_id == 18) // 0010010
				return {10, 12};
			else if (motif_id - absolute_motif_id == 36) // 0100100
				return {10, 13};
			else if (motif_id - absolute_motif_id == 34) // 0100010
				return {10, 14};
		} 
		else if (absolute_motif_id == 74) { // 1001010
			if (motif_id - absolute_motif_id == 33) // 0100001
				return {10, 12};
			else if (motif_id - absolute_motif_id == 20) // 0010100
				return {10, 13};
			else if (motif_id - absolute_motif_id == 17) // 0010001
				return {10, 14};
		} 
		else if (absolute_motif_id == 82) { // 1010010 
			if (motif_id - absolute_motif_id == 36) // 0100100
				return {10, 12};
			else if (motif_id - absolute_motif_id == 9) // 0001001
				return {10, 13};
			else if (motif_id - absolute_motif_id == 12) // 0001100
				return {10, 14};
		}
		else if (absolute_motif_id == 84) { // 1010100
			if (motif_id - absolute_motif_id == 10) // 0001010
				return {10, 12};
			else if (motif_id - absolute_motif_id == 33) // 0100001
				return {10, 13};
			else if (motif_id - absolute_motif_id == 34) // 0100010
				return {10, 14};
		}
		else if (absolute_motif_id == 97) { // 1100001
			if (motif_id - absolute_motif_id == 20) // 0010100
				return {10, 12};
			else if (motif_id - absolute_motif_id == 10) // 0001010
				return {10, 13};
			else if (motif_id - absolute_motif_id == 12) // 0001100
				return {10, 14};
		} 
		else if (absolute_motif_id == 100) { // 1100100
			if (motif_id - absolute_motif_id == 9) // 0001001
				return {10, 12};
			else if (motif_id - absolute_motif_id == 18) // 0010010
				return {10, 13};
			else if (motif_id - absolute_motif_id == 17) // 0010001
				return {10, 14};
		}      
		else if (absolute_motif_id == 76 || absolute_motif_id == 81 || absolute_motif_id == 98) // 1001100 1010001 1100010
			return {10, 15};
		else if (absolute_motif_id == 9) { // 0001001
			if (motif_id - absolute_motif_id == 82) // 1010010
				return {10, 16};
			else if (motif_id - absolute_motif_id == 100) // 1100100
				return {10, 17};
			else if (motif_id - absolute_motif_id == 98) // 1100010
				return {10, 18};
		} 
		else if (absolute_motif_id == 10) { // 0001010
			if (motif_id - absolute_motif_id == 97) // 1100001
				return {10, 16};
			else if (motif_id - absolute_motif_id == 84) // 1010100
				return {10, 17};
			else if (motif_id - absolute_motif_id == 81) // 1010001
				return {10, 18};
		} 
		else if (absolute_motif_id == 18) { // 0010010 
			if (motif_id - absolute_motif_id == 100) // 1100100
				return {10, 16};
			else if (motif_id - absolute_motif_id == 73) // 1001001
				return {10, 17};
			else if (motif_id - absolute_motif_id == 76) // 1001100
				return {10, 18};
		}
		else if (absolute_motif_id == 20) { // 0010100
			if (motif_id - absolute_motif_id == 74) // 1001010
				return {10, 16};
			else if (motif_id - absolute_motif_id == 97) // 1100001
				return {10, 17};
			else if (motif_id - absolute_motif_id == 98) // 1100010
				return {10, 18};
		}
		else if (absolute_motif_id == 33) { // 0100001
			if (motif_id - absolute_motif_id == 84) // 1010100
				return {10, 16};
			else if (motif_id - absolute_motif_id == 74) // 1001010
				return {10, 17};
			else if (motif_id - absolute_motif_id == 76) // 1001100
				return {10, 18};
		} 
		else if (absolute_motif_id == 36) { // 0100100
			if (motif_id - absolute_motif_id == 73) // 1001001
				return {10, 16};
			else if (motif_id - absolute_motif_id == 82) // 1010010
				return {10, 17};
			else if (motif_id - absolute_motif_id == 81) // 1010001
				return {10, 18};
		}      
		else if (absolute_motif_id == 12 || absolute_motif_id == 17 || absolute_motif_id == 34) // 0001100 0010001 0100010
			return {10, 19};
		else if (motif_id - absolute_motif_id == 25 || motif_id - absolute_motif_id == 28 || motif_id - absolute_motif_id == 42 || 
				 motif_id - absolute_motif_id == 44 || motif_id - absolute_motif_id == 49 || motif_id - absolute_motif_id == 50) // (motif_id - absolute_motif_id) 0011001 0011100 0101010 0101100 0110001 0110010  
			return {10, 20};
		else if (motif_id - absolute_motif_id == 26 || motif_id - absolute_motif_id == 41 || motif_id - absolute_motif_id == 52) // (motif_id - absolute_motif_id) 0011010 0101001 0110100
			return {10, 21};
		else if (motif_id - absolute_motif_id == 89 || motif_id - absolute_motif_id == 92 || motif_id - absolute_motif_id == 106 || 
				 motif_id - absolute_motif_id == 108 || motif_id - absolute_motif_id == 113 || motif_id - absolute_motif_id == 114) // (motif_id - absolute_motif_id) 1011001 1011100 1101010 1101100 1110001 1110010  
			return {10, 22};
		else if (motif_id - absolute_motif_id == 90 || motif_id - absolute_motif_id == 105 || motif_id - absolute_motif_id == 116) // (motif_id - absolute_motif_id) 1011010 1101001 1110100
			return {10, 23};
		else if (absolute_motif_id == 88 || absolute_motif_id == 104 || absolute_motif_id == 112) // 1011000 1101000 1110000 : 1 2 0
			return {10, 24};
		else if (absolute_motif_id == 24 || absolute_motif_id == 40 || absolute_motif_id == 48) // 0011000 0101000 0110000 : 0 2 0
			return {10, 25};
		else if (motif_id - absolute_motif_id == 13 || motif_id - absolute_motif_id == 14 || motif_id - absolute_motif_id == 19 || 
				 motif_id - absolute_motif_id == 21 || motif_id - absolute_motif_id == 35 || motif_id - absolute_motif_id == 38) // (motif_id - absolute_motif_id) 0001101 0001110 0010011 0010101 0100011 0100110
			return {10, 26};
		else if (motif_id - absolute_motif_id == 11 || motif_id - absolute_motif_id == 22 || motif_id - absolute_motif_id == 37) // (motif_id - absolute_motif_id) 0001011 0010110 0100101
			return {10, 27};
		else if (motif_id - absolute_motif_id == 77 || motif_id - absolute_motif_id == 78 || motif_id - absolute_motif_id == 83 || 
				 motif_id - absolute_motif_id == 85 || motif_id - absolute_motif_id == 99 || motif_id - absolute_motif_id == 102) // (motif_id - absolute_motif_id) 1001101 1001110 1010011 1010101 1100011 1100110
			return {10, 28};
		else if (motif_id - absolute_motif_id == 75 || motif_id - absolute_motif_id == 86 || motif_id - absolute_motif_id == 101) // (motif_id - absolute_motif_id) 1001011 1010110 1100101
			return {10, 29};
		else if (absolute_motif_id == 64) // 1000000 : 1 0 0
			return {10, 30};
		else if (absolute_motif_id == 0) // 0000000 : 0 0 0
			return {10, 31};
	}
	else if (id_to_index[motif_id] - 4 == 11) { // h-motif 11 : not symmetric
		if (absolute_motif_id == 91 || absolute_motif_id == 93 || absolute_motif_id == 94 || 
			absolute_motif_id == 107 || absolute_motif_id == 109 || absolute_motif_id == 110 || 
			absolute_motif_id == 115 || absolute_motif_id == 117 || absolute_motif_id == 118) // 1011011 1011101 1011110 1101011 1101101 1101110 1110011 1110101 1110110 : 1 2 2
			return {11, 0};
		else if (absolute_motif_id == 27 || absolute_motif_id == 29 || absolute_motif_id == 30 || 
				 absolute_motif_id == 43 || absolute_motif_id == 45 || absolute_motif_id == 46 || 
				 absolute_motif_id == 51 || absolute_motif_id == 53 || absolute_motif_id == 54) // 0011011 0011101 0011110 0101011 0101101 0101110 0110011 0110101 0110110 : 0 2 2
			return {11, 1};
		else if (absolute_motif_id == 75 || absolute_motif_id == 77 || absolute_motif_id == 78 || 
				 absolute_motif_id == 83 || absolute_motif_id == 85 || absolute_motif_id == 86 || 
				 absolute_motif_id == 99 || absolute_motif_id == 101 || absolute_motif_id == 102) // 1001011 1001101 1001110 1010011 1010101 1010110 1100011 1100101 1100110 : 1 1 2
			return {11, 2};
		else if (absolute_motif_id == 11 || absolute_motif_id == 13 || absolute_motif_id == 14 || 
				 absolute_motif_id == 19 || absolute_motif_id == 21 || absolute_motif_id == 22 || 
				 absolute_motif_id == 35 || absolute_motif_id == 37 || absolute_motif_id == 38) // 0001011 0001101 0001110 0010011 0010101 0010110 0100011 0100101 0100110 : 0 1 2
			return {11, 3};
		else if (absolute_motif_id == 67 || absolute_motif_id == 69 || absolute_motif_id == 70) // 1000011 1000101 1000110 : 1 0 2
			return {11, 4};
		else if (absolute_motif_id == 3 || absolute_motif_id == 5 || absolute_motif_id == 6) // 0000011 0000101 0000110 : 0 0 2
			return {11, 5};
		else if (absolute_motif_id == 89 || absolute_motif_id == 105 || absolute_motif_id == 113 || 
				 absolute_motif_id == 90 || absolute_motif_id == 106 || absolute_motif_id == 114 || 
				 absolute_motif_id == 92 || absolute_motif_id == 108 || absolute_motif_id == 116) // 1011001 1101001 1110001 1011010 1101010 1110010 1011100 1101100 1110100 : 1 2 1
			return {11, 6};
		else if (absolute_motif_id == 25 || absolute_motif_id == 41 || absolute_motif_id == 49 || 
				 absolute_motif_id == 26 || absolute_motif_id == 42 || absolute_motif_id == 50 || 
				 absolute_motif_id == 28 || absolute_motif_id == 44 || absolute_motif_id == 52) // 0011001 0101001 0110001 0011010 0101010 0110010 0011100 0101100 0110100 : 0 2 1
			return {11, 7};
		else if (absolute_motif_id == 76 || absolute_motif_id == 81 || absolute_motif_id == 98) // 1001100 1010001 1100010 
			return {11, 8};
		else if (absolute_motif_id == 73 || absolute_motif_id == 74 || absolute_motif_id == 82 || 
				 absolute_motif_id == 84 || absolute_motif_id == 97 || absolute_motif_id == 100) // 1001001 1001010 1010010 1010100 1100001 1100100
			return {11, 9};
		else if (absolute_motif_id == 12 || absolute_motif_id == 17 || absolute_motif_id == 34) // 0001100 0010001 0100010 
			return {11, 10};
		else if (absolute_motif_id == 9 || absolute_motif_id == 10 || absolute_motif_id == 18 || 
				 absolute_motif_id == 20 || absolute_motif_id == 33 || absolute_motif_id == 36) // 0001001 0001010 0010010 0010100 0100001 0100100
			return {11, 11};
		else if (absolute_motif_id == 65 || absolute_motif_id == 66 || absolute_motif_id == 68) // 1000001 1000010 1000100 : 1 0 1
			return {11, 12};
		else if (absolute_motif_id == 1 || absolute_motif_id == 2 || absolute_motif_id == 4) // 0000001 0000010 0000100 : 0 0 1
			return {11, 13};
		else if (absolute_motif_id == 88 || absolute_motif_id == 104 || absolute_motif_id == 112) // 1011000 1101000 1110000 : 1 2 0
			return {11, 14};
		else if (absolute_motif_id == 24 || absolute_motif_id == 40 || absolute_motif_id == 48) // 0011000 0101000 0110000 : 0 2 0
			return {11, 15};
		else if (absolute_motif_id == 72 || absolute_motif_id == 80 || absolute_motif_id == 96) // 1001000 1010000 1100000 : 1 1 0
			return {11, 16};
		else if (absolute_motif_id == 8 || absolute_motif_id == 16 || absolute_motif_id == 32) // 0001000 0010000 0100000 : 0 1 0
			return {11, 17};
		else if (absolute_motif_id == 64) // 1000000 : 1 0 0
			return {11, 18};
		else if (absolute_motif_id == 0) // 0000000 : 0 0 0
			return {11, 19};
	}
	else if (id_to_index[motif_id] - 4 == 12) { // h-motif 12 : not symmetric
		if (absolute_motif_id == 95 || absolute_motif_id == 111 || absolute_motif_id == 119) // 1011111 1101111 1110111 : 1 2 3
			return {12, 0};
		else if (absolute_motif_id == 31 || absolute_motif_id == 47 || absolute_motif_id == 55) // 0011111 0101111 0110111 : 0 2 3
			return {12, 1};
		else if (absolute_motif_id == 79 || absolute_motif_id == 87 || absolute_motif_id == 103) // 1001111 1010111 1100111 : 1 1 3
			return {12, 2};
		else if (absolute_motif_id == 15 || absolute_motif_id == 23 || absolute_motif_id == 39) // 0001111 0010111 0100111 : 0 1 3
			return {12, 3};
		else if (absolute_motif_id == 71) // 1000111 : 1 0 3
			return {12, 4};
		else if (absolute_motif_id == 7) // 0000111 : 0 0 3
			return {12, 5};
		else if (absolute_motif_id == 93 || absolute_motif_id == 110 || absolute_motif_id == 115) // 1011101 1101110 1110011
			return {12, 6};
		else if (absolute_motif_id == 91 || absolute_motif_id == 94 || absolute_motif_id == 107 || 
				 absolute_motif_id == 109 || absolute_motif_id == 117 || absolute_motif_id == 118) // 1011011 1011110 1101011 1101101 1110101 1110110
			return {12, 7};
		else if (absolute_motif_id == 29 || absolute_motif_id == 46 || absolute_motif_id == 51) // 0011101 0101110 0110011
			return {12, 8};
		else if (absolute_motif_id == 27 || absolute_motif_id == 30 || absolute_motif_id == 43 || 
				 absolute_motif_id == 45 || absolute_motif_id == 53 || absolute_motif_id == 54) // 0011011 0011110 0101011 0101101 0110101 0110110
			return {12, 9};
		else if (absolute_motif_id == 75 || absolute_motif_id == 86 || absolute_motif_id == 101) // 1001011 1010110 1100101
			return {12, 10};
		else if (absolute_motif_id == 77 || absolute_motif_id == 78 || absolute_motif_id == 83 || 
				 absolute_motif_id == 85 || absolute_motif_id == 99 || absolute_motif_id == 102) { // 1001101 1001110 1010011 1010101 1100011 1100110
			if (motif_id - absolute_motif_id == 12 || motif_id - absolute_motif_id == 17 || motif_id - absolute_motif_id == 34)
				return {12, 11};
			else
				return {12, 12};	
		}
		else if (absolute_motif_id == 11 || absolute_motif_id == 22 || absolute_motif_id == 37) // 0001011 0010110 0100101
			return {12, 13};
		else if (absolute_motif_id == 13 || absolute_motif_id == 14 || absolute_motif_id == 19 || 
				 absolute_motif_id == 21 || absolute_motif_id == 35 || absolute_motif_id == 38) { // 0001101 0001110 0010011 0010101 0100011 0100110
			if (motif_id - absolute_motif_id == 76 || motif_id - absolute_motif_id == 81 || motif_id - absolute_motif_id == 98)
				return {12, 14};
			else
				return {12, 15};	
		}
		else if (motif_id - absolute_motif_id == 26 || motif_id - absolute_motif_id == 41 || motif_id - absolute_motif_id == 52) // (motif_id - absolute_motif_id) 0011010 0101001 0110100
			return {12, 16};
		else if (motif_id - absolute_motif_id == 25 || motif_id - absolute_motif_id == 28 || motif_id - absolute_motif_id == 42 || 
				 motif_id - absolute_motif_id == 44 || motif_id - absolute_motif_id == 49 || motif_id - absolute_motif_id == 50) // (motif_id - absolute_motif_id) 0011001 0011100 0101010 0101100 0110001 0110010  
			return {12, 17};
		else if (motif_id - absolute_motif_id == 90 || motif_id - absolute_motif_id == 105 || motif_id - absolute_motif_id == 116) // (motif_id - absolute_motif_id) 1011010 1101001 1110100
			return {12, 18};
		else if (motif_id - absolute_motif_id == 89 || motif_id - absolute_motif_id == 92 || motif_id - absolute_motif_id == 106 || 
				 motif_id - absolute_motif_id == 108 || motif_id - absolute_motif_id == 113 || motif_id - absolute_motif_id == 114) // (motif_id - absolute_motif_id) 1011001 1011100 1101010 1101100 1110001 1110010  
			return {12, 19};
		else if (absolute_motif_id == 90 || absolute_motif_id == 105 || absolute_motif_id == 116) // 1011010 1101001 1110100
			return {12, 20};
		else if (absolute_motif_id == 89 || absolute_motif_id == 92 || absolute_motif_id == 106 || 
				 absolute_motif_id == 108 || absolute_motif_id == 113 || absolute_motif_id == 114) // 1011001 1011100 1101010 1101100 1110001 1110010  
			return {12, 21};
		else if (absolute_motif_id == 26 || absolute_motif_id == 41 || absolute_motif_id == 52) // 0011010 0101001 0110100
			return {12, 22};
		else if (absolute_motif_id == 25 || absolute_motif_id == 28 || absolute_motif_id == 42 || 
				 absolute_motif_id == 44 || absolute_motif_id == 49 || absolute_motif_id == 50) // 0011001 0011100 0101010 0101100 0110001 0110010  
			return {12, 23};
		else if (motif_id - absolute_motif_id == 11 || motif_id - absolute_motif_id == 22 || motif_id - absolute_motif_id == 37) // (motif_id - absolute_motif_id) 0001011 0010110 0100101
			return {12, 24};
		else if (motif_id - absolute_motif_id == 13 || motif_id - absolute_motif_id == 14 || motif_id - absolute_motif_id == 19 || 
				 motif_id - absolute_motif_id == 21 || motif_id - absolute_motif_id == 35 || motif_id - absolute_motif_id == 38) { // (motif_id - absolute_motif_id) 0001101 0001110 0010011 0010101 0100011 0100110
			if (absolute_motif_id == 76 || absolute_motif_id == 81 || absolute_motif_id == 98)
				return {12, 25};
			else 
				return {12, 26};	
		}
		else if (motif_id - absolute_motif_id == 75 || motif_id - absolute_motif_id == 86 || motif_id - absolute_motif_id == 101) // (motif_id - absolute_motif_id) 1001011 1010110 1100101
			return {12, 27};
		else if (motif_id - absolute_motif_id == 77 || motif_id - absolute_motif_id == 78 || motif_id - absolute_motif_id == 83 || 
				 motif_id - absolute_motif_id == 85 || motif_id - absolute_motif_id == 99 || motif_id - absolute_motif_id == 102) {// (motif_id - absolute_motif_id) 1001101 1001110 1010011 1010101 1100011 1100110
			if (absolute_motif_id == 12 || absolute_motif_id == 17 || absolute_motif_id == 34)
				return {12, 28};
			else 
				return {12, 29};
		}
		else if (motif_id - absolute_motif_id == 29 || motif_id - absolute_motif_id == 46 || motif_id - absolute_motif_id == 51) // (motif_id - absolute_motif_id) 0011101 0101110 0110011
			return {12, 30};
		else if (motif_id - absolute_motif_id == 27 || motif_id - absolute_motif_id == 30 || motif_id - absolute_motif_id == 43 || 
				 motif_id - absolute_motif_id == 45 || motif_id - absolute_motif_id == 53 || motif_id - absolute_motif_id == 54) // (motif_id - absolute_motif_id) 0011011 0011110 0101011 0101101 0110101 0110110
			return {12, 31};
		else if (motif_id - absolute_motif_id == 93 || motif_id - absolute_motif_id == 110 || motif_id - absolute_motif_id == 115) // (motif_id - absolute_motif_id) 1011101 1101110 1110011
			return {12, 32};
		else if (motif_id - absolute_motif_id == 91 || motif_id - absolute_motif_id == 94 || motif_id - absolute_motif_id == 107 || 
				 motif_id - absolute_motif_id == 109 || motif_id - absolute_motif_id == 117 || motif_id - absolute_motif_id == 118) // (motif_id - absolute_motif_id) 1011011 1011110 1101011 1101101 1110101 1110110
			return {12, 33};
		else if (absolute_motif_id == 88 || absolute_motif_id == 104 || absolute_motif_id == 112) // 1011000 1101000 1110000 : 1 2 0
			return {12, 34};
		else if (absolute_motif_id == 24 || absolute_motif_id == 40 || absolute_motif_id == 48) // 0011000 0101000 0110000 : 0 2 0
			return {12, 35};
		else if (absolute_motif_id == 72 || absolute_motif_id == 80 || absolute_motif_id == 96) // 1001000 1010000 1100000 : 1 1 0
			return {12, 36};
		else if (absolute_motif_id == 8 || absolute_motif_id == 16 || absolute_motif_id == 32) // 0001000 0010000 0100000 : 0 1 0
			return {12, 37};
		else if (absolute_motif_id == 64) // 1000000 : 1 0 0
			return {12, 38};
		else if (absolute_motif_id == 0) // 0000000 : 0 0 0
			return {12, 39};
	}
	else if (id_to_index[motif_id] - 4 == 13) { // h-motif 13
		if (absolute_motif_id == 120) // 1111000: 1 3 0
			return {13, 0};
		else if (absolute_motif_id == 56) // 0111000: 0 3 0
			return {13, 1};
		else if (absolute_motif_id == 88 || absolute_motif_id == 104 || absolute_motif_id == 112) // 1011000 1101000 1110000 : 1 2 0
			return {13, 2};
		else if (absolute_motif_id == 24 || absolute_motif_id == 40 || absolute_motif_id == 48) // 0011000 0101000 0110000 : 0 2 0
			return {13, 3};
		else if (absolute_motif_id == 72 || absolute_motif_id == 80 || absolute_motif_id == 96) // 1001000 1010000 1100000 : 1 1 0
			return {13, 4};
		else if (absolute_motif_id == 8 || absolute_motif_id == 16 || absolute_motif_id == 32) // 0001000 0010000 0100000 : 0 1 0
			return {13, 5};
		else if (absolute_motif_id == 64) // 1000000 : 1 0 0
			return {13, 6};
		else if (absolute_motif_id == 0) // : 0000000 : 0 0 0
			return {13, 7};
	}
	else if (id_to_index[motif_id] - 4 == 14) { // h-motif 14 : not symmetric
		if (absolute_motif_id == 121 || absolute_motif_id == 122 || absolute_motif_id == 124) // 1111001 1111010 1111100 : 1 3 1
			return {14, 0};
		else if (absolute_motif_id == 57 || absolute_motif_id == 58 || absolute_motif_id == 60) // 0111001 0111010 0111100 : 0 3 1
			return {14, 1};
		else if (absolute_motif_id == 90 || absolute_motif_id == 105 || absolute_motif_id == 116) // 1011010 1101001 1110100
			return {14, 2};
		else if (absolute_motif_id == 89 || absolute_motif_id == 92 || absolute_motif_id == 106 || 
				 absolute_motif_id == 108 || absolute_motif_id == 113 || absolute_motif_id == 114) // 1011001 1011100 1101010 1101100 1110001 1110010  
			return {14, 3};
		else if (absolute_motif_id == 26 || absolute_motif_id == 41 || absolute_motif_id == 52) // 0011010 0101001 0110100
			return {14, 4};
		else if (absolute_motif_id == 25 || absolute_motif_id == 28 || absolute_motif_id == 42 || 
				 absolute_motif_id == 44 || absolute_motif_id == 49 || absolute_motif_id == 50) // 0011001 0011100 0101010 0101100 0110001 0110010  
			return {14, 5};
		else if (absolute_motif_id == 76 || absolute_motif_id == 81 || absolute_motif_id == 98) // 1001100 1010001 1100010 
			return {14, 6};
		else if (absolute_motif_id == 73 || absolute_motif_id == 74 || absolute_motif_id == 82 || 
				 absolute_motif_id == 84 || absolute_motif_id == 97 || absolute_motif_id == 100) // 1001001 1001010 1010010 1010100 1100001 1100100
			return {14, 7};
		else if (absolute_motif_id == 12 || absolute_motif_id == 17 || absolute_motif_id == 34) // 0001100 0010001 0100010 
			return {14, 8};
		else if (absolute_motif_id == 9 || absolute_motif_id == 10 || absolute_motif_id == 18 || 
				 absolute_motif_id == 20 || absolute_motif_id == 33 || absolute_motif_id == 36) // 0001001 0001010 0010010 0010100 0100001 0100100
			return {14, 9};
		else if (absolute_motif_id == 65 || absolute_motif_id == 66 || absolute_motif_id == 68) // 1000001 1000010 1000100 : 1 0 1
			return {14, 10};
		else if (absolute_motif_id == 1 || absolute_motif_id == 2 || absolute_motif_id == 4) // 0000001 0000010 0000100 : 0 0 1
			return {14, 11};
		else if (absolute_motif_id == 120) // 1111000 : 1 3 0
			return {14, 12};
		else if (absolute_motif_id == 56) // 0111000 : 0 3 0
			return {14, 13};
		else if (motif_id - absolute_motif_id == 12 || motif_id - absolute_motif_id == 17 || motif_id - absolute_motif_id == 34) // (motif_id - absolute_motif_id) 0001100 0010001 0100010 
			return {14, 14};
		else if (motif_id - absolute_motif_id == 9 || motif_id - absolute_motif_id == 10 || motif_id - absolute_motif_id == 18 || 
				 motif_id - absolute_motif_id == 20 || motif_id - absolute_motif_id == 33 || motif_id - absolute_motif_id == 36) // (motif_id - absolute_motif_id) 0001001 0001010 0010010 0010100 0100001 0100100
			return {14, 15};
		else if (motif_id - absolute_motif_id == 76 || motif_id - absolute_motif_id == 81 || motif_id - absolute_motif_id == 98) // (motif_id - absolute_motif_id) 1001100 1010001 1100010 
			return {14, 16};
		else if (motif_id - absolute_motif_id == 73 || motif_id - absolute_motif_id == 74 || motif_id - absolute_motif_id == 82 || 
				 motif_id - absolute_motif_id == 84 || motif_id - absolute_motif_id == 97 || motif_id - absolute_motif_id == 100) // (motif_id - absolute_motif_id) 1001001 1001010 1010010 1010100 1100001 1100100
			return {14, 17};
		else if (motif_id - absolute_motif_id == 26 || motif_id - absolute_motif_id == 41 || motif_id - absolute_motif_id == 52) // (motif_id - absolute_motif_id) 0011010 0101001 0110100
			return {14, 18};
		else if (motif_id - absolute_motif_id == 25 || motif_id - absolute_motif_id == 28 || motif_id - absolute_motif_id == 42 || 
				 motif_id - absolute_motif_id == 44 || motif_id - absolute_motif_id == 49 || motif_id - absolute_motif_id == 50) // (motif_id - absolute_motif_id) 0011001 0011100 0101010 0101100 0110001 0110010  
			return {14, 19};
		else if (motif_id - absolute_motif_id == 90 || motif_id - absolute_motif_id == 105 || motif_id - absolute_motif_id == 116) // (motif_id - absolute_motif_id) 1011010 1101001 1110100
			return {14, 20};
		else if (motif_id - absolute_motif_id == 89 || motif_id - absolute_motif_id == 92 || motif_id - absolute_motif_id == 106 || 
				 motif_id - absolute_motif_id == 108 || motif_id - absolute_motif_id == 113 || motif_id - absolute_motif_id == 114) // (motif_id - absolute_motif_id) 1011001 1011100 1101010 1101100 1110001 1110010  
			return {14, 21};
		else if (absolute_motif_id == 64) // 1000000 : 1 0 0
			return {14, 22};
		else if (absolute_motif_id == 0) // 0000000 : 0 0 0
			return {14, 23};
	}
	else if (id_to_index[motif_id] - 4 == 15) { // h-motif 15 : not symmetric
		if (absolute_motif_id == 123 || absolute_motif_id == 125 || absolute_motif_id == 126) // 1111011 1111101 1111110 : 1 3 2
			return {15, 0};
		else if (absolute_motif_id == 59 || absolute_motif_id == 61 || absolute_motif_id == 62) // 0111011 0111101 0111110 : 0 3 2
			return {15, 1};
		else if (absolute_motif_id == 93 || absolute_motif_id == 110 || absolute_motif_id == 115) // 1011101 1101110 1110011
			return {15, 2};
		else if (absolute_motif_id == 91 || absolute_motif_id == 94 || absolute_motif_id == 107 || 
				 absolute_motif_id == 109 || absolute_motif_id == 117 || absolute_motif_id == 118) // 1011011 1011110 1101011 1101101 1110101 1110110
			return {15, 3};
		else if (absolute_motif_id == 29 || absolute_motif_id == 46 || absolute_motif_id == 51) // 0011101 0101110 0110011
			return {15, 4};
		else if (absolute_motif_id == 27 || absolute_motif_id == 30 || absolute_motif_id == 43 || 
				 absolute_motif_id == 45 || absolute_motif_id == 53 || absolute_motif_id == 54) // 0011011 0011110 0101011 0101101 0110101 0110110
			return {15, 5};
		else if (absolute_motif_id == 75 || absolute_motif_id == 86 || absolute_motif_id == 101) // 1001011 1010110 1100101
			return {15, 6};
		else if (absolute_motif_id == 77 || absolute_motif_id == 78 || absolute_motif_id == 83 || 
				 absolute_motif_id == 85 || absolute_motif_id == 99 || absolute_motif_id == 102) // 1001101 1001110 1010011 1010101 1100011 1100110
			return {15, 7};
		else if (absolute_motif_id == 11 || absolute_motif_id == 22 || absolute_motif_id == 37) // 0001011 0010110 0100101
			return {15, 8};
		else if (absolute_motif_id == 13 || absolute_motif_id == 14 || absolute_motif_id == 19 || 
				 absolute_motif_id == 21 || absolute_motif_id == 35 || absolute_motif_id == 38) // 0001101 0001110 0010011 0010101 0100011 0100110
			return {15, 9};
		else if (absolute_motif_id == 67 || absolute_motif_id == 69 || absolute_motif_id == 70) // 1000011 1000101 1000110 : 1 0 2
			return {15, 10};
		else if (absolute_motif_id == 3 || absolute_motif_id == 5 || absolute_motif_id == 6) // 0000011 0000101 0000110 : 0 0 2
			return {15, 11};
		else if (absolute_motif_id == 121 || absolute_motif_id == 122 || absolute_motif_id == 124) // 1111001 1111010 1111100 : 1 3 1
			return {15, 12};
		else if (absolute_motif_id == 57 || absolute_motif_id == 58 || absolute_motif_id == 60) // 0111001 0111010 0111100 : 0 3 1
			return {15, 13};
		else if (absolute_motif_id == 90 || absolute_motif_id == 105 || absolute_motif_id == 116) // 1011010 1101001 1110100
			return {15, 14};
		else if (absolute_motif_id == 89) { // 1011001
			if (motif_id - absolute_motif_id == 34) // 0100010 
				return {15, 15};
			else if (motif_id - absolute_motif_id == 36) // 0100100
				return {15, 16};
		} 
		else if (absolute_motif_id == 92) { // 1011100
			if (motif_id - absolute_motif_id == 34) // 0100010 
				return {15, 15};
			else if (motif_id - absolute_motif_id == 33) // 0100001
				return {15, 16};
		}
		else if (absolute_motif_id == 106) { // 1101010
			if (motif_id - absolute_motif_id == 17) // 0010001
				return {15, 15};
			else if (motif_id - absolute_motif_id == 20) // 0010100
				return {15, 16};
		}
		else if (absolute_motif_id == 108) { // 1101100
			if (motif_id - absolute_motif_id == 17) // 0010001
				return {15, 15};
			else if (motif_id - absolute_motif_id == 18) // 0010010
				return {15, 16};
		}
		else if (absolute_motif_id == 113) { // 1110001
			if (motif_id - absolute_motif_id == 12) // 0001100
				return {15, 15};
			else if (motif_id - absolute_motif_id == 10) // 0001010
				return {15, 16};
		}
		else if (absolute_motif_id == 114) { // 1110010
			if (motif_id - absolute_motif_id == 12) // 0001100 
				return {15, 15};
			else if (motif_id - absolute_motif_id == 9) // 0001001
				return {15, 16};
		}
		else if (absolute_motif_id == 26 || absolute_motif_id == 41 || absolute_motif_id == 52) // 0011010 0101001 0110100
			return {15, 17};
		else if (absolute_motif_id == 25) { // 0011001
			if (motif_id - absolute_motif_id == 98) // 1100010 
				return {15, 18};
			else if (motif_id - absolute_motif_id == 100) // 1100100
				return {15, 19};
		} 
		else if (absolute_motif_id == 28) { // 0011100
			if (motif_id - absolute_motif_id == 98) // 1100010 
				return {15, 18};
			else if (motif_id - absolute_motif_id == 97) // 1100001
				return {15, 19};
		}
		else if (absolute_motif_id == 42) { // 0101010
			if (motif_id - absolute_motif_id == 81) // 1010001
				return {15, 18};
			else if (motif_id - absolute_motif_id == 84) // 1010100
				return {15, 19};
		}
		else if (absolute_motif_id == 44) { // 0101100
			if (motif_id - absolute_motif_id == 81) // 1010001
				return {15, 18};
			else if (motif_id - absolute_motif_id == 82) // 1010010
				return {15, 19};
		}
		else if (absolute_motif_id == 49) { // 0110001
			if (motif_id - absolute_motif_id == 76) // 1001100
				return {15, 18};
			else if (motif_id - absolute_motif_id == 74) // 1001010
				return {15, 19};
		}
		else if (absolute_motif_id == 50) { // 0110010
			if (motif_id - absolute_motif_id == 76) // 1001100 
				return {15, 18};
			else if (motif_id - absolute_motif_id == 73) // 1001001
				return {15, 19};
		}
		else if (motif_id - absolute_motif_id == 26 || motif_id - absolute_motif_id == 41 || motif_id - absolute_motif_id == 52) // 0011010 0101001 0110100
			return {15, 20};
		else if (motif_id - absolute_motif_id == 25) { // 0011001
			if (absolute_motif_id == 98) // 1100010 
				return {15, 21};
			else if (absolute_motif_id == 100) // 1100100
				return {15, 22};
		} 
		else if (motif_id - absolute_motif_id == 28) { // 0011100
			if (absolute_motif_id == 98) // 1100010 
				return {15, 21};
			else if (absolute_motif_id == 97) // 1100001
				return {15, 22};
		}
		else if (motif_id - absolute_motif_id == 42) { // 0101010
			if (absolute_motif_id == 81) // 1010001
				return {15, 21};
			else if (absolute_motif_id == 84) // 1010100
				return {15, 22};
		}
		else if (motif_id - absolute_motif_id == 44) { // 0101100
			if (absolute_motif_id == 81) // 1010001
				return {15, 21};
			else if (absolute_motif_id == 82) // 1010010
				return {15, 22};
		}
		else if (motif_id - absolute_motif_id == 49) { // 0110001
			if (absolute_motif_id == 76) // 1001100
				return {15, 21};
			else if (absolute_motif_id == 74) // 1001010
				return {15, 22};
		}
		else if (motif_id - absolute_motif_id == 50) { // 0110010
			if (absolute_motif_id == 76) // 1001100 
				return {15, 21};
			else if (absolute_motif_id == 73) // 1001001
				return {15, 22};
		}
		else if (motif_id - absolute_motif_id == 90 || motif_id - absolute_motif_id == 105 || motif_id - absolute_motif_id == 116) // 1011010 1101001 1110100
			return {15, 23};
		else if (motif_id - absolute_motif_id == 89) { // 1011001
			if (absolute_motif_id == 34) // 0100010 
				return {15, 24};
			else if (absolute_motif_id == 36) // 0100100
				return {15, 25};
		} 
		else if (motif_id - absolute_motif_id == 92) { // 1011100
			if (absolute_motif_id == 34) // 0100010 
				return {15, 24};
			else if (absolute_motif_id == 33) // 0100001
				return {15, 25};
		}
		else if (motif_id - absolute_motif_id == 106) { // 1101010
			if (absolute_motif_id == 17) // 0010001
				return {15, 24};
			else if (absolute_motif_id == 20) // 0010100
				return {15, 25};
		}
		else if (motif_id - absolute_motif_id == 108) { // 1101100
			if (absolute_motif_id == 17) // 0010001
				return {15, 24};
			else if (absolute_motif_id == 18) // 0010010
				return {15, 25};
		}
		else if (motif_id - absolute_motif_id == 113) { // 1110001
			if (absolute_motif_id == 12) // 0001100
				return {15, 24};
			else if (absolute_motif_id == 10) // 0001010
				return {15, 25};
		}
		else if (motif_id - absolute_motif_id == 114) { // 1110010
			if (absolute_motif_id == 12) // 0001100 
				return {15, 24};
			else if (absolute_motif_id == 9) // 0001001
				return {15, 25};
		}
		else if (absolute_motif_id == 65 || absolute_motif_id == 66 || absolute_motif_id == 68) // 1000001 1000010 1000100 : 1 0 1
			return {15, 26};
		else if (absolute_motif_id == 1 || absolute_motif_id == 2 || absolute_motif_id == 4) // 0000001 0000010 0000100 : 0 0 1
			return {15, 27};
		else if (absolute_motif_id == 120) // 1111000 : 1 3 0
			return {15, 28};
		else if (absolute_motif_id == 56) // 0111000 : 0 3 0
			return {15, 29};
		else if (motif_id - absolute_motif_id == 11 || motif_id - absolute_motif_id == 22 || motif_id - absolute_motif_id == 37) // (motif_id - absolute_motif_id) 0001011 0010110 0100101
			return {15, 30};
		else if (motif_id - absolute_motif_id == 13 || motif_id - absolute_motif_id == 14 || motif_id - absolute_motif_id == 19 || 
				 motif_id - absolute_motif_id == 21 || motif_id - absolute_motif_id == 35 || motif_id - absolute_motif_id == 38) // (motif_id - absolute_motif_id) 0001101 0001110 0010011 0010101 0100011 0100110
			return {15, 31};
		else if (motif_id - absolute_motif_id == 75 || motif_id - absolute_motif_id == 86 || motif_id - absolute_motif_id == 101) // (motif_id - absolute_motif_id) 1001011 1010110 1100101
			return {15, 32};
		else if (motif_id - absolute_motif_id == 77 || motif_id - absolute_motif_id == 78 || motif_id - absolute_motif_id == 83 || 
				 motif_id - absolute_motif_id == 85 || motif_id - absolute_motif_id == 99 || motif_id - absolute_motif_id == 102) // (motif_id - absolute_motif_id) 1001101 1001110 1010011 1010101 1100011 1100110
			return {15, 33};
		else if (motif_id - absolute_motif_id == 29 || motif_id - absolute_motif_id == 46 || motif_id - absolute_motif_id == 51) // (motif_id - absolute_motif_id) 0011101 0101110 0110011
			return {15, 34};
		else if (motif_id - absolute_motif_id == 27 || motif_id - absolute_motif_id == 30 || motif_id - absolute_motif_id == 43 || 
				 motif_id - absolute_motif_id == 45 || motif_id - absolute_motif_id == 53 || motif_id - absolute_motif_id == 54) // (motif_id - absolute_motif_id) 0011011 0011110 0101011 0101101 0110101 0110110
			return {15, 35};
		else if (motif_id - absolute_motif_id == 93 || motif_id - absolute_motif_id == 110 || motif_id - absolute_motif_id == 115) // (motif_id - absolute_motif_id) 1011101 1101110 1110011
			return {15, 36};
		else if (motif_id - absolute_motif_id == 91 || motif_id - absolute_motif_id == 94 || motif_id - absolute_motif_id == 107 || 
				 motif_id - absolute_motif_id == 109 || motif_id - absolute_motif_id == 117 || motif_id - absolute_motif_id == 118) // (motif_id - absolute_motif_id) 1011011 1011110 1101011 1101101 1110101 1110110
			return {15, 37};
		else if (absolute_motif_id == 64) // 1000000 : 1 0 0
			return {15, 38};
		else if (absolute_motif_id == 0) // 0000000 : 0 0 0
			return {15, 39};
	}
	else if (id_to_index[motif_id] - 4 == 16) { // h-motif 16 : not symmetric
		if (absolute_motif_id == 127) // 1111111 : 1 3 3
			return {16, 0};
		else if (absolute_motif_id == 63) // 0111111 : 0 3 3
			return {16, 1};
		else if (absolute_motif_id == 95 || absolute_motif_id == 111 || absolute_motif_id == 119) // 1011111 1101111 1110111 : 1 2 3
			return {16, 2};
		else if (absolute_motif_id == 31 || absolute_motif_id == 47 || absolute_motif_id == 55) // 0011111 0101111 0110111 : 0 2 3
			return {16, 3};
		else if (absolute_motif_id == 79 || absolute_motif_id == 87 || absolute_motif_id == 103) // 1001111 1010111 1100111 : 1 1 3
			return {16, 4};
		else if (absolute_motif_id == 15 || absolute_motif_id == 23 || absolute_motif_id == 39) // 0001111 0010111 0100111 : 0 1 3
			return {16, 5};
		else if (absolute_motif_id == 71) // 1000111 : 1 0 3
			return {16, 6};
		else if (absolute_motif_id == 7) // 0000111 : 0 0 3
			return {16, 7};
		else if (absolute_motif_id == 123 || absolute_motif_id == 125 || absolute_motif_id == 126) // 1111011 1111101 1111110 : 1 3 2
			return {16, 8};
		else if (absolute_motif_id == 59 || absolute_motif_id == 61 || absolute_motif_id == 62) // 0111011 0111101 0111110 : 0 3 2
			return {16, 9};
		else if (absolute_motif_id == 93 || absolute_motif_id == 110 || absolute_motif_id == 115) // 1011101 1101110 1110011
			return {16, 10};
		else if (absolute_motif_id == 91 || absolute_motif_id == 94 || absolute_motif_id == 107 || 
				 absolute_motif_id == 109 || absolute_motif_id == 117 || absolute_motif_id == 118) // 1011011 1011110 1101011 1101101 1110101 1110110
			return {16, 11};
		else if (absolute_motif_id == 29 || absolute_motif_id == 46 || absolute_motif_id == 51) // 0011101 0101110 0110011
			return {16, 12};
		else if (absolute_motif_id == 27 || absolute_motif_id == 30 || absolute_motif_id == 43 || 
				 absolute_motif_id == 45 || absolute_motif_id == 53 || absolute_motif_id == 54) // 0011011 0011110 0101011 0101101 0110101 0110110
			return {16, 13};
		else if (absolute_motif_id == 75 || absolute_motif_id == 86 || absolute_motif_id == 101) // 1001011 1010110 1100101
			return {16, 14};
		else if (absolute_motif_id == 77 || absolute_motif_id == 78 || absolute_motif_id == 83 || 
				 absolute_motif_id == 85 || absolute_motif_id == 99 || absolute_motif_id == 102) // 1001101 1001110 1010011 1010101 1100011 1100110
			return {16, 15};
		else if (absolute_motif_id == 11 || absolute_motif_id == 22 || absolute_motif_id == 37) // 0001011 0010110 0100101
			return {16, 16};
		else if (absolute_motif_id == 13 || absolute_motif_id == 14 || absolute_motif_id == 19 || 
				 absolute_motif_id == 21 || absolute_motif_id == 35 || absolute_motif_id == 38) // 0001101 0001110 0010011 0010101 0100011 0100110
			return {16, 17};
		else if (absolute_motif_id == 67 || absolute_motif_id == 69 || absolute_motif_id == 70) // 1000011 1000101 1000110 : 1 0 2
			return {16, 18};
		else if (absolute_motif_id == 3 || absolute_motif_id == 5 || absolute_motif_id == 6) // 0000011 0000101 0000110 : 0 0 2
			return {16, 19};
		else if (absolute_motif_id == 121 || absolute_motif_id == 122 || absolute_motif_id == 124) // 1111001 1111010 1111100 : 1 3 1
			return {16, 20};
		else if (absolute_motif_id == 57 || absolute_motif_id == 58 || absolute_motif_id == 60) // 0111001 0111010 0111100 : 0 3 1
			return {16, 21};
		else if (absolute_motif_id == 90 || absolute_motif_id == 105 || absolute_motif_id == 116) // 1011010 1101001 1110100
			return {16, 22};
		else if (absolute_motif_id == 89 || absolute_motif_id == 92 || absolute_motif_id == 106 || 
				 absolute_motif_id == 108 || absolute_motif_id == 113 || absolute_motif_id == 114) // 1011001 1011100 1101010 1101100 1110001 1110010  
			return {16, 23};
		else if (absolute_motif_id == 26 || absolute_motif_id == 41 || absolute_motif_id == 52) // 0011010 0101001 0110100
			return {16, 24};
		else if (absolute_motif_id == 25 || absolute_motif_id == 28 || absolute_motif_id == 42 || 
				 absolute_motif_id == 44 || absolute_motif_id == 49 || absolute_motif_id == 50) // 0011001 0011100 0101010 0101100 0110001 0110010  
			return {16, 25};
		else if (absolute_motif_id == 76 || absolute_motif_id == 81 || absolute_motif_id == 98) // 1001100 1010001 1100010 
			return {16, 26};
		else if (absolute_motif_id == 73 || absolute_motif_id == 74 || absolute_motif_id == 82 || 
				 absolute_motif_id == 84 || absolute_motif_id == 97 || absolute_motif_id == 100) // 1001001 1001010 1010010 1010100 1100001 1100100
			return {16, 27};
		else if (absolute_motif_id == 12 || absolute_motif_id == 17 || absolute_motif_id == 34) // 0001100 0010001 0100010 
			return {16, 28};
		else if (absolute_motif_id == 9 || absolute_motif_id == 10 || absolute_motif_id == 18 || 
				 absolute_motif_id == 20 || absolute_motif_id == 33 || absolute_motif_id == 36) // 0001001 0001010 0010010 0010100 0100001 0100100
			return {16, 29};
		else if (absolute_motif_id == 65 || absolute_motif_id == 66 || absolute_motif_id == 68) // 1000001 1000010 1000100 : 1 0 1
			return {16, 30};
		else if (absolute_motif_id == 1 || absolute_motif_id == 2 || absolute_motif_id == 4) // 0000001 0000010 0000100 : 0 0 1
			return {16, 31};
		else if (absolute_motif_id == 120) // 1111000: 1 3 0
			return {16, 32};
		else if (absolute_motif_id == 56) // 0111000: 0 3 0
			return {16, 33};
		else if (absolute_motif_id == 88 || absolute_motif_id == 104 || absolute_motif_id == 112) // 1011000 1101000 1110000 : 1 2 0
			return {16, 34};
		else if (absolute_motif_id == 24 || absolute_motif_id == 40 || absolute_motif_id == 48) // 0011000 0101000 0110000 : 0 2 0
			return {16, 35};
		else if (absolute_motif_id == 72 || absolute_motif_id == 80 || absolute_motif_id == 96) // 1001000 1010000 1100000 : 1 1 0
			return {16, 36};
		else if (absolute_motif_id == 8 || absolute_motif_id == 16 || absolute_motif_id == 32) // 0001000 0010000 0100000 : 0 1 0
			return {16, 37};
		else if (absolute_motif_id == 64) // 1000000 : 1 0 0
			return {16, 38};
		else if (absolute_motif_id == 0) // : 0000000 : 0 0 0
			return {16, 39};
	}
	else if (id_to_index[motif_id] - 4 == 17) { // h-motif 17
		if (absolute_motif_id == 24 || absolute_motif_id == 40 || absolute_motif_id == 48) // 0011000 0101000 0110000 : 0 2 0
			return {17, 0};
		else if (absolute_motif_id == 8 || absolute_motif_id == 16 || absolute_motif_id == 32) // 0001000 0010000 0100000 : 0 1 0
			return {17, 1};
		else if (absolute_motif_id == 0) // : 0000000 : 0 0 0
			return {17, 2};
	}
	else if (id_to_index[motif_id] - 4 == 18) { // h-motif 18
		if (absolute_motif_id == 25 || absolute_motif_id == 41 || absolute_motif_id == 49 || 
			absolute_motif_id == 26 || absolute_motif_id == 42 || absolute_motif_id == 50 || 
			absolute_motif_id == 28 || absolute_motif_id == 44 || absolute_motif_id == 52) // 0011001 0101001 0110001 0011010 0101010 0110010 0011100 0101100 0110100 : 0 2 1
			return {18, 0};
		else if (absolute_motif_id == 9 || absolute_motif_id == 10 || absolute_motif_id == 12 || 
				 absolute_motif_id == 17 || absolute_motif_id == 18 || absolute_motif_id == 20 || 
				 absolute_motif_id == 33 || absolute_motif_id == 34 || absolute_motif_id == 36) // 0001001 0001010 0001100 0010001 0010010 0010100 0100001 0100010 0100100 : 0 1 1
			return {18, 1};
		else if (absolute_motif_id == 1 || absolute_motif_id == 2 || absolute_motif_id == 4) // 0000001 0000010 0000100 : 0 0 1
			return {18, 2};
		else if (absolute_motif_id == 24 || absolute_motif_id == 40 || absolute_motif_id == 48) // 0011000 0101000 0110000 : 0 2 0
			return {18, 3};
		else if (absolute_motif_id == 8 || absolute_motif_id == 16 || absolute_motif_id == 32) // 0001000 0010000 0100000 : 0 1 0
			return {18, 4};
		else if (absolute_motif_id == 0) // : 0000000 : 0 0 0
			return {18, 5};
	}
	else if (id_to_index[motif_id] - 4 == 19) { // h-motif 19 : not symmetric
		if (absolute_motif_id == 25 || absolute_motif_id == 41 || absolute_motif_id == 49 || 
			absolute_motif_id == 26 || absolute_motif_id == 42 || absolute_motif_id == 50 || 
			absolute_motif_id == 28 || absolute_motif_id == 44 || absolute_motif_id == 52) // 0011001 0101001 0110001 0011010 0101010 0110010 0011100 0101100 0110100 : 0 2 1
			return {19, 0};
		else if (absolute_motif_id == 9 || absolute_motif_id == 10 || absolute_motif_id == 18 || 
				 absolute_motif_id == 20 || absolute_motif_id == 33 || absolute_motif_id == 36) // 0001001 0001010 0010010 0010100 0100001 0100100
			return {19, 1};
		else if (absolute_motif_id == 12 || absolute_motif_id == 17 || absolute_motif_id == 34) // 0001100 0010001 0100010
			return {19, 2};
		else if (absolute_motif_id == 1 || absolute_motif_id == 2 || absolute_motif_id == 4) // 0000001 0000010 0000100 : 0 0 1
			return {19, 3};
		else if (absolute_motif_id == 24 || absolute_motif_id == 40 || absolute_motif_id == 48) // 0011000 0101000 0110000 : 0 2 0
			return {19, 4};
		else if (motif_id - absolute_motif_id == 12 || motif_id - absolute_motif_id == 17 || motif_id - absolute_motif_id == 34) // (motif_id - absolute_motif_id) 0001100 0010001 0100010
			return {19, 5};
		else if (motif_id - absolute_motif_id == 9 || motif_id - absolute_motif_id == 10 || motif_id - absolute_motif_id == 18 || 
				 motif_id - absolute_motif_id == 20 || motif_id - absolute_motif_id == 33 || motif_id - absolute_motif_id == 36) // (motif_id - absolute_motif_id) 0001001 0001010 0010010 0010100 0100001 0100100
			return {19, 6};
		else if (absolute_motif_id == 0) // 0000000 : 0 0 0
			return {19, 7};
	}
	else if (id_to_index[motif_id] - 4 == 20) { // h-motif 20 : not symmetric
		if (absolute_motif_id == 27 || absolute_motif_id == 29 || absolute_motif_id == 30 || 
			absolute_motif_id == 43 || absolute_motif_id == 45 || absolute_motif_id == 46 || 
			absolute_motif_id == 51 || absolute_motif_id == 53 || absolute_motif_id == 54) // 0011011 0011101 0011110 0101011 0101101 0101110 0110011 0110101 0110110 : 0 2 2
			return {20, 0};
		else if (absolute_motif_id == 11 || absolute_motif_id == 22 || absolute_motif_id == 37) // 0001011 0010110 0100101
			return {20, 1};
		else if (absolute_motif_id == 13 || absolute_motif_id == 14 || absolute_motif_id == 19 || 
				 absolute_motif_id == 21 || absolute_motif_id == 35 || absolute_motif_id == 38) // 0001101 0001110 0010011 0010101 0100011 0100110
			return {20, 2};
		else if (absolute_motif_id == 3 || absolute_motif_id == 5 || absolute_motif_id == 6) // 0000011 0000101 0000110 : 0 0 2
			return {20, 3};
		else if (absolute_motif_id == 26 || absolute_motif_id == 41 || absolute_motif_id == 52) // 0011010 0101001 0110100
			return {20, 4};
		else if (absolute_motif_id == 25 || absolute_motif_id == 28 || absolute_motif_id == 42 || 
				 absolute_motif_id == 44 || absolute_motif_id == 49 || absolute_motif_id == 50) // 0011001 0011100 0101010 0101100 0110001 0110010  
			return {20, 5};
		else if (absolute_motif_id == 9) { // 0001001
			if (motif_id - absolute_motif_id == 18) // 0010010
				return {20, 6};
			else if (motif_id - absolute_motif_id == 36) // 0100100
				return {20, 7};
			else if (motif_id - absolute_motif_id == 34) // 0100010
				return {20, 8};
		} 
		else if (absolute_motif_id == 10) { // 0001010
			if (motif_id - absolute_motif_id == 33) // 0100001
				return {20, 6};
			else if (motif_id - absolute_motif_id == 20) // 0010100
				return {20, 7};
			else if (motif_id - absolute_motif_id == 17) // 0010001
				return {20, 8};
		} 
		else if (absolute_motif_id == 18) { // 0010010 
			if (motif_id - absolute_motif_id == 36) // 0100100
				return {20, 6};
			else if (motif_id - absolute_motif_id == 9) // 0001001
				return {20, 7};
			else if (motif_id - absolute_motif_id == 12) // 0001100
				return {20, 8};
		}
		else if (absolute_motif_id == 20) { // 0010100
			if (motif_id - absolute_motif_id == 10) // 0001010
				return {20, 6};
			else if (motif_id - absolute_motif_id == 33) // 0100001
				return {20, 7};
			else if (motif_id - absolute_motif_id == 34) // 0100010
				return {20, 8};
		}
		else if (absolute_motif_id == 33) { // 0100001
			if (motif_id - absolute_motif_id == 20) // 0010100
				return {20, 6};
			else if (motif_id - absolute_motif_id == 10) // 0001010
				return {20, 7};
			else if (motif_id - absolute_motif_id == 12) // 0001100
				return {20, 8};
		} 
		else if (absolute_motif_id == 36) { // 0100100
			if (motif_id - absolute_motif_id == 9) // 0001001
				return {20, 6};
			else if (motif_id - absolute_motif_id == 18) // 0010010
				return {20, 7};
			else if (motif_id - absolute_motif_id == 17) // 0010001
				return {20, 8};
		}      
		else if (absolute_motif_id == 12 || absolute_motif_id == 17 || absolute_motif_id == 34) // 0001100 0010001 0100010
			return {20, 9};
		else if (motif_id - absolute_motif_id == 25 || motif_id - absolute_motif_id == 28 || motif_id - absolute_motif_id == 42 || 
				 motif_id - absolute_motif_id == 44 || motif_id - absolute_motif_id == 49 || motif_id - absolute_motif_id == 50) // (motif_id - absolute_motif_id) 0011001 0011100 0101010 0101100 0110001 0110010  
			return {20, 10};
		else if (motif_id - absolute_motif_id == 26 || motif_id - absolute_motif_id == 41 || motif_id - absolute_motif_id == 52) // (motif_id - absolute_motif_id) 0011010 0101001 0110100
			return {20, 11};
		else if (absolute_motif_id == 24 || absolute_motif_id == 40 || absolute_motif_id == 48) // 0011000 0101000 0110000 : 0 2 0
			return {20, 12};
		else if (motif_id - absolute_motif_id == 13 || motif_id - absolute_motif_id == 14 || motif_id - absolute_motif_id == 19 || 
				 motif_id - absolute_motif_id == 21 || motif_id - absolute_motif_id == 35 || motif_id - absolute_motif_id == 38) // (motif_id - absolute_motif_id) 0001101 0001110 0010011 0010101 0100011 0100110
			return {20, 13};
		else if (motif_id - absolute_motif_id == 11 || motif_id - absolute_motif_id == 22 || motif_id - absolute_motif_id == 37) // (motif_id - absolute_motif_id) 0001011 0010110 0100101
			return {20, 14};
		else if (absolute_motif_id == 0) // 0000000 : 0 0 0
			return {20, 15};
	}
	else if (id_to_index[motif_id] - 4 == 21) { // h-motif 21 : not symmetric
		if (absolute_motif_id == 27 || absolute_motif_id == 29 || absolute_motif_id == 30 || 
			absolute_motif_id == 43 || absolute_motif_id == 45 || absolute_motif_id == 46 || 
			absolute_motif_id == 51 || absolute_motif_id == 53 || absolute_motif_id == 54) // 0011011 0011101 0011110 0101011 0101101 0101110 0110011 0110101 0110110 : 0 2 2
			return {21, 0};
		else if (absolute_motif_id == 11 || absolute_motif_id == 13 || absolute_motif_id == 14 || 
				 absolute_motif_id == 19 || absolute_motif_id == 21 || absolute_motif_id == 22 || 
				 absolute_motif_id == 35 || absolute_motif_id == 37 || absolute_motif_id == 38) // 0001011 0001101 0001110 0010011 0010101 0010110 0100011 0100101 0100110 : 0 1 2
			return {21, 1};
		else if (absolute_motif_id == 3 || absolute_motif_id == 5 || absolute_motif_id == 6) // 0000011 0000101 0000110 : 0 0 2
			return {21, 2};
		else if (absolute_motif_id == 25 || absolute_motif_id == 41 || absolute_motif_id == 49 || 
				 absolute_motif_id == 26 || absolute_motif_id == 42 || absolute_motif_id == 50 || 
				 absolute_motif_id == 28 || absolute_motif_id == 44 || absolute_motif_id == 52) // 0011001 0101001 0110001 0011010 0101010 0110010 0011100 0101100 0110100 : 0 2 1
			return {21, 3};
		else if (absolute_motif_id == 12 || absolute_motif_id == 17 || absolute_motif_id == 34) // 0001100 0010001 0100010 
			return {21, 4};
		else if (absolute_motif_id == 9 || absolute_motif_id == 10 || absolute_motif_id == 18 || 
				 absolute_motif_id == 20 || absolute_motif_id == 33 || absolute_motif_id == 36) // 0001001 0001010 0010010 0010100 0100001 0100100
			return {21, 5};
		else if (absolute_motif_id == 1 || absolute_motif_id == 2 || absolute_motif_id == 4) // 0000001 0000010 0000100 : 0 0 1
			return {21, 6};
		else if (absolute_motif_id == 24 || absolute_motif_id == 40 || absolute_motif_id == 48) // 0011000 0101000 0110000 : 0 2 0
			return {21, 7};
		else if (absolute_motif_id == 8 || absolute_motif_id == 16 || absolute_motif_id == 32) // 0001000 0010000 0100000 : 0 1 0
			return {21, 8};
		else if (absolute_motif_id == 0) // 0000000 : 0 0 0
			return {21, 9};
	}
	else if (id_to_index[motif_id] - 4 == 22) { // h-motif 22 : not symmetric
		if (absolute_motif_id == 31 || absolute_motif_id == 47 || absolute_motif_id == 55) // 0011111 0101111 0110111 : 0 2 3
			return {22, 0};
		else if (absolute_motif_id == 15 || absolute_motif_id == 23 || absolute_motif_id == 39) // 0001111 0010111 0100111 : 0 1 3
			return {22, 1};
		else if (absolute_motif_id == 7) // 0000111 : 0 0 3
			return {22, 2};
		else if (absolute_motif_id == 29 || absolute_motif_id == 46 || absolute_motif_id == 51) // 0011101 0101110 0110011
			return {22, 3};
		else if (absolute_motif_id == 27 || absolute_motif_id == 30 || absolute_motif_id == 43 || 
				 absolute_motif_id == 45 || absolute_motif_id == 53 || absolute_motif_id == 54) // 0011011 0011110 0101011 0101101 0110101 0110110
			return {22, 4};
		else if (absolute_motif_id == 11 || absolute_motif_id == 22 || absolute_motif_id == 37) // 0001011 0010110 0100101
			return {22, 5};
		else if (absolute_motif_id == 13 || absolute_motif_id == 14 || absolute_motif_id == 19 || 
				 absolute_motif_id == 21 || absolute_motif_id == 35 || absolute_motif_id == 38) { // 0001101 0001110 0010011 0010101 0100011 0100110
			if (motif_id - absolute_motif_id == 12 || motif_id - absolute_motif_id == 17 || motif_id - absolute_motif_id == 34)
				return {22, 6};
			else
				return {22, 7};
		}
		else if (motif_id - absolute_motif_id == 26 || motif_id - absolute_motif_id == 41 || motif_id - absolute_motif_id == 52) // (motif_id - absolute_motif_id) 0011010 0101001 0110100
			return {22, 8};
		else if (motif_id - absolute_motif_id == 25 || motif_id - absolute_motif_id == 28 || motif_id - absolute_motif_id == 42 || 
				 motif_id - absolute_motif_id == 44 || motif_id - absolute_motif_id == 49 || motif_id - absolute_motif_id == 50) // (motif_id - absolute_motif_id) 0011001 0011100 0101010 0101100 0110001 0110010  
			return {22, 9};
		else if (absolute_motif_id == 26 || absolute_motif_id == 41 || absolute_motif_id == 52) // 0011010 0101001 0110100
			return {22, 10};
		else if (absolute_motif_id == 25 || absolute_motif_id == 28 || absolute_motif_id == 42 || 
				 absolute_motif_id == 44 || absolute_motif_id == 49 || absolute_motif_id == 50) // 0011001 0011100 0101010 0101100 0110001 0110010  
			return {22, 11};
		else if (motif_id - absolute_motif_id == 11 || motif_id - absolute_motif_id == 22 || motif_id - absolute_motif_id == 37) // (motif_id - absolute_motif_id) 0001011 0010110 0100101
			return {22, 12};
		else if (motif_id - absolute_motif_id == 13 || motif_id - absolute_motif_id == 14 || motif_id - absolute_motif_id == 19 || 
				 motif_id - absolute_motif_id == 21 || motif_id - absolute_motif_id == 35 || motif_id - absolute_motif_id == 38) { // (motif_id - absolute_motif_id) 0001101 0001110 0010011 0010101 0100011 0100110
			if (absolute_motif_id == 12 || absolute_motif_id == 17 || absolute_motif_id == 34)
				return {22, 13};
			else
				return {22, 14};	
		}
		else if (motif_id - absolute_motif_id == 29 || motif_id - absolute_motif_id == 46 || motif_id - absolute_motif_id == 51) // (motif_id - absolute_motif_id) 0011101 0101110 0110011
			return {22, 15};
		else if (motif_id - absolute_motif_id == 27 || motif_id - absolute_motif_id == 30 || motif_id - absolute_motif_id == 43 || 
				 motif_id - absolute_motif_id == 45 || motif_id - absolute_motif_id == 53 || motif_id - absolute_motif_id == 54) // (motif_id - absolute_motif_id) 0011011 0011110 0101011 0101101 0110101 0110110
			return {22, 16};
		else if (absolute_motif_id == 24 || absolute_motif_id == 40 || absolute_motif_id == 48) // 0011000 0101000 0110000 : 0 2 0
			return {22, 17};
		else if (absolute_motif_id == 8 || absolute_motif_id == 16 || absolute_motif_id == 32) // 0001000 0010000 0100000 : 0 1 0
			return {22, 18};
		else if (absolute_motif_id == 0) // 0000000 : 0 0 0
			return {22, 19};
	}
	else if (id_to_index[motif_id] - 4 == 23) { // h-motif 23
		if (absolute_motif_id == 56) // 0111000: 0 3 0
			return {23, 0};
		else if (absolute_motif_id == 24 || absolute_motif_id == 40 || absolute_motif_id == 48) // 0011000 0101000 0110000 : 0 2 0
			return {23, 1};
		else if (absolute_motif_id == 8 || absolute_motif_id == 16 || absolute_motif_id == 32) // 0001000 0010000 0100000 : 0 1 0
			return {23, 2};
		else if (absolute_motif_id == 0) // : 0000000 : 0 0 0
			return {23, 3};
	}
	else if (id_to_index[motif_id] - 4 == 24) { // h-motif 24 : not symmetric
		if (absolute_motif_id == 57 || absolute_motif_id == 58 || absolute_motif_id == 60) // 0111001 0111010 0111100 : 0 3 1
			return {24, 0};
		else if (absolute_motif_id == 26 || absolute_motif_id == 41 || absolute_motif_id == 52) // 0011010 0101001 0110100
			return {24, 1};
		else if (absolute_motif_id == 25 || absolute_motif_id == 28 || absolute_motif_id == 42 || 
				 absolute_motif_id == 44 || absolute_motif_id == 49 || absolute_motif_id == 50) // 0011001 0011100 0101010 0101100 0110001 0110010  
			return {24, 2};
		else if (absolute_motif_id == 12 || absolute_motif_id == 17 || absolute_motif_id == 34) // 0001100 0010001 0100010 
			return {24, 3};
		else if (absolute_motif_id == 9 || absolute_motif_id == 10 || absolute_motif_id == 18 || 
				 absolute_motif_id == 20 || absolute_motif_id == 33 || absolute_motif_id == 36) // 0001001 0001010 0010010 0010100 0100001 0100100
			return {24, 4};
		else if (absolute_motif_id == 1 || absolute_motif_id == 2 || absolute_motif_id == 4) // 0000001 0000010 0000100 : 0 0 1
			return {24, 5};
		else if (absolute_motif_id == 56) // 0111000 : 0 3 0
			return {24, 6};
		else if (motif_id - absolute_motif_id == 12 || motif_id - absolute_motif_id == 17 || motif_id - absolute_motif_id == 34) // (motif_id - absolute_motif_id) 0001100 0010001 0100010 
			return {24, 7};
		else if (motif_id - absolute_motif_id == 9 || motif_id - absolute_motif_id == 10 || motif_id - absolute_motif_id == 18 || 
				 motif_id - absolute_motif_id == 20 || motif_id - absolute_motif_id == 33 || motif_id - absolute_motif_id == 36) // (motif_id - absolute_motif_id) 0001001 0001010 0010010 0010100 0100001 0100100
			return {24, 8};
		else if (motif_id - absolute_motif_id == 26 || motif_id - absolute_motif_id == 41 || motif_id - absolute_motif_id == 52) // (motif_id - absolute_motif_id) 0011010 0101001 0110100
			return {24, 9};
		else if (motif_id - absolute_motif_id == 25 || motif_id - absolute_motif_id == 28 || motif_id - absolute_motif_id == 42 || 
				 motif_id - absolute_motif_id == 44 || motif_id - absolute_motif_id == 49 || motif_id - absolute_motif_id == 50) // (motif_id - absolute_motif_id) 0011001 0011100 0101010 0101100 0110001 0110010  
			return {24, 10};
		else if (absolute_motif_id == 0) // 0000000 : 0 0 0
			return {24, 11};
	}
	else if (id_to_index[motif_id] - 4 == 25) { // h-motif 25 : not symmetric
		if (absolute_motif_id == 59 || absolute_motif_id == 61 || absolute_motif_id == 62) // 0111011 0111101 0111110 : 0 3 2
			return {25, 0};
		else if (absolute_motif_id == 29 || absolute_motif_id == 46 || absolute_motif_id == 51) // 0011101 0101110 0110011
			return {25, 1};
		else if (absolute_motif_id == 27 || absolute_motif_id == 30 || absolute_motif_id == 43 || 
				 absolute_motif_id == 45 || absolute_motif_id == 53 || absolute_motif_id == 54) // 0011011 0011110 0101011 0101101 0110101 0110110
			return {25, 2};
		else if (absolute_motif_id == 11 || absolute_motif_id == 22 || absolute_motif_id == 37) // 0001011 0010110 0100101
			return {25, 3};
		else if (absolute_motif_id == 13 || absolute_motif_id == 14 || absolute_motif_id == 19 || 
				 absolute_motif_id == 21 || absolute_motif_id == 35 || absolute_motif_id == 38) // 0001101 0001110 0010011 0010101 0100011 0100110
			return {25, 4};
		else if (absolute_motif_id == 3 || absolute_motif_id == 5 || absolute_motif_id == 6) // 0000011 0000101 0000110 : 0 0 2
			return {25, 5};
		else if (absolute_motif_id == 57 || absolute_motif_id == 58 || absolute_motif_id == 60) // 0111001 0111010 0111100 : 0 3 1
			return {25, 6};
		else if (absolute_motif_id == 26 || absolute_motif_id == 41 || absolute_motif_id == 52) // 0011010 0101001 0110100
			return {25, 7};
		else if (absolute_motif_id == 25) { // 0011001
			if (motif_id - absolute_motif_id == 34) // 0100010 
				return {25, 8};
			else if (motif_id - absolute_motif_id == 36) // 0100100
				return {25, 9};
		} 
		else if (absolute_motif_id == 28) { // 0011100
			if (motif_id - absolute_motif_id == 34) // 0100010 
				return {25, 8};
			else if (motif_id - absolute_motif_id == 33) // 0100001
				return {25, 9};
		}
		else if (absolute_motif_id == 42) { // 0101010
			if (motif_id - absolute_motif_id == 17) // 0010001
				return {25, 8};
			else if (motif_id - absolute_motif_id == 20) // 0010100
				return {25, 9};
		}
		else if (absolute_motif_id == 44) { // 0101100
			if (motif_id - absolute_motif_id == 17) // 0010001
				return {25, 8};
			else if (motif_id - absolute_motif_id == 18) // 0010010
				return {25, 9};
		}
		else if (absolute_motif_id == 49) { // 0110001
			if (motif_id - absolute_motif_id == 12) // 0001100
				return {25, 8};
			else if (motif_id - absolute_motif_id == 10) // 0001010
				return {25, 9};
		}
		else if (absolute_motif_id == 50) { // 0110010
			if (motif_id - absolute_motif_id == 12) // 0001100 
				return {25, 8};
			else if (motif_id - absolute_motif_id == 9) // 0001001
				return {25, 9};
		}
		else if (motif_id - absolute_motif_id == 26 || motif_id - absolute_motif_id == 41 || motif_id - absolute_motif_id == 52) // 0011010 0101001 0110100
			return {25, 10};
		else if (motif_id - absolute_motif_id == 25) { // 0011001
			if (absolute_motif_id == 34) // 0100010 
				return {25, 11};
			else if (absolute_motif_id == 36) // 0100100
				return {25, 12};
		} 
		else if (motif_id - absolute_motif_id == 28) { // 0011100
			if (absolute_motif_id == 34) // 0100010 
				return {25, 11};
			else if (absolute_motif_id == 33) // 0100001
				return {25, 12};
		}
		else if (motif_id - absolute_motif_id == 42) { // 0101010
			if (absolute_motif_id == 17) // 0010001
				return {25, 11};
			else if (absolute_motif_id == 20) // 0010100
				return {25, 12};
		}
		else if (motif_id - absolute_motif_id == 44) { // 0101100
			if (absolute_motif_id == 17) // 0010001
				return {25, 11};
			else if (absolute_motif_id == 18) // 0010010
				return {25, 12};
		}
		else if (motif_id - absolute_motif_id == 49) { // 0110001
			if (absolute_motif_id == 12) // 0001100
				return {25, 11};
			else if (absolute_motif_id == 10) // 0001010
				return {25, 12};
		}
		else if (motif_id - absolute_motif_id == 50) { // 0110010
			if (absolute_motif_id == 12) // 0001100 
				return {25, 11};
			else if (absolute_motif_id == 9) // 0001001
				return {25, 12};
		}
		else if (absolute_motif_id == 1 || absolute_motif_id == 2 || absolute_motif_id == 4) // 0000001 0000010 0000100 : 0 0 1
			return {25, 13};
		else if (absolute_motif_id == 56) // 0111000 : 0 3 0
			return {25, 14};
		else if (motif_id - absolute_motif_id == 11 || motif_id - absolute_motif_id == 22 || motif_id - absolute_motif_id == 37) // (motif_id - absolute_motif_id) 0001011 0010110 0100101
			return {25, 15};
		else if (motif_id - absolute_motif_id == 13 || motif_id - absolute_motif_id == 14 || motif_id - absolute_motif_id == 19 || 
				 motif_id - absolute_motif_id == 21 || motif_id - absolute_motif_id == 35 || motif_id - absolute_motif_id == 38) // (motif_id - absolute_motif_id) 0001101 0001110 0010011 0010101 0100011 0100110
			return {25, 16};
		else if (motif_id - absolute_motif_id == 29 || motif_id - absolute_motif_id == 46 || motif_id - absolute_motif_id == 51) // (motif_id - absolute_motif_id) 0011101 0101110 0110011
			return {25, 17};
		else if (motif_id - absolute_motif_id == 27 || motif_id - absolute_motif_id == 30 || motif_id - absolute_motif_id == 43 || 
				 motif_id - absolute_motif_id == 45 || motif_id - absolute_motif_id == 53 || motif_id - absolute_motif_id == 54) // (motif_id - absolute_motif_id) 0011011 0011110 0101011 0101101 0110101 0110110
			return {25, 18};
		else if (absolute_motif_id == 0) // 0000000 : 0 0 0
			return {25, 19};
	}
	else if (id_to_index[motif_id] - 4 == 26) { // h-motif 26 : not symmetric
		if (absolute_motif_id == 63) // 0111111 : 0 3 3
			return {26, 0};
		else if (absolute_motif_id == 31 || absolute_motif_id == 47 || absolute_motif_id == 55) // 0011111 0101111 0110111 : 0 2 3
			return {26, 1};
		else if (absolute_motif_id == 15 || absolute_motif_id == 23 || absolute_motif_id == 39) // 0001111 0010111 0100111 : 0 1 3
			return {26, 2};
		else if (absolute_motif_id == 7) // 0000111 : 0 0 3
			return {26, 3};
		else if (absolute_motif_id == 59 || absolute_motif_id == 61 || absolute_motif_id == 62) // 0111011 0111101 0111110 : 0 3 2
			return {26, 4};
		else if (absolute_motif_id == 29 || absolute_motif_id == 46 || absolute_motif_id == 51) // 0011101 0101110 0110011
			return {26, 5};
		else if (absolute_motif_id == 27 || absolute_motif_id == 30 || absolute_motif_id == 43 || 
				 absolute_motif_id == 45 || absolute_motif_id == 53 || absolute_motif_id == 54) // 0011011 0011110 0101011 0101101 0110101 0110110
			return {26, 6};
		else if (absolute_motif_id == 11 || absolute_motif_id == 22 || absolute_motif_id == 37) // 0001011 0010110 0100101
			return {26, 7};
		else if (absolute_motif_id == 13 || absolute_motif_id == 14 || absolute_motif_id == 19 || 
				 absolute_motif_id == 21 || absolute_motif_id == 35 || absolute_motif_id == 38) // 0001101 0001110 0010011 0010101 0100011 0100110
			return {26, 8};
		else if (absolute_motif_id == 3 || absolute_motif_id == 5 || absolute_motif_id == 6) // 0000011 0000101 0000110 : 0 0 2
			return {26, 9};
		else if (absolute_motif_id == 57 || absolute_motif_id == 58 || absolute_motif_id == 60) // 0111001 0111010 0111100 : 0 3 1
			return {26, 10};
		else if (absolute_motif_id == 26 || absolute_motif_id == 41 || absolute_motif_id == 52) // 0011010 0101001 0110100
			return {26, 11};
		else if (absolute_motif_id == 25 || absolute_motif_id == 28 || absolute_motif_id == 42 || 
				 absolute_motif_id == 44 || absolute_motif_id == 49 || absolute_motif_id == 50) // 0011001 0011100 0101010 0101100 0110001 0110010  
			return {26, 12};
		else if (absolute_motif_id == 12 || absolute_motif_id == 17 || absolute_motif_id == 34) // 0001100 0010001 0100010 
			return {26, 13};
		else if (absolute_motif_id == 9 || absolute_motif_id == 10 || absolute_motif_id == 18 || 
				 absolute_motif_id == 20 || absolute_motif_id == 33 || absolute_motif_id == 36) // 0001001 0001010 0010010 0010100 0100001 0100100
			return {26, 14};
		else if (absolute_motif_id == 1 || absolute_motif_id == 2 || absolute_motif_id == 4) // 0000001 0000010 0000100 : 0 0 1
			return {26, 15};
		else if (absolute_motif_id == 56) // 0111000: 0 3 0
			return {26, 16};
		else if (absolute_motif_id == 24 || absolute_motif_id == 40 || absolute_motif_id == 48) // 0011000 0101000 0110000 : 0 2 0
			return {26, 17};
		else if (absolute_motif_id == 8 || absolute_motif_id == 16 || absolute_motif_id == 32) // 0001000 0010000 0100000 : 0 1 0
			return {26, 18};
		else if (absolute_motif_id == 0) // : 0000000 : 0 0 0
			return {26, 19};
	}
	else 
		return {0, 0};
}
