all:
	#g++ gism.cpp -o gism -std=c++11 

	g++ -std=c++11 -O3 -msse4.2 -fomit-frame-pointer -funroll-loops -DNDEBUG -I ~/include -L ~/lib gism.cpp -o gism -lsdsl -ldivsufsort -ldivsufsort64
