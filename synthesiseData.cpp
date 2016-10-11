#include <iostream> //cout,endl
#include <string> //string,
#include <vector> //vector,push_back
#include <fstream> //ifstream,is_open,good,getline,close
#include <sstream> //
#include <ctime> //
#include <cstdlib> //

int main()
{

std::srand(std::time(NULL));

return 0;
}

/*

This script takes parameters...
- n : number of positions in text T
- d% : percentage of positions in text T which are degenerate i.e. represent indels
- S_max : maximum size of set at any position T[i] i.e. maximum number of S_j
- L_max : upper bound on length of any string S_j in T[i]
...and makes use of...
- DNA : a vector of chars holding the DNA alphabet
...in order to...
1. define an int vector D of random positions in T (totalling d% of T) which should be sets instead of singleton chars
2. from 0 to n-1:
	if pos in D then:
		createSet()
	else:
		pickRandomBase()
...where...
createSet():
1. set random s <= S_max (range begins from 2)
2. for each s, set random l <= L_max (range begins from 0, where 0 is epsilon)
3. pickRandomBase() l times, s times 
pickRandomBase():
1. randomNumberGenerator for int x = {0, 1, 2, 3}
2. return DNA[x]
*/
