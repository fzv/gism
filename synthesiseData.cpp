// g++ synthesiseData.cpp -o synthesise -std=c++11

#include <iostream> //cout,endl
#include <string> //string,
#include <vector> //vector,push_back
#include <fstream> //ifstream,is_open,good,getline,close
#include <sstream> //
#include <ctime> //
#include <cstdlib> //
#include <math.h> //ceil

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

std::string pickRandomBase(std::vector<std::string> *dna);
void createSet();

int main(int argc, char* argv[])
{

/* user given parameters */
//...-n: number of positions in text T
int n = atoi(argv[2]); 
std::cout << "n is " << n << std::endl;
//...- d : percentage of positions in text T which are degenerate i.e. represent indels
int d = atoi(argv[4]); 
std::cout << "d is " << d << "%" << std::endl;
//...- Smax : maximum size of set at any position T[i] i.e. maximum number of S_j
int Smax = atoi(argv[6]); 
std::cout << "Smax is " << Smax << std::endl;
//...- Lmax : upper bound on length of any string S_j in T[i]
int Lmax = atoi(argv[8]); 
std::cout << "Lmax is " << Lmax << std::endl;

/* declare other variables */
// a vector of strings holding the DNA alphabet
std::vector<std::string> dna = {"A", "G", "C", "T"};
// int vector D of random positions in T (totalling d% of T) which should be sets instead of singleton chars
std::vector<int> D;

/* build D */
//convert % into number dprime out of n
int dprime = n / 100;
std::cout << "dprime is " << dprime << std::endl;
dprime = dprime * d;
std::cout << "dprime is " << dprime << std::endl;
dprime = ceil(dprime);
std::cout << "dprime is " << dprime << std::endl;
//for 0 to dprime-1
//...push back random number between 0 to n-1 to D
for (int i=0; i<dprime; i++){
	D.push_back( std::rand() % n );
}
//print D
std::cout << "printing vector D" << std::endl;
for (std::vector<int>::iterator it=D.begin(); it!=D.end(); it++){
	std::cout << *it << " ";
}
std::cout << std::endl;




//set seed for pickRandomBase()
std::srand(std::time(NULL));
//test adding char to stringstream
std::stringstream s;
s << pickRandomBase(&dna);
std::cout << s.str() << std::endl;



return 0;
}

void createSet()
{
//set random s <= Smax (range begins from 2)
///int s = ( std::rand() % (Smax-2) ) + 2;

}

std::string pickRandomBase(std::vector<std::string> *dna)
{
int number = ( std::rand() % 4 );
std::cout << (*dna)[number] << std::endl;
return (*dna)[number];
}
