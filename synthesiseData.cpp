// TO COMPILE
// g++ synthesiseData.cpp -o synthesise -std=c++11
// TO RUN
//fatima@fatima-VirtualBox:~/gism$ ./synthesise -n 30 -d 10 -Smax 5 -Lmax 5 -o text

#include <iostream> //cout,endl
#include <string> //string,
#include <vector> //vector,push_back
#include <fstream> //ifstream,is_open,good,getline,close
#include <sstream> //
#include <ctime> //
#include <cstdlib> //
#include <math.h> //ceil
#include <algorithm> //sort, find

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
2. for each i in s, set random l <= L_max (range begins from 0, where 0 is epsilon)
3. pickRandomBase() l times, s times 
pickRandomBase():
1. randomNumberGenerator for int x = {0, 1, 2, 3}
2. return DNA[x]
*/

std::string pickRandomBase(std::vector<std::string> *dna);
std::string createSet(int *Smax, int *Lmax, std::vector<std::string> *dna);

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
//...- o : name of output file
std::string o = argv[10]; 
std::cout << "output file name: " << o << std::endl;

/* declare other variables */
// a vector of strings holding the DNA alphabet
std::vector<std::string> dna = {"A", "G", "C", "T"};
// int vector D of random positions in T (totalling d% of T) which should be sets instead of singleton chars
std::vector<int> D;

/* set seed for random number generation */
std::srand(std::time(NULL));

/* build D */
//convert % into number dprime out of n
float dprime = ceil( ( float(n) / 100 ) * d );
std::cout << "dprime is " << dprime << std::endl;
//for 0 to dprime-1
//...push back random number between 0 to n-1 to D
for (int i=0; i<dprime; i++){
	D.push_back( std::rand() % n );
}
//sort D
std::sort(D.begin(), D.end());
//print D
std::cout << "vector D: ";
for (std::vector<int>::iterator it=D.begin(); it!=D.end(); it++){
	std::cout << *it << " ";
}
std::cout << std::endl;

/* build T */
std::ofstream file;
file.open(o);
for (int x=0; x<n; x++){
	if ( std::find(D.begin(), D.end(), x) != D.end() ){ //present in D
		file << createSet(&Smax, &Lmax, &dna);
	} else {
		file << pickRandomBase(&dna);
	}
}
file.close();


return 0;
}

std::string createSet(int *Smax, int *Lmax, std::vector<std::string> *dna)
{
//set random s <= Smax (range begins from 2)
int s = ( std::rand() % ( (*Smax) -1) ) + 2;
if (s >= (*Smax)) s=(*Smax);
//
std::stringstream ss;
ss << "{";
for (int i = 0; i<s; i++){
	int l = std::rand() % (*Lmax);
	if (l==0){
		ss << "E";
	} else {
		for (int j = 0; j<l; j++) ss << pickRandomBase(dna);
	}
	ss << ",";
}
std::string Ti = ss.str();
Ti.pop_back();
Ti += "}";
return Ti;
}

std::string pickRandomBase(std::vector<std::string> *dna)
{
int number = ( std::rand() % 4 );
return (*dna)[number];
}
