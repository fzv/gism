#include <iostream> //cout,endl
#include <string> //string,
#include <vector> //vector,push_back
#include <fstream> //ifstream,is_open,good,getline,close
#include <sstream> //
#include "sdsl/suffix_trees.hpp"
#include "sdsl/util.hpp"
#include <iterator>

/************************************************************************************/
/******************************* FUNCTION DECLARATIONS ******************************/
/************************************************************************************/
std::vector<int> computeBorderTable(std::string X, std::vector<int> B);
std::vector<int> computeBorder(std::string temp, std::vector<int> B);
void preKMP(std::string pattern, int f[]);
bool KMP(std::string needle, std::string haystack);
std::list<std::vector<std::vector<int>>> computeBps(std::list<std::vector<std::vector<int>>> L, std::vector<int> report, std::vector<int> B, std::vector<int> Bprime, std::string P);

/*********************************************************************/
/******************************* GISM ********************************/
/*********************************************************************/
int main()
{

//Welcome message

std::cout << "Welcome to GISM - Generalised Indeterminate String Matching\n" << std::endl;

//Parse input file

std::string line;
std::vector<std::string> lines;
std::ifstream inputFile("testdata");

if (inputFile.is_open()){
	if (inputFile.good()){
		while (getline(inputFile, line)){
			lines.push_back(line);
		}
	}
	inputFile.close();
} else {
	std::cout << "Unable to open file." << std::endl;
}

std::string P = lines[1];
std::string t = lines[3];

std::list<std::vector<std::string>> T;
std::vector<std::string> tempVector;
std::stringstream tempString;

for (int i=0; i<t.length(); i++){ //loop through text string
	if (t[i]=='{'){ //if new pos in T
		tempVector.clear(); //clear vector to hold all s in t[i]
		tempString.str(""); //clear string stream to hold all a in s[j]
		tempString.clear();//clear string stream to hold all a in s[j]
	} else if (t[i]=='}'){ //if reach end of t[i]
		tempVector.push_back(tempString.str()); //add previous s to tempVector
		T.push_back(tempVector); //fill current pos in T with tempVector
	} else if (t[i]==','){ //if new s in t[i]
		tempVector.push_back(tempString.str()); //add previous s to tempVector
		tempString.str(""); //clear string stream to hold all a in s[j]
		tempString.clear(); //clear string stream to hold all a in s[j]
	} else { //if next a in s
		tempString << t[i]; //add a to string stream
	}
}

tempVector.clear();
tempString.str("");
tempString.clear();

//Print input sequences
std::cout << "string T:" << std::endl;
for (std::list<std::vector<std::string>>::iterator i=T.begin(); i!=T.end(); i++){
	tempVector = *i;
	for (std::vector<std::string>::iterator j=tempVector.begin(); j!=tempVector.end(); j++){
		std::cout << *j << " ";
	}
	std::cout << std::endl;
}
std::cout << "\nstring P:" << std::endl;
std::cout << P << std::endl;
std::cout << std::endl;

//Construct Suffix Tree of pattern P
std::string file = "pattern";
sdsl::cst_sct3<> cst;
construct(cst, file, 1);

//Do stuff with STp
std::cout << "number of nodes in suffix tree " << cst.nodes() << std::endl << std::endl;

sdsl::cst_sct3<>::size_type d;
sdsl::cst_sct3<>::node_type v;
sdsl::cst_sct3<>::size_type s;
sdsl::cst_sct3<>::size_type sn;
bool l;
sdsl::cst_sct3<>::size_type lb;
sdsl::cst_sct3<>::size_type rb;
sdsl::cst_sct3<>::size_type c;
sdsl::cst_sct3<>::char_type a;

for (sdsl::cst_sct3<>::const_iterator it = cst.begin(); it!=cst.end(); it++)
{
if(it.visit()==1) //if we have not traversed the subtree rooted at v
{
	v = *it;

	sn = cst.sn(v);
	std::cout << "Suffix number " << sn << std::endl;

	d = cst.node_depth(v);
	std::cout << "Node depth " << d << std::endl;

	s = cst.size(v);
	std::cout << s << " leaves in subtree rooted at v" << std::endl;

	l = cst.is_leaf(v);
	std::cout << "I am a leaf: " << l << std::endl;

	lb = cst.lb(v);
	std::cout << "Index of leftmost leaf in SA " << lb << std::endl;

	rb = cst.rb(v);
	std::cout << "Index of rightmost leaf in SA " << rb << std::endl;

	c = cst.degree(v);
	std::cout << "Number of children " << c << std::endl;
	
	//a = cst.edge(v,1);
	//std::cout << "First letter on edge label from root to v: " << a << std::endl;

	std::cout << "L(v) : ";
	for (int i = 1; i <= cst.depth(v); i++)
	{
		std::cout << cst.edge(v,i);
	}
	std::cout << std::endl;

	std::cout << std::endl;
}
}

//Lemma 2
std::stringstream x;
std::string X;
std::string alpha = "mnoqrsvwxyz"; //asuming no more than 12 s in T[i]
int unique = 0;
std::vector<int> B;
std::vector<int> Bprime;
std::vector<int> report;
std::list<std::vector<std::vector<int>>> L;



for (std::list<std::vector<std::string>>::iterator i=T.begin(); i!=T.end(); i++){
	tempVector = *i;
	x << P << "p";
	for (std::vector<std::string>::iterator j=tempVector.begin(); j!=tempVector.end(); j++){
		x << *j << alpha[unique];
		Bprime.push_back(x.str().length()-2);
		unique++;
	}
	X = x.str();
	X.pop_back();
	//string X is ready
	std::cout << X << " : " << std::endl;
	int m = X.length();
	for (int b = 0; b<Bprime.size(); b++) std::cout << Bprime[b] << " ";
	std::cout << std::endl << "check above indexes of below border table" << std::endl;
	if (i==T.begin())
		{
		B = computeBorderTable(X, B);
		L = computeBps(L, report, B, Bprime, P);
		//if |S| >= m ...
		}
	else
		{
		B = computeBorderTable(X, B);
		L = computeBps(L, report, B, Bprime, P);
		//if |S| < m ...
		//compute Bsp
		//if there exists ...
		//if |S| >= m ...
		}

	//clean up
	Bprime.clear();
	B.clear();
	report.clear();
	x.str("");
	x.clear();
	unique = 0;
	std::cout << std::endl << std::endl;
}







return 0;
}
/************************************************************************************/
/******************************* FUNCTION DEFINITIONS *******************************/
/************************************************************************************/
std::list<std::vector<std::vector<int>>> computeBps(std::list<std::vector<std::vector<int>>> L, std::vector<int> report, std::vector<int> B, std::vector<int> Bprime, std::string P)
{
	int Bi;	
	std::vector<int> Sj;
	std::vector<std::vector<int>> Bps;
	for (int i = 0; i != Bprime.size(); i++)
	{
		Bi = Bprime[i];
		if (B[Bi] != 0)
		{
			std::cout << "looking at " << Bi <<"th pos in B: " << B[Bi] << std::endl;
		}
		else
		{
			Sj.push_back(-1);
		}
		Bps.push_back(Sj);
	}
	L.push_back(Bps);

	//print L
	for (int i = 0; i != L.size(); i++)
	{
		std::cout << "L[" << i << "]" << std::endl;
		for (int j = 0; j != i.size(); j++)
		{
			for (int k = 0; k != j.size(); k++)
			{
				
			}
		}
	}

	return L;
}

std::vector<int> computeBorderTable(std::string X, std::vector<int> B)
{
	int m = X.length();

	for (int b = 0; b < m; b++) B.push_back(0);

	//algorithm from Maxime's book
	int p = 0;
	for (int q = 1; q < m; q++){
		B[q-1] = p;
		while (p >= 0 & X[q]!=X[p]){ //what does E match with??
			if (p==0){
				p = -1;
			} else {
				p = B[p-1];
			}
		}
		p++;
	}
	B[m-1] = p;

	//print border table
	for (std::vector<int>::iterator it = B.begin(); it != B.end(); it++){
		std::cout << *it << " ";
	}
	std::cout << std::endl << std::endl;

	return B;
}

//credit to http://www.sanfoundry.com/cpp-program-implement-kruth-morris-patt-algorithm-kmp/

void preKMP(std::string pattern, int f[])
{
	int m = pattern.length();
	int k;
	f[0] = -1;
	for (int i = 1; i < m; i++)
	{
		k = f[i - 1];
		while (k >= 0)
		{
			if (pattern[k] == pattern[i - 1])
			{
				break;
			}
			else
			{
				k = f[k];
			}
		}
		f[i] = k + 1;
	}

}

bool KMP(std::string needle, std::string haystack)
{
	int m = needle.length();
	int n = haystack.length();
	int f[m];
	preKMP(needle, f);
	int i = 0;
	int k = 0;
	while (i<n)
	{
		if (k==-1)
		{
			i++;
			k = 0;
		}
        	else if (haystack[i] == needle[k])
		{
			i++;
			k++;
			if (k==m) return 1;
		}
		else
		{
			k = f[k];
		}
	}
	return 0;
}
