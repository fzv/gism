#include <iostream> //cout,endl
#include <string> //string,
#include <vector> //vector,push_back
#include <fstream> //ifstream,is_open,good,getline,close
#include <sstream> //
#include "sdsl/suffix_trees.hpp"
#include "sdsl/suffix_arrays.hpp"
#include "sdsl/lcp.hpp"
#include "sdsl/util.hpp"
#include <iterator>
#include "sdsl/rmq_support.hpp"

/************************************************************************************/
/******************************* FUNCTION DECLARATIONS ******************************/
/************************************************************************************/

std::vector<int> computeBorderTable(std::string X, std::vector<int> B);
std::vector<int> computeBorder(std::string temp, std::vector<int> B);
void preKMP(std::string pattern, int f[]);
bool KMP(std::string needle, std::string haystack);
std::vector<std::vector<std::vector<int>>> computeBps(std::vector<std::vector<std::vector<int>>> L, std::vector<int> report, std::vector<int> B, std::vector<int> Bprime, std::string P);
sdsl::csa_bitcompressed<> computeSuffixArray(std::string s);
void printSuffixArray(sdsl::csa_bitcompressed<> SA);
void printVector(std::vector<int> vector);
std::vector<int> computeLCParray(std::string s, sdsl::csa_bitcompressed<> SA, std::vector<int> iSA, std::vector<int> LCP);
int getlcp(int suffx, int suffy, std::vector<int> iSA, std::vector<int> LCP, sdsl::rmq_succinct_sct<> rmq);
bool checkL(int value, std::vector<std::vector<std::vector<int>>> L, int i);
void printL(std::vector<std::vector<std::vector<int>>> L);
std::vector<std::vector<std::vector<int>>> insertL(int value, std::vector<std::vector<std::vector<int>>> L, int i, int S_j);

/***********************************************************************************/
/************************************ GISM *****************************************/
/***********************************************************************************/
int main()
{
{
std::string X = "banana";
		sdsl::csa_bitcompressed<> SA = computeSuffixArray(X);
		printSuffixArray(SA);
		int size = SA.size();
		std::vector<int> iSA(size, 0);
		for (int i = 0; i != size; i ++) iSA[SA[i]] = i;
		printVector(iSA);
		std::vector<int> LCP(size, 0);
		LCP = computeLCParray(X, SA, iSA, LCP);
		printVector(LCP);
}
/*********************            Welcome message         **************************/

std::cout << "GISM - Generalised Indeterminate String Matching" << std::endl << std::endl;

/*********************            Parse input file          **************************/

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

/*********************            Print input sequences         **************************/
std::cout << "string T:" << std::endl;
for (std::list<std::vector<std::string>>::iterator i=T.begin(); i!=T.end(); i++){
	tempVector = *i;
	for (std::vector<std::string>::iterator j=tempVector.begin(); j!=tempVector.end(); j++){
		std::cout << *j << " ";
	}
	std::cout << std::endl;
}
std::cout << std::endl << "string P:" << std::endl;
std::cout << P << std::endl << std::endl;

/*********************                       GISM                  **************************/
std::vector<int> report;

{ //gism
std::stringstream x;
std::string X; // concatenation of P and all S_j, separated by unique chars
std::vector<int> B; // border table
std::vector<int> Bprime; //B'[j] = i s.t. i is ending pos of S_j in X
std::vector<std::vector<std::vector<int>>> L; //as defined in paper

for (std::list<std::vector<std::string>>::iterator i=T.begin(); i!=T.end(); i++){ //for each pos T[i] in T 
	std::vector<std::string> tempVector = *i; //tempVector holds list of all S_j in T[i]
	x << P << "$"; //concatenate unique symbol to P to form string X
	for (std::vector<std::string>::iterator j=tempVector.begin(); j!=tempVector.end(); j++){ //for each S_j in T[i]
		std::string S_j = *j;
		x << S_j << "$"; //concatenate S_j and unique letter to string X
		if (S_j.length() >= P.length()){ //then P could occur in S_j
			if (KMP(P, S_j)==1) report.push_back(std::distance(T.begin(),i)); //if P occurs in S_j, report pos T[i]
		}
		Bprime.push_back(x.str().length()-2); //in B': store ending pos of S_j in X
	}
	X = x.str(); //convert stringstream x into string X
	X.pop_back(); //remove unecessary unique letter at end pos of X
	std::cout << "String X: " << X << std::endl; //print string X
	for (int b = 0; b<Bprime.size(); b++) std::cout << Bprime[b] << " "; //print vector B'
	std::cout << std::endl << "check above indexes of below border table" << std::endl;
	if (i==T.begin())
		{
		B = computeBorderTable(X, B);
		L = computeBps(L, report, B, Bprime, P);
		}
	else
		{
		B = computeBorderTable(X, B);
		L = computeBps(L, report, B, Bprime, P);
		//
		sdsl::csa_bitcompressed<> SA = computeSuffixArray(X);
		printSuffixArray(SA);
		int size = SA.size();
		std::vector<int> iSA(size, 0);
		for (int i = 0; i != size; i ++) iSA[SA[i]] = i;
		std::cout << "inverse suffix array" << std::endl; printVector(iSA);
		std::vector<int> LCP(size, 0);
		LCP = computeLCParray(X, SA, iSA, LCP);
		//
		sdsl::rmq_succinct_sct<> rmq;
		rmq = sdsl::rmq_succinct_sct<>(&LCP);
		//
		std::cout << "longest common prefix array" << std::endl; printVector(LCP);
		//
		std::vector<bool> A(P.length(),true);
		//
		int len;
		int cumulative_len = 0;
		for (int b = 0; b<Bprime.size(); b++){ //for all S_j in T[i]
			len = Bprime[b] - P.length() - cumulative_len - b;
			int suffs = Bprime[b]-len+1;
			std::cout << std::endl << X.substr(suffs,len) << std::endl;
			std::cout << "is of length " << len << std::endl;
			cumulative_len += len;
			if (len < P.length()){
				for (int suffp = 0; suffp < P.length(); suffp++){
					int lcp = getlcp(suffp, suffs, iSA, LCP, rmq);
					std::cout << "lcp of suffixes " << suffp << " and " << suffs << " is " << lcp << std::endl;
					if (lcp == 0){
						//do nothing
					} else if (lcp >= len){ //S_j occurs in P
						lcp = len;
						int Li = std::distance(T.begin(),i);
						if (checkL(suffp-1, L, Li-1)){
							std::cout << "can extend to pos ";
							int endpos = suffp+lcp-1;
							std::cout << endpos << std::endl;
							if (A[endpos]==true){
								std::cout << "not already added to L" << std::endl;
								A[endpos]=false;
								L[Li][b].push_back(endpos);
							}
						}
						//check L_i-1 for value suffp-1
							//if yes 1) turn off A[endpos of p] 2) add endpos of p to L_i
					} else { //prefix of S_j is a suffix of P
						//
					}
				}
			}
		}
		printL(L);
		//check all j in A
		//compute Bsp
		//if there exists ..
		}
	//** report **//
	for(std::vector<int>::iterator it = report.begin(); it != report.end(); it++) std::cout << *it << " ";
	//** clean up **//
	Bprime.clear();
	B.clear();
	report.clear();
	x.str("");
	x.clear();
	std::cout << std::endl << std::endl;
}

} //end_gism

return 0;

} //end_main




/************************************************************************************/
/******************************* FUNCTION DEFINITIONS *******************************/
/************************************************************************************/


/*********************                 Print L             **************************/
void printL(std::vector<std::vector<std::vector<int>>> L)
{
	int inti = 0;
	int intj = 0;
	for (std::vector<std::vector<std::vector<int>>>::iterator i = L.begin(); i != L.end(); i++)
	{
		std::cout << std::endl << "L[" << inti << "]" << std::endl;
		std::vector<std::vector<int>> tempVecVec = *i;
		inti++;
		for (std::vector<std::vector<int>>::iterator j = tempVecVec.begin(); j != tempVecVec.end(); j++)
		{
			std::cout << "  S_" << intj << " : ";
			std::vector<int> tempVec = *j;
			intj++;
			for (std::vector<int>::iterator k = tempVec.begin(); k != tempVec.end(); k++)
			{
				int tempInt = *k;
				std::cout << tempInt << " "; //ending pos of pref(P) = suff(Sj)
			}
		}
	}
}

/*********************          Return true if int exists in L_i          **************************/
bool checkL(int value, std::vector<std::vector<std::vector<int>>> L, int i)
{
	std::vector<std::vector<int>> tempVecVec = L[i];
	for (std::vector<std::vector<int>>::iterator j = tempVecVec.begin(); j != tempVecVec.end(); j++)
	{
		std::vector<int> tempVec = *j;
		for (std::vector<int>::iterator k = tempVec.begin(); k != tempVec.end(); k++)
		{
			int tempInt = *k;
			if (tempInt==value) return 1;
		}
	}
return 0;
}



/*********************          Compute lcp(suffix x, suffix y)          **************************/
int getlcp(int suffx, int suffy, std::vector<int> iSA, std::vector<int> LCP, sdsl::rmq_succinct_sct<> rmq)
{
int i;
int j;
if (iSA[suffx] < iSA[suffy]){
	i = iSA[suffx];
	j = iSA[suffy];
} else {
	i = iSA[suffy];
	j = iSA[suffx];
}
std::cout << "i = " << i << std::endl;
std::cout << "j = " << j << std::endl;
auto min_idx = rmq(i+1,j); 
int lcp = LCP[min_idx];
return lcp;
}

/*********************          Compute LCP array of string s         **************************/
std::vector<int> computeLCParray(std::string s, sdsl::csa_bitcompressed<> SA, std::vector<int> iSA, std::vector<int> LCP)
{
s = s.append("$");
int n = s.length(); 
int lcp = 0;
for (int i = 0; i < n; i++){
	if (iSA[i] == n-1){
		lcp = 0;
		continue;
	}
	int j = SA[iSA[i]+1];
	while ( ( (i+lcp) < n ) && ( (j+lcp) < n ) && ( s[i+lcp]==s[j+lcp] ) ) lcp++;
	LCP[iSA[i]] = lcp;
	if (lcp > 0) lcp--;
}
std::vector<int>::iterator it = LCP.begin();
LCP.insert(it, 0);
LCP.pop_back();

return LCP;
}

/*********************          Print any int vector         **************************/
void printVector(std::vector<int> vector)
{
for (std::vector<int>::iterator it = vector.begin(); it != vector.end(); it ++){
	std::cout << *it << " "; 
}
std::cout << std::endl << std::endl;
}

/*********************          Compute suffix array of string s         **************************/
sdsl::csa_bitcompressed<> computeSuffixArray(std::string s)
{
sdsl::csa_bitcompressed<> SA;
construct_im(SA, s, 1);
return SA;
}

/*********************          Print suffix array SA         **************************/
void printSuffixArray(sdsl::csa_bitcompressed<> SA){
std::cout << std::endl << "Suffix array" << std::endl;
for (sdsl::csa_bitcompressed<>::iterator it = SA.begin(); it != SA.end(); it ++){
	std::cout << *it << " "; 
}
std::cout << std::endl << std::endl;
}

/*********************                     Compute B_p,s               **************************/
// B_p,s stores all prefixes of pattern P that are suffixes of S_j, given parameters:
// L: stores all B_p,s (for all pos in T)
// report: stores pos T[i] to be reported as output
// B: border table
// B': B'[j] = i s.t. i is ending pos of S_j in X
// P: pattern string
std::vector<std::vector<std::vector<int>>> computeBps(std::vector<std::vector<std::vector<int>>> L, std::vector<int> report, std::vector<int> B, std::vector<int> Bprime, std::string P)
{	
	std::vector<int> Sj;
	std::vector<std::vector<int>> Bps;
	for (int i = 0; i != Bprime.size(); i++)
	{
		int Bi = Bprime[i];
		if (B[Bi] != 0)
		{
			std::cout << "looking at " << Bi << "th pos in B: " << B[Bi] << std::endl;
			for(int x = B[Bi]; x != 0; x = B[x-1]){
				Sj.push_back(x-1);
			}
		}
		else
		{
			Sj.push_back(-1);
		}
		Bps.push_back(Sj);
		Sj.clear();
	}
	L.push_back(Bps);

	//** print L **//
	int inti = 0;
	int intj = 0;
	for (std::vector<std::vector<std::vector<int>>>::iterator i = L.begin(); i != L.end(); i++)
	{
		std::cout << std::endl << "L[" << inti << "]" << std::endl;
		std::vector<std::vector<int>> tempVecVec = *i;
		inti++;
		for (std::vector<std::vector<int>>::iterator j = tempVecVec.begin(); j != tempVecVec.end(); j++)
		{
			std::cout << "  S_" << intj << " : ";
			std::vector<int> tempVec = *j;
			intj++;
			for (std::vector<int>::iterator k = tempVec.begin(); k != tempVec.end(); k++)
			{
				int tempInt = *k;
				std::cout << tempInt << " "; //ending pos of pref(P) = suff(Sj)
			}
		}
	}

	return L;
}

/*********************                     Compute border table               **************************/
// given parameters:
// X: string formed from concatenation of P and all S_j in T[i], separated by unique chars
// B: vector to hold border table
std::vector<int> computeBorderTable(std::string X, std::vector<int> B)
{
	int m = X.length();

	for (int b = 0; b < m; b++) B.push_back(0);

	//** algorithm from Maxime's book **//
	int p = 0;
	for (int q = 1; q < m; q++){
		B[q-1] = p;
		while (p >= 0 & X[q]!=X[p]){ //what does E match with?? what about reporting occ where E in T[i...j]?
			if (p==0){
				p = -1;
			} else {
				p = B[p-1];
			}
		}
		p++;
	}
	B[m-1] = p;

	//** print border table **//
	for (std::vector<int>::iterator it = B.begin(); it != B.end(); it++){
		std::cout << *it << " ";
	}
	std::cout << std::endl << std::endl;

	return B;
}

/*********************                    KMP string matching algorithm              **************************/
//credit to http://www.sanfoundry.com/cpp-program-implement-kruth-morris-patt-algorithm-kmp/
//returns 1 if needle found in haystack, else returns 0

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






















/*********************                Do stuff with SA(C)          **************************/
/*
std::cout << std::endl << "sigma" << std::endl;
std::cout << SA.sigma << std::endl;
std::cout << std::endl << "text" << std::endl;
std::cout << SA.text << std::endl;
*/




/*
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
*/


/*********************            (C)oncatenate all strings         **************************/
/*
std::string C; 
std::vector<int> Cprime;
// C'[i] = k s.t. C[k] is starting pos of S_k, where max(k) = total number of S_j in all T 

{ //concatenation

std::stringstream Cstream;
int k = 0;
std::string unique = "bdfhijklmnopqrsuvwxyz"; //assumption: no more than 20 S_k in T

Cstream << P;
for (std::list<std::vector<std::string>>::iterator i = T.begin(); i != T.end(); i++){
	std::vector<std::string> T_i = *i;
	for (std::vector<std::string>::iterator j = T_i.begin(); j != T_i.end(); j++){
		std::string S_j = *j;
		Cprime.push_back(Cstream.str().length()+1);
		Cstream << unique[k] << S_j;
		k++;
	}
}
C = Cstream.str();

} //end_concatenation

std::cout << C << std::endl;
*/

/*********************                     Print C'               **************************/
/*
std::cout << std::endl << std::endl;
for (std::vector<int>::iterator it = Cprime.begin(); it != Cprime.end(); it ++){
	std::cout << *it << " "; 
}
std::cout << std::endl << std::endl;
*/

/*********************           Construct suffix array of C       **************************/
//sdsl::csa_bitcompressed<> SA = computeSuffixArray(C);

/*********************                    Print SA(C)             **************************/
//printSuffixArray(SA);

/*********************      Construct inverse suffix array of C   **************************/
/*
int size = SA.size();
std::vector<int> iSA(size, 0);
for (int i = 0; i != size; i ++) iSA[SA[i]] = i;
*/

/*********************                    Print iSA(C)             **************************/
//std::cout << std::endl << "iSA of C" << std::endl;
//printVector(iSA);

/*********************   Construct LCP array of C - Kasai's algorithm  **************************/
//std::vector<int> LCP(size, 0);
//LCP = computeLCParray(C, SA, iSA, LCP);

/*********************                    Print LCP(C)             **************************/
//std::cout << std::endl << "LCP array" << std::endl;
//printVector(LCP);
