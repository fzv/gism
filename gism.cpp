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
#include <algorithm> //std::find

/************************************************************************************/
/******************************* FUNCTION DECLARATIONS ******************************/
/************************************************************************************/

void computeBorderTable(std::string *X, std::vector<int> *B);
std::vector<int> computeBorder(std::string temp, std::vector<int> B);
void preKMP(std::string pattern, int f[]);
bool KMP(std::string needle, std::string haystack);
void computeBps(std::vector<std::vector<int>> *L, std::vector<int> *B, std::vector<int> *Bprime, std::string *P);
sdsl::csa_bitcompressed<> computeSuffixArray(std::string s);
void printSuffixArray(sdsl::csa_bitcompressed<> SA);
void printVector(std::vector<int> *vector);
std::vector<int> computeLCParray(std::string s, sdsl::csa_bitcompressed<> SA, std::vector<int> iSA, std::vector<int> LCP);
int getlcp(int suffx, int suffy, std::vector<int> iSA, std::vector<int> LCP, sdsl::rmq_succinct_sct<> rmq);
bool checkL(int value, std::vector<std::vector<int>> *L, int i);
void printL(std::vector<std::vector<int>> *L);
std::vector<std::vector<std::vector<int>>> insertL(int value, std::vector<std::vector<int>> L, int i, int S_j);
void maintainL(std::vector<std::vector<int>> *L, int i);
void printSeqs(std::list<std::vector<std::string>> *T, std::string *P);
void updateR(std::vector<bool> *R, std::vector<std::vector<int>> *L, int m, int i);

/***********************************************************************************/
/************************************ GISM *****************************************/
/***********************************************************************************/
int main()
{
/*********************            Welcome message         **************************/

std::cout << "GISM - Generalised Indeterminate String Matching" << std::endl << std::endl;

/*********************            Parse input file          **************************/

std::string P;
std::list<std::vector<std::string>> T;

{ //block

std::string line;
std::vector<std::string> lines;
std::ifstream inputFile("testdata1");

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

P = lines[1];
std::string t = lines[3];

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
	} //end_if
} //end_for

} //end_block

/*********************            Print input sequences         **************************/
printSeqs(&T, &P);

/*********************                       GISM                  **************************/

std::vector<int> report;
std::vector<std::vector<int>> L; //as defined in paper

for (std::list<std::vector<std::string>>::iterator it=T.begin(); it!=T.end(); it++){ //for each pos T[i] in T 
	int i = std::distance(T.begin(),it); //we are in T[i]
	std::cout << "\n\nwe are in pos " << i << " of T......" << std::endl;
	std::stringstream x; //stringstream used to create string X
	std::vector<int> B; //border table
	std::vector<int> Bprime; //B'[j] = i s.t. i is ending pos of S_j in X
	bool epsilon = false; //flag empty string in T[i]
	x << P << "$"; //concatenate unique symbol to P to initiaise string X
	for (std::vector<std::string>::iterator j=(*it).begin(); j!=(*it).end(); j++){ //for each S_j in T[i]
		if ((*j) == "E") epsilon=true; //if S_j is empty string, set flag to true 
		x << (*j) << "$"; //concatenate S_j and unique letter to string X
		if ((*j).length() >= P.length()){ //if P could occur in S_j
			if (KMP(P, (*j))){ //if P occurs in S_j
				report.push_back(i); //report pos T[i]
				std::cout << "reporting " << i << std::endl;
			}
		}
		Bprime.push_back(x.str().length()-2); //in B': store ending pos of S_j in X
	}
	std::string X = x.str(); //concatenation of P and all S_j, separated by unique chars
	X.pop_back(); //remove unecessary unique letter at end pos of X
	std::cout << "String X: " << X << std::endl; //print string X
	printVector(&Bprime); //print vector B'
	std::cout << "check above indexes of below border table" << std::endl;
	if (it==T.begin()) { //only for T[0] do:
		computeBorderTable(&X, &B);
		computeBps(&L, &B, &Bprime, &P);
	} else {
		computeBorderTable(&X, &B);
		computeBps(&L, &B, &Bprime, &P);
		//construct suffix array of X
		sdsl::csa_bitcompressed<> SA = computeSuffixArray(X);
		printSuffixArray(SA);
		int size = SA.size();
		//construct inverse suffix array of X
		std::vector<int> iSA(size, 0);
		for (int r = 0; r != size; r ++) iSA[SA[r]] = r;
		std::cout << "inverse suffix array" << std::endl; printVector(&iSA);
		//construct longest common prefix array of X + prepare for rmq
		std::vector<int> LCP(size, 0);
		LCP = computeLCParray(X, SA, iSA, LCP);
		sdsl::rmq_succinct_sct<> rmq;
		rmq = sdsl::rmq_succinct_sct<>(&LCP);
		std::cout << "longest common prefix array" << std::endl; printVector(&LCP);
		// initialise bitvector A to aid extension of prefixes of P
		std::vector<bool> A(P.length(),true);
		// initialise bitvector R to aid reporting of occurences of P in T
		std::vector<bool> R(P.length(),false);
		//
		int len; //to store length of S_j
		int cumulative_len = 0; //length of X minus length of S_j
		for (int b = 0; b<Bprime.size(); b++){ //for all S_j in T[i]

			len = Bprime[b] - P.length() - cumulative_len - b; //length of S_j
			int suffs = Bprime[b]-len+1; //start pos of S_j in X
			std::cout << std::endl << X.substr(suffs,len) << std::endl; //extract+print S_j from X
			std::cout << "is of length " << len << std::endl;
			cumulative_len += len; //update cum. length in preparation for next S_j

			if (len < P.length()) { //if S_j could occur in P
				std::cout << "length of S_j is less than P" << std::endl;
				for (int suffp = 1; suffp < P.length(); suffp++){ //for each suffix of P
					int lcp = getlcp(suffp, suffs, iSA, LCP, rmq); //lcp of S_j and suffix of P
					std::cout << "\nlcp of suffixes " << suffp << " and " << suffs << " is " << lcp << std::endl;
					if (lcp >= len){
						std::cout << "S_j occurs in P" << std::endl;
						std::cout << "checking previous pos of P in L[i-1]" << std::endl;
						if (checkL(suffp-1, &L, i-1)){ //check if can extend prefix of P from L[i-1]
							std::cout << "can extend to pos ";
							int endpos = suffp+len-1; //can be extended to endpos in P
							std::cout << endpos << std::endl;
							if (A[endpos]==true){ //if endpos not in L[i]
								std::cout << "not already added to L" << std::endl;
								A[endpos]=false; //do not allow to add endpos to L[i] again
								L[i].push_back(endpos); //add endpos to L[i]
							}
						} else { std::cout << "unable to extend a bit" << std::endl; }
					} //end_if lcp>=len
				} //end_for each suffix of P
			} //end_if S_j could occur in P
			//for all S_j in T[i]:
			updateR(&R, &L, P.length(), i);
/*
			for (int suffp = 1; suffp < P.length(); suffp++){ //for each suffix of P
				int lcp = getlcp(suffp, suffs, iSA, LCP, rmq); //lcp of S_j and suffix of P
				std::cout << "\nlcp of suffixes " << suffp << " and " << suffs << " is " << lcp << std::endl;
				int p = P.length() - lcp - 1; // 
				std::cout << "checking p = " << p << " in L[i-1]" << std::endl;
				if (checkL(p, &L, i-1)){ //if p occurs in L[i-1]
					std::cout << "reporting " << i << std::endl;
					report.push_back(i); //report T[i]
				} else { std::cout << "unable to extend to end" << std::endl; }
			} //end_for each suffix of P
*/
		} //end_for all S_j in T[i]
	} //end_if T[0]
	if (epsilon==true) L[i].insert(L[i].end(), L[i-1].begin(), L[i-1].end());
	std::cout << std::endl << "printing L outside function" << std::endl;
	printL(&L);
	std::cout << std::endl << std::endl;
} //end_GISM


//** report **//
std::cout << std::endl <<"pattern occurs in text, ending at the following positions" << std::endl;
printVector(&report);
std::cout << std::endl;

return 0;

} //end_main




/************************************************************************************/
/******************************* FUNCTION DEFINITIONS *******************************/
/************************************************************************************/

void updateR(std::vector<bool> *R, std::vector<std::vector<int>> *L, int m, int i)
{
std::cout << "inside updateR function" << std::endl;
for (std::vector<int>::iterator p = (*L)[i-1].begin(); p != (*L)[i-1].end(); p++){
		std::cout << "for p=" << (*p) << "desired j is " << m-(*p)-2 << std::endl;
		if ((*p)!=-1000){
			(*R)[m-(*p)-2] = true;
		}
}
}

/*******************           Print (T)ext and (P)attern       ************************/
void printSeqs(std::list<std::vector<std::string>> *T, std::string *P)
{
std::cout << "string T:" << std::endl;
for (std::list<std::vector<std::string>>::iterator i=(*T).begin(); i!=(*T).end(); i++){
	for (std::vector<std::string>::iterator j=(*i).begin(); j!=(*i).end(); j++){
		std::cout << *j << " ";
	}
	std::cout << std::endl;
}
std::cout << std::endl << "string P:" << std::endl;
std::cout << (*P) << std::endl << std::endl;
}



/*********************                 Print L             **************************/
void printL(std::vector<std::vector<int>> *L)
{
	for (std::vector<std::vector<int>>::iterator i = (*L).begin(); i != (*L).end(); i++)
	{
		std::cout << std::endl << "L[" << std::distance((*L).begin(), i) << "]: ";
		std::vector<int> tempVec = *i;
		printVector(&tempVec);
	}
}

/*********************          Return true if value exists in L_i          **************************/
bool checkL(int value, std::vector<std::vector<int>> *L, int i)
{
	for (std::vector<int>::iterator j = (*L)[i].begin(); j != (*L)[i].end(); j++)
	{
		if ((*j)==value) return 1;
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
//std::cout << "i = " << i << std::endl;
//std::cout << "j = " << j << std::endl;
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
void printVector(std::vector<int> *vector)
{
for (std::vector<int>::iterator it = (*vector).begin(); it != (*vector).end(); it ++){
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
// B: border table
// B': B'[j] = i s.t. i is ending pos of S_j in X
// P: pattern string
void computeBps(std::vector<std::vector<int>> *L, std::vector<int> *B, std::vector<int> *Bprime, std::string *P)
{	
	//std::vector<int> Sj;
	std::vector<int> Bps;
	for (int i = 0; i != (*Bprime).size(); i++)
	{
		int Bi = (*Bprime)[i];
		if ((*B)[Bi] != 0)
		{
			std::cout << "looking at " << Bi << "th pos in B: " << (*B)[Bi] << std::endl;
			for(int x = (*B)[Bi]; x != 0; x = (*B)[x-1]){
				//Sj.push_back(x-1);
				Bps.push_back(x-1);
			}
		}
		else
		{
			//Sj.push_back(-1000);
			Bps.push_back(-1000);
		}
		//Bps.push_back(Sj);
		//Sj.clear();
	}
	(*L).push_back(Bps);
}

/*********************                     Compute border table               **************************/
// given parameters:
// X: string formed from concatenation of P and all S_j in T[i], separated by unique chars
// B: vector to hold border table
void computeBorderTable(std::string *X, std::vector<int> *B)
{
	int m = (*X).length();

	for (int b = 0; b < m; b++) (*B).push_back(0);

	//** algorithm from Maxime's book **//
	int p = 0;
	for (int q = 1; q < m; q++){
		(*B)[q-1] = p;
		while (p >= 0 & (*X)[q]!=(*X)[p]){ 
			if (p==0){
				p = -1;
			} else {
				p = (*B)[p-1];
			}
		}
		p++;
	}
	(*B)[m-1] = p;

	//** print border table **//
	for (std::vector<int>::iterator it = (*B).begin(); it != (*B).end(); it++){
		std::cout << *it << " ";
	}
	std::cout << std::endl << std::endl;
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
