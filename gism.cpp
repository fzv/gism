//USAGE
//fatima@fatima-VirtualBox:~/gism$ ./gism testdata1


//security of the program
//emory consumption
//speed of program

//indentation style
//blocks of code beginning with comments
//consistent naming scheme
//avoid deep nesting
//limit line length
//consistent names for temp variables e.g. i and j
//SYMBOLIC_NAMES instead of numbers
//to use or not to use: using namespace std; and sdsl;
//dont use global variables except to communicate between functions
//classes needed?

//README
/*
install SDSL lite library
*/


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
#include <ctime> 

/************************************************************************************/
/******************************* FUNCTION DECLARATIONS ******************************/
/************************************************************************************/

void parseInput(std::string *P, std::list<std::vector<std::string>> *T, std::string textfile, std::string patfile, std::ofstream *tempfile);
void prepareX(std::stringstream *x, std::vector<int> *Bprime, bool *epsilon, std::string *P, std::vector<int> *report, int *i, std::string *X, std::list<std::vector<std::string>>::iterator it, std::vector<int> *f, int *prev_position, bool *deg, bool *reporter);
void computeBorderTable(std::string *X, std::vector<int> *B);
void preKMP(std::string *pattern, std::vector<int> *f);
//bool KMP(std::string *needle, std::string *haystack);
void KMP(std::string *needle, std::string *haystack, std::vector<int> *f, bool *deg, std::vector<int> *report, int *prev_position, bool *reporter);
void computeBps(std::vector<int> *Li, std::vector<int> *B, std::vector<int> *Bprime, std::string *P);
sdsl::csa_bitcompressed<> computeSuffixArray(std::string s);
//void computeLCParray(std::string *s, sdsl::csa_bitcompressed<> *SA, std::vector<int> *iSA, std::vector<int> *LCP);
//int getlcp(int *suffx, int *suffy, std::vector<int> *iSA, std::vector<int> *LCP, sdsl::rmq_succinct_sct<> *rmq);
//////////////////////int getlcp(int *suffx, int *suffy, std::vector<int> *iSA, std::vector<int> *LCP);
void updateBitVector(std::vector<bool> *BV, std::vector<int> *Li_1, int m);
void reporting(std::vector<int> *vector, std::ofstream *outfile);
void extend(std::string *sj, sdsl::cst_sada<> *stp, sdsl::csa_bitcompressed<> *sap, std::vector<int> *Li, std::vector<bool> *E, std::vector<bool> *PREV, int *len);
void extendToEnd(std::string *sj, sdsl::cst_sada<> *stp, sdsl::csa_bitcompressed<> *sap, std::vector<int> *Li, std::vector<bool> *R, int *len, std::vector<int> *report, int *prev_position, std::string *P, bool *deg, bool *reporter);

// testing purposes only //
void printSeqs(std::list<std::vector<std::string>> *T, std::string *P);
void printVector(std::vector<int> *vector);
void printSuffixArray(sdsl::csa_bitcompressed<> SA);



/***********************************************************************************/
/************************************ GISM *****************************************/
/***********************************************************************************/
int main(int argc, char* argv[])
{

/* Welcome Message */

std::cout << "GISM - Generalised Indeterminate String Matching" << std::endl;

/* Declare & Assign Seq Variables */

std::string P;
std::ofstream tempfile; 
std::list<std::vector<std::string>> T;
std::string textfile = argv[1];
///std::cout << "Text file name: " << textfile << std::endl;
std::string patfile = argv[2];
///std::cout << "Pattern file name: " << patfile << std::endl;
parseInput(&P, &T, textfile, patfile, &tempfile);
///std::cout << "there are " << T.size() << " positions in T" << std::endl;
//printSeqs(&T, &P);
//pre-process pattern ready for KMPs


std::clock_t start_time = clock();


std::vector<int> f(P.length(),0);
//int f[P.length()];
preKMP(&P, &f);


//Construct Suffix Tree of pattern P
//std::string file = "pattern";
sdsl::cst_sada<> stp; //documentation says construction is slow but fast operations, compared to other = vice versa
construct(stp, "temporary.gism", 1); //build suffix tree using pattern file
std::remove("temporary.gism");

/*
std::cout << "number of nodes in suffix tree " << stp.nodes() << std::endl << std::endl;

sdsl::cst_sada<>::size_type d;
sdsl::cst_sada<>::node_type v;
sdsl::cst_sada<>::size_type s;
sdsl::cst_sada<>::size_type sn;
bool l;
sdsl::cst_sada<>::size_type lb;
sdsl::cst_sada<>::size_type rb;
sdsl::cst_sada<>::size_type c;
sdsl::cst_sada<>::char_type a;

for (sdsl::cst_sada<>::const_iterator it = stp.begin(); it!=stp.end(); it++)
{
if(it.visit()==1) //if we have not traversed the subtree rooted at v
{
	v = *it;

	sn = stp.sn(v);
	std::cout << "Suffix number " << sn << std::endl;

	d = stp.node_depth(v);
	std::cout << "Node depth " << d << std::endl;

	s = stp.size(v);
	std::cout << s << " leaves in subtree rooted at v" << std::endl;

	l = stp.is_leaf(v);
	std::cout << "I am a leaf: " << l << std::endl;

	lb = stp.lb(v);
	std::cout << "Index of leftmost leaf in SA " << lb << std::endl;

	rb = stp.rb(v);
	std::cout << "Index of rightmost leaf in SA " << rb << std::endl;

	c = stp.degree(v);
	std::cout << "Number of children " << c << std::endl;
	
	a = stp.edge(v,1);
	std::cout << "First letter on edge label from root to v: " << a << std::endl;

	std::cout << "L(v) : ";
	for (int i = 1; i <= stp.depth(v); i++)
	{
		std::cout << stp.edge(v,i);
	}
	std::cout << std::endl;

	std::cout << std::endl;
}
}
*/



//Construct Suffix Array of pattern P
sdsl::csa_bitcompressed<> sap = computeSuffixArray(P);
///std::cout << "printing suffix array" << std::endl;
///printSuffixArray(sap);


/* GISM */
std::vector<int> report; //vector storing all reported i of T
std::vector<int> Li;
std::vector<int> Li_1;
int prev_position = -1;

/* Loop Through Each T[i] */
for (std::list<std::vector<std::string>>::iterator it=T.begin(); it!=T.end(); it++){
	/* declare i */
	int i = std::distance(T.begin(),it);
	bool deg = false;
	if ( (*it).size() > 1) deg = true;
	bool reporter = false;
	//std::cout << deg << std::endl;
	std::cout << "\n\nwe are in pos " << i << " of T......" << std::endl;
	/* prepare all S_j in T[i] */
	std::stringstream x; //stringstream used to create string X
	std::vector<int> Bprime; //B'[j] = i s.t. i is ending pos of S_j in X
	bool epsilon = false; //flag empty string in T[i]
	std::vector<int> B; //border table
	std::string X; //as defined in paper
	prepareX(&x, &Bprime, &epsilon, &P, &report, &i, &X, it, &f, &prev_position, &deg, &reporter); // O(X + S_j)
	/* begin GISM algorithm */
	if (it==T.begin()) { //only for T[0] do:
		/* STEP 1: FIND PREFIXES OF P */
		///std::cout << "/* STEP 1: FIND PREFIXES OF P */" << std::endl;
		computeBorderTable(&X, &B); // O(X)
		computeBps(&Li, &B, &Bprime, &P);
		///printVector(&Li);
	} else {
		Li_1 = Li;
		Li.clear();
		/* STEP 1: FIND PREFIXES OF P */
		///std::cout << "/* STEP 1: FIND PREFIXES OF P */" << std::endl;
		computeBorderTable(&X, &B); // O(X)
		computeBps(&Li, &B, &Bprime, &P);
		///std::cout << "printing Li after step 1" << std::endl;
		///printVector(&Li);
		// initialise bitvector (E)xtend to aid extension of prefixes of P
		std::vector<bool> E(P.length(),false);
		updateBitVector(&E, &Li_1, P.length());
		// initialise bitvector (PREV)iously extended to aid extension of prefixes of P
		std::vector<bool> PREV(P.length(),true);
		/* STEP 2: EXTEND PREFIXES OF P */
		///std::cout << "/* STEP 2: EXTEND PREFIXES OF P */" << std::endl;
		int cumulative_len = 0; //length of X minus length of S_j
		for (int b = 0; b<Bprime.size(); b++){ //for all S_j in T[i]
			int len = Bprime[b] - P.length() - cumulative_len - b; //length of S_j
			int startpos = Bprime[b]-len+1; //start pos of S_j in X
			cumulative_len += len; //update cum. length in preparation for next S_j
			if (len < P.length()) { //if S_j could occur in P
				std::string sj = X.substr(startpos,len);
				///std::cout << "\nnow looking at sj = " << sj << std::endl;
				extend(&sj, &stp, &sap, &Li, &E, &PREV, &len);
			} //end_if S_j could occur in P
		} //end_for all S_j in T[i]
		///std::cout << "printing Li after step 2" << std::endl;
		///printVector(&Li);
		/* STEP 3: EXTEND PREFIXES OF P TO END OF P*/
		///std::cout << "/* STEP 3: EXTEND PREFIXES OF P TO END OF P*/" << std::endl;
		std::vector<bool> R(P.length(),false); //bitvector (R)eport to aid reporting of occurences of P in T
		updateBitVector(&R, &Li_1, P.length());
		cumulative_len = 0; //length of X minus length of S_j
		for (int b = 0; b<Bprime.size(); b++){ //for all S_j in T[i]
			int len = Bprime[b] - P.length() - cumulative_len - b; //length of S_j
			int startpos = Bprime[b]-len+1; //start pos of S_j in X
			cumulative_len += len; //update cum. length in preparation for next S_j
			std::string sj = X.substr(startpos,len);
			///std::cout << "now looking at sj = " << sj << std::endl;
			extendToEnd(&sj, &stp, &sap, &Li, &R, &len, &report, &prev_position, &P, &deg, &reporter);
		} //end_for all S_j in T[i]
		///std::cout << "printing Li after step 3" << std::endl;
		///printVector(&Li);



		if (epsilon==true) Li.insert(Li.end(), Li_1.begin(), Li_1.end());
	} //end_if T[0]

	if ((*it).size() > 1) {
		prev_position++;
	} else {
		prev_position += (*it)[0].size();
	}

} //end_GISM

std::clock_t end_time = clock();
double elapsed_secs = double(end_time - start_time) / CLOCKS_PER_SEC;


//** report **//
std::ofstream outfile; 
std::string output_file = argv[3];
outfile.open(output_file);
outfile << elapsed_secs << std::endl;
outfile << "#" << report.size() << std::endl;
//outfile << std::endl << "pattern occurs in text, ending at the following positions" << std::endl;
reporting(&report, &outfile);
outfile << std::endl;
outfile.close();

return 0;

} //end_main




/************************************************************************************/
/******************************* FUNCTION DEFINITIONS *******************************/
/************************************************************************************/

/*******************           extend to end       ************************/


void extendToEnd(std::string *sj, sdsl::cst_sada<> *stp, sdsl::csa_bitcompressed<> *sap, std::vector<int> *Li, std::vector<bool> *R, int *len, std::vector<int> *report, int *prev_position, std::string *P, bool *deg, bool *reporter)
{

std::string s = (*sj);
int occ = -1;
sdsl::cst_sada<>::size_type lb;
sdsl::cst_sada<>::size_type rb;
sdsl::cst_sada<>::node_type v = (*stp).root();
auto it = s.begin();
for (uint64_t char_pos = 0; it != s.end(); ++it){
	if ( forward_search( (*stp), v, it-s.begin(), *it, char_pos) > 0 ){
		occ = it-s.begin();
		///std::cout << "occ = " << occ << std::endl;
		lb = (*stp).lb(v);
		rb = (*stp).rb(v);

		//lb--;
		///std::cout << "lb = " << lb << std::endl;
		//rb--;
		///std::cout << "rb = " << rb << std::endl;
		///std::cout <<  "found " << (*sj).substr(0,occ+1) << std::endl;
		if (occ > -1){ //prefix of S_j found in P, but where in P? Must be suffix!
			for (int i = lb; i < rb+1; i++){
				///std::cout << "occurs at pos " << (*sap)[i] << " in pattern" << std::endl;	
				int startpos = (*sap)[i];
				int endpos = startpos+occ; //can be extended to P[endpos]
				///std::cout << "ends at pos " << endpos << " in the pattern" << std::endl;
				if ( (*R)[startpos] && endpos == ( (*P).length()-1 ) ){
				//if ((*R)[startpos] && (*stp).is_leaf(v)){
					if ((*deg)==true) {
						//if ((*reporter)==false){
							(*report).push_back((*prev_position)+1);
							std::cout << "!!extend reporting pos " << occ << std::endl;
							//(*reporter)==true;
						//}
					} else {
						(*report).push_back((*prev_position)+occ+1);
						std::cout << "!!extend reporting pos " << occ << std::endl;
						
					}
					///std::cout << "reporting pos " << (*tpos) << std::endl;
				}

			}
		} else {
			///std::cout << "not suffix of p" << std::endl;
		}

	} else {
		break;
	}
}


}

/*******************           extend        ************************/


void extend(std::string *sj, sdsl::cst_sada<> *stp, sdsl::csa_bitcompressed<> *sap, std::vector<int> *Li, std::vector<bool> *E, std::vector<bool> *PREV, int *len)
{

std::string s = (*sj);
int occ = -1;
sdsl::cst_sada<>::size_type lb;
sdsl::cst_sada<>::size_type rb;
sdsl::cst_sada<>::node_type v = (*stp).root();
auto it = s.begin();
for (uint64_t char_pos = 0; it != s.end(); ++it){
	if ( forward_search( (*stp), v, it-s.begin(), *it, char_pos) > 0 ){
		occ = it-s.begin();
		lb = (*stp).lb(v);
		rb = (*stp).rb(v);
	} else {
		break;
	}
}
//lb--;
//rb--;
if (occ == (*len)-1){ // whole sj found in P, but where?
	for (int i = lb; i < rb+1; i++){
		///std::cout << "sj occurs at pos " << (*sap)[i] << " in pattern" << std::endl;	
		int startpos = (*sap)[i];
		//int endpos = startpos+(*len)-1; //can be extended to P[endpos]
		int endpos = startpos+occ;
		if ((*E)[startpos] && (*PREV)[endpos]){
			(*PREV)[endpos]=false;
			(*Li).push_back(endpos);
			///std::cout << "adding to Li" << std::endl;
		}

	}
} else {
	///std::cout << "sj DOES NOT occur in p" << std::endl;
}

}

/*******************           report        ************************/
// post-processing
/*
void reporting(std::vector<int> *vector,
		std::ofstream *outfile
		)
{
if ( (*vector).size()==0 ){
	std::cout << "Nothing to report." << std::endl;
} else {
	(*outfile) << (*vector)[0] << std::endl;
	for (int i=1; i<(*vector).size(); i++){
		if((*vector)[i] != (*vector)[i-1]) (*outfile) << (*vector)[i] << std::endl;
	}
}
}
*/
/////FOR COMPARSON WITH RITU ONLY:
void reporting(std::vector<int> *vector,
		std::ofstream *outfile
		)
{
if ( (*vector).size()==0 ){
	std::cout << "Nothing to report." << std::endl;
} else {
	(*outfile) << (*vector)[0];
	for (int i=1; i<(*vector).size(); i++){
		(*outfile) << std::endl << (*vector)[i];
	}
}
}


/*******************           prepare X        ************************/
// O(X)
void prepareX(std::stringstream *x,
		std::vector<int> *Bprime,
		bool *epsilon,
		std::string *P,
		std::vector<int> *report,
		int *i,
		std::string *X,
		std::list<std::vector<std::string>>::iterator it,
		std::vector<int> *f,
		int *prev_position,
		bool *deg,
		bool *reporter
		)
{
(*x) << (*P) << "$"; //concatenate unique symbol to P to initiaise string X



for (std::vector<std::string>::iterator j=(*it).begin(); j!=(*it).end(); j++){ //for each S_j in T[i]
	if ((*j) == "E"){
		(*epsilon) = true; //if S_j is empty string, set flag to true 
	} else {
		std::cout << (*j) << " ";
		(*x) << (*j) << "$"; //concatenate S_j and unique letter to string X
		if ((*j).length() >= (*P).length()){ //if P could occur in S_j
			//std::cout << "calling kmp" << std::endl;
			//if ( (*reporter)==false ){
			KMP( &(*P) , &(*j) , f, &(*deg), &(*report), &(*prev_position), &(*reporter));
			//}
		}
		(*Bprime).push_back((*x).str().length()-2); //in B': store ending pos of S_j in X
	}
}
(*X) = (*x).str(); //concatenation of P and all S_j, separated by unique chars
(*X).pop_back(); //remove unecessary unique letter at end pos of X
///std::cout << "String X: " << (*X) << std::endl; 
///printVector(Bprime); 
///std::cout << "check above indexes of below border table" << std::endl;
}

/*******************           parse input file       ************************/
// pre-processing
void parseInput(std::string *P,
		std::list<std::vector<std::string>> *T,
		std::string textfile,
		std::string patfile,
		std::ofstream *tempfile
		)
{
//text
std::string tline;
std::vector<std::string> tlines;
std::ifstream tFile(textfile);

if (tFile.is_open()){
	if (tFile.good()){
		while (getline(tFile, tline)){
			tlines.push_back(tline);
		}
	}
	tFile.close();
} else {
	std::cout << "Unable to open file." << std::endl;
}

std::string t = tlines[0];


std::vector<std::string> tempVector;
std::stringstream tempString;

for (int i=0; i<t.length(); i++){ //loop through text string
	if (t[i]=='{'){ 
		if (i!=0 && t[i-1] != '}'){
			tempVector.push_back(tempString.str());
			(*T).push_back(tempVector);
			tempVector.clear();
			tempString.str("");
			tempString.clear();
		}
	} else if (t[i]=='}'){ //if reach end of t[i]
		tempVector.push_back(tempString.str());
		(*T).push_back(tempVector);
		tempVector.clear();
		tempString.str("");
		tempString.clear();
	} else if (t[i]==','){ //if new s in t[i]
		tempVector.push_back(tempString.str());
		tempString.str("");
		tempString.clear();
	} else { //if t[i] == a 
		tempString << t[i]; //add a to string stream
	} //end_if
} //end_for


if (t[t.length()-1] != '}'){
		tempVector.push_back(tempString.str());
		(*T).push_back(tempVector);
		tempVector.clear();
		tempString.str("");
		tempString.clear();
}



//pattern
std::string pline;
std::vector<std::string> plines;
std::ifstream pFile(patfile);

if (pFile.is_open()){
	if (pFile.good()){
		while (getline(pFile, pline)){
			plines.push_back(pline);
		}
	}
	pFile.close();
} else {
	std::cout << "Unable to open file." << std::endl;
}

(*P) = plines[0];

(*tempfile).open("temporary.gism");
(*tempfile) << (*P);
(*tempfile).close();






/*
for (int i=0; i<t.length(); i++){ //loop through text string
	if (t[i]=='{'){ //if new pos in T
		tempVector.clear(); //clear vector to hold all s in t[i]
		tempString.str(""); //clear string stream to hold all a in s[j]
		tempString.clear();//clear string stream to hold all a in s[j]
	} else if (t[i]=='}'){ //if reach end of t[i]
		tempVector.push_back(tempString.str()); //add previous s to tempVector
		(*T).push_back(tempVector); //fill current pos in T with tempVector
	} else if (t[i]==','){ //if new s in t[i]
		tempVector.push_back(tempString.str()); //add previous s to tempVector
		tempString.str(""); //clear string stream to hold all a in s[j]
		tempString.clear(); //clear string stream to hold all a in s[j]
	} else { //if next a in s
		tempString << t[i]; //add a to string stream
	} //end_if
} //end_for
*/
}

/*******************           Update bit vector       ************************/
// O(Li_1)
void updateBitVector(std::vector<bool> *BV,
			std::vector<int> *Li_1,
			int m
			)
{
///std::cout << std::endl << "inside updateBV function" << std::endl;
for (std::vector<int>::iterator p = (*Li_1).begin(); p != (*Li_1).end(); p++){
	///std::cout << "for p=" << (*p) << " desired j is " << ((*p)+1) << std::endl;
	if ((*p)!=-1000){
		(*BV)[((*p)+1)] = true;
	}
}
}

/*******************           Print (T)ext and (P)attern       ************************/
// pre-processing
void printSeqs(std::list<std::vector<std::string>> *T,
		std::string *P
		)
{
std::cout << "string T:" << std::endl;
for (std::list<std::vector<std::string>>::iterator i=(*T).begin(); i!=(*T).end(); i++){
	std::cout << "pos " << std::distance((*T).begin(),i) << std::endl;
	for (std::vector<std::string>::iterator j=(*i).begin(); j!=(*i).end(); j++){
		std::cout << *j << " ";
	}
	std::cout << std::endl;
}
std::cout << std::endl << "string P:" << std::endl;
std::cout << (*P) << std::endl << std::endl;
}



/*********************                 Print L             **************************/
// post-processing
void printL(std::vector<std::vector<int>> *L)
{
for (std::vector<std::vector<int>>::iterator i = (*L).begin(); i != (*L).end(); i++)
{
	std::cout << std::endl << "L[" << std::distance((*L).begin(), i) << "]: ";
	std::vector<int> tempVec = *i;
	printVector(&tempVec);
}
}


/*********************          Compute lcp(suffix x, suffix y)          **************************/

// O(1)
int getlcp(int *suffx,
	int *suffy,
	std::vector<int> *iSA,
	std::vector<int> *LCP,
	sdsl::rmq_succinct_sct<> *rmq)
{
int i;
int j;
if ((*iSA)[(*suffx)] < (*iSA)[(*suffy)]){
	i = (*iSA)[(*suffx)];
	j = (*iSA)[(*suffy)];
} else {
	i = (*iSA)[(*suffy)];
	j = (*iSA)[(*suffx)];
}
//std::cout << "i = " << i << std::endl;
//std::cout << "j = " << j << std::endl;
auto min_idx = (*rmq)(i+1,j); 
int lcp = (*LCP)[min_idx];
return lcp;
}


/*********************          Compute lcp(suffix x, suffix y)          **************************/
/*
// O(1)
int getlcp(int *suffx,
	int *suffy,
	std::vector<int> *iSA,
	std::vector<int> *LCP
	)
{
int i;
int j;
int lcp;
if ((*iSA)[(*suffx)] < (*iSA)[(*suffy)]){
	i = (*iSA)[(*suffx)];
	j = (*iSA)[(*suffy)];
} else {
	i = (*iSA)[(*suffy)];
	j = (*iSA)[(*suffx)];
}

lcp = (*LCP)[i+1];
for (int k = i+2; k <= j; k++){
	if ((*LCP)[k] < lcp) lcp = (*LCP)[k];
}
return lcp;

}
*/
/*********************          Compute LCP array of string s         **************************/
// Kasai's algorithm
// O(n)

void computeLCParray(std::string *s,
				sdsl::csa_bitcompressed<> *SA,
				std::vector<int> *iSA,
				std::vector<int> *LCP
				)
{
(*s) = (*s).append("$");
int n = (*s).length(); 
int lcp = 0;
for (int i = 0; i < n; i++){
	if ((*iSA)[i] == n-1){
		lcp = 0;
		continue;
	}
	int j = (*SA)[(*iSA)[i]+1];
	while ( ( (i+lcp) < n ) && ( (j+lcp) < n ) && ( (*s)[i+lcp]==(*s)[j+lcp] ) ) lcp++;
	(*LCP)[(*iSA)[i]] = lcp;
	if (lcp > 0) lcp--;
}
std::vector<int>::iterator it = (*LCP).begin();
(*LCP).insert(it, 0);
(*LCP).pop_back();

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
void printSuffixArray(sdsl::csa_bitcompressed<> SA)
{
std::cout << std::endl << "Suffix array" << std::endl;
for (sdsl::csa_bitcompressed<>::iterator it = SA.begin(); it != SA.end(); it ++){
	std::cout << *it << " "; 
}
std::cout << std::endl << std::endl;
}

/*********************                     Compute B_p,s               **************************/
//O( B' log(B) )
void computeBps(std::vector<int> *Li, //stores all B_p,s (for all pos in T)
				std::vector<int> *B, //border table
				std::vector<int> *Bprime, //B'[j] = i s.t. i is ending pos of S_j in X
				std::string *P //pattern string
				)
{	
///std::cout << "computing Bp,s" << std::endl;
for (int i = 0; i != (*Bprime).size(); i++)
{
	int Bi = (*Bprime)[i];
	if ((*B)[Bi] != 0)
	{
		///std::cout << "looking at " << Bi << "th pos in B: " << (*B)[Bi] << std::endl;
		for(int x = (*B)[Bi]; x != 0; x = (*B)[x-1]){
			(*Li).push_back(x-1);
		}
	}
	else
	{
		(*Li).push_back(-1000);
	}
}
}


/*********************                     Compute border table               **************************/
// O(X)
// given parameters:
// X: string formed from concatenation of P and all S_j in T[i], separated by unique chars
// B: vector to hold border table
void computeBorderTable(std::string *X,
			std::vector<int> *B
			)
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
/*
std::cout << "border table" << std::endl;
for (std::vector<int>::iterator it = (*B).begin(); it != (*B).end(); it++){
	std::cout << *it << " ";
}
std::cout << std::endl << std::endl;
*/
}

/*********************                    KMP string matching algorithm              **************************/
// O(haystack)
//credit to http://www.sanfoundry.com/cpp-program-implement-kruth-morris-patt-algorithm-kmp/
//returns 1 if needle found in haystack, else returns 0

/*
void preKMP(std::string *pattern, std::vector<int> *f)
{
//std::cout << "doig nprekmp" << std::endl;
	int m = (*pattern).length();
	int j;
	(*f)[0] = -1;
	for (int i = 1; i < m; i++)
	{
		j = (*f)[i - 1];
		while (j >= 0)
		{
			if ((*pattern)[j] == (*pattern)[i - 1])
			{
				break;
			}
			else
			{
				j = (*f)[j];
			}
		}
		(*f)[i] = j + 1;
	}

}
*/
void preKMP(std::string *pattern, std::vector<int> *f)
{
int m = (*pattern).length();
int j = 0;
int i = 1;
(*f)[0] = 0;

while (i < m)
{
	if ( (*pattern)[i] ==  (*pattern)[j]){
		j++;
		(*f)[i] = j;
		i++;
	} else {
		if (j != 0){
			j = (*f)[j-1];
		} else {
			(*f)[i] = 0;
			i++;
		}
	}
}

}



/*
bool KMP(std::string *needle, std::string *haystack)
{
	int m = (*needle).length();
	int n = (*haystack).length();
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
        	else if ((*haystack)[i] == (*needle)[k])
		{
			k++;
			if (k==m) return 1;
			i++;
		}
		else
		{
			k = f[k];
		}
	}
	return 0;
}
*/

/*
int KMP(std::string *needle, std::string *haystack, std::vector<int> *f)
{
//std::cout << "using kmp" << std::endl;
	int m = (*needle).length();
	int n = (*haystack).length();
	//int f[m];
	//////////////preKMP(needle, f);
	int i = 0;
	int j = 0;
	while (i<n)
	{
		if (j==-1)
		{
			i++;
			j = 0;
		}
        	else if ((*haystack)[i] == (*needle)[j])
		{
			j++;
			if (j==m) return i;
			i++;
		}
		else
		{
			j = (*f)[j];
		}
	}
	return -1;
}
*/

void KMP(std::string *pattern, std::string *text, std::vector<int> *f, bool *deg, std::vector<int> *report, int *prev_position, bool *reporter)
{
//std::cout << "using kmp" << std::endl;
	int m = (*pattern).length();
	int n = (*text).length();

//std::cout << "pattern is " << (*pattern) << std::endl;
//std::cout << "S_j is " << (*text) << std::endl;

	int i = 0;
	int j = 0;
	while (i<n)
	{
		if ((*text)[i] == (*pattern)[j])
		{
			//std::cout << (*text)[i] << " = " << (*pattern)[j] << std::endl;
			i++;
			j++;
			
		}
		//std::cout << "j=" << j << "   " << "i=" << i << std::endl;
		if (j==m)
		{
			//std::cout << "found P in S_j ending at S_j[" << i-1 << "]" << std::endl;

				if ((*deg)==true) {
					(*report).push_back((*prev_position)+1);
					std::cout << "!!kmp reporting " << (*text) << std::endl;
					//(*reporter) = true;
					//break;
				} else {
					(*report).push_back((*prev_position) + i-1 + 1 );
					std::cout << "!!kmp reporting pos " << ( i-1 )  << " in " << (*text) << std::endl;
				}
			j = (*f)[j-1];
			//std::cout << "j=" << j << "   " << "i=" << i << std::endl;
		}
		if  ( (*text)[i] != (*pattern)[j] )
		{
			//std::cout << "moving along" << std::endl;
			if (j != 0){
				j = (*f)[j-1];
			} else {
				i++;
			}
			//std::cout << "j=" << j << "   " << "i=" << i << std::endl;
		}
	}
}













































/*
		//construct suffix array of X
		sdsl::csa_bitcompressed<> SA = computeSuffixArray(X);
		///printSuffixArray(SA);
		int size = SA.size();
		//construct inverse suffix array of X
		std::vector<int> iSA(size, 0);
		for (int r = 0; r != size; r ++) iSA[SA[r]] = r;
		///std::cout << "inverse suffix array" << std::endl; printVector(&iSA);
		//construct longest common prefix array of X + prepare for rmq
		std::vector<int> LCP(size, 0);
		computeLCParray(&X, &SA, &iSA, &LCP);
		sdsl::rmq_succinct_sct<> rmq; ////////////////////////////////////////////////////////////////////////////////////////////
		rmq = sdsl::rmq_succinct_sct<>(&LCP); ////////////////////////////////////////////////////////////////////////////////////////////
		///std::cout << "longest common prefix array" << std::endl; printVector(&LCP);
*/

//old step 3
/*
			for (int suffp = 1; suffp < P.length(); suffp++){ //for each suffix of P
				//////////////////////////////////////////////////////////////////int lcp = getlcp(&suffp, &suffs, &iSA, &LCP);
				int lcp = getlcp(&suffp, &suffs, &iSA, &LCP, &rmq); //lcp of S_j and suffix of P
				///std::cout << "\nlcp of suffixes " << suffp << " and " << suffs << " is " << lcp << std::endl;
				if (lcp > len) lcp = len;
				if (R[suffp] && lcp == (P.length()-suffp) ){
					///std::cout << "reporting " << i << std::endl;
					report.push_back(i); //report T[i]
				}
			} //end_for each suffix of P
*/



/*
// search for CTA inside ST(P)
std::string test_sj = "CTA";
std::vector<bool> test_occ( test_sj.length() , false ); //if all on then sj occurs in P
std::cout << "printing bool vector" << std::endl;
for (int i = 0; i < test_sj.length(); i++) std::cout << test_occ[i];
std::cout << std::endl;
//bitwise AND bool vector of same length, all true > sj found in p if all true
int test_int=0; //????????

sdsl::cst_sada<>::node_type v;
sdsl::cst_sada<>::char_type a;
//traversing entire ST
for (sdsl::cst_sada<>::const_iterator it = cst.begin(); it!=cst.end(); it++)
{
if(it.visit()==1) //if we have not traversed the subtree rooted at v
{

	v = *it; //node

	a = cst.edge(v,1);
	std::cout << "First letter on edge label from root to v: " << a << std::endl;

	if (a == test_sj[test_int]){ //if letter is letter of sj
		test_occ[test_int] = true; //set bit to on
		test_int++; //now look for next letter in sj
		for (int i = 1; i <= cst.depth(v); i++) //check rest of edge label
		{
			if ( cst.edge(v,i) == test_sj[test_int] ){
				test_occ[test_int] = true;
				test_int++;
			} else {
				continue;
			}
		}
		std::cout << std::endl;
	} else {
		std::cout << "failed to find - moving to next node" << std::endl;
	}

}
}

std::cout << "printing bool vector" << std::endl;
for (int i = 0; i < test_sj.length(); i++) std::cout << test_occ[i];
std::cout << std::endl;

if (test_occ) std::cout << "found!" << std::endl;


/********************************************************************************************
//Do stuff with STp
std::cout << "number of nodes in suffix tree " << cst.nodes() << std::endl << std::endl;

sdsl::cst_sada<>::size_type d;
sdsl::cst_sada<>::node_type v;
sdsl::cst_sada<>::size_type s;
sdsl::cst_sada<>::size_type sn;
bool l;
sdsl::cst_sada<>::size_type lb;
sdsl::cst_sada<>::size_type rb;
sdsl::cst_sada<>::size_type c;
sdsl::cst_sada<>::char_type a;

for (sdsl::cst_sada<>::const_iterator it = cst.begin(); it!=cst.end(); it++)
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
	
	a = cst.edge(v,1);
	std::cout << "First letter on edge label from root to v: " << a << std::endl;

	std::cout << "L(v) : ";
	for (int i = 1; i <= cst.depth(v); i++)
	{
		std::cout << cst.edge(v,i);
	}
	std::cout << std::endl;

	std::cout << std::endl;
}
}

auto testv = cst.select_child( cst.child( cst.root(),'A' ),2 );
std::cout << "TESTING TRAVERSAL ";
for (int i = 1; i <= cst.depth(testv); i++){
	std::cout << cst.edge(testv,i);
}
std::cout << std::endl;

std::cout << "I am a leaf: " << cst.is_leaf(testv) << std::endl;

auto testvv = cst.select_child( cst.child( testv,'T' ),1 );

std::cout << "TESTING TRAVERSAL ";
for (int i = 1; i <= cst.depth(testv); i++){
	std::cout << cst.edge(testv,i);
}
std::cout << std::endl;
**************************************************************************/


//old step 2
/*
				///std::cout << "length of S_j is less than P" << std::endl;
				for (int suffp = 1; suffp < P.length(); suffp++){ //for each suffix of P
					//////////////////////////////////////////////////////////////////int lcp = getlcp(&suffp, &suffs, &iSA, &LCP);
					int lcp = getlcp(&suffp, &suffs, &iSA, &LCP, &rmq); //lcp of S_j and suffix of P
					///std::cout << "\nlcp of suffixes " << suffp << " and " << suffs << " is " << lcp << std::endl;
					if (lcp >= len && E[suffp]){ //check if can extend prefix of P from L[i-1]
						///std::cout << "S_j occurs in P" << std::endl;
						///std::cout << "checking previous pos of P in L[i-1]" << std::endl;
						///std::cout << "can extend to pos ";
						int endpos = suffp+len-1; //can be extended to endpos in P
						///std::cout << endpos << " of P" << std::endl;
						if (PREV[endpos]){ //if endpos not in L[i]
							///std::cout << "not already added to L" << std::endl;
							PREV[endpos]=false; //do not allow to add endpos to L[i] again
							////L[i].push_back(endpos); //add endpos to L[i]
							Li.push_back(endpos); //add endpos to L[i]
						}
					} //end_if lcp>=len
				} //end_for each suffix of P
*/

/*
		int endpos;
		int cumulative_len = 0;
		for (int b = 0; b<Bprime.size(); b++){ //for all S_j in T[i]
			int len = Bprime[b] - P.length() - cumulative_len - b; //length of S_j
			int sj = Bprime[b]-len+1; //start pos of S_j in X
			cumulative_len += len; //update cum. length in preparation for next S_j
			if (len < P.length()) { //if S_j could occur in P
				std::string S_j = X.substr(sj,len);
				std::cout << std::endl << "following string can occur inside P " << S_j << std::endl;

				do {
					endpos = KMP(&S_j, &P);// return either -1 or pos in P where S_j ends
					std::cout << "result of KMP = " << endpos << std::endl;
					if ( endpos != -1 && E[endpos]==true ){
						if (PREV[endpos]){ //if endpos not in L[i]
							PREV[endpos]=false;
							Li.push_back(endpos); //add endpos to L[i]
							std::cout << "adding endpos to Li" << std::endl;
						}
					}
				} while ( (endpos != -1) || (endpos != (P.length()-1)) );

			} //end_if S_j could occur in P
		} //end_for all S_j in T[i]

*/







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
