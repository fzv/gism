#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream> 
#include "sdsl/suffix_trees.hpp"
#include "sdsl/suffix_arrays.hpp"
#include "sdsl/lcp.hpp"
#include "sdsl/util.hpp"
#include <iterator>
#include "sdsl/rmq_support.hpp"
#include <algorithm>
#include <ctime> 

/************************************************************************************/
/******************************* FUNCTION DECLARATIONS ******************************/
/************************************************************************************/

void parseInput(
std::string *P
, std::list<std::vector<std::string>> *T
, std::string textfile
, std::string patfile
, std::ofstream *tempfile
);
void prepareX(
std::vector<int> *Bprime
, bool *epsilon
, std::string *P
, std::vector<int> *report
, int *i
, std::string *X
, std::list<std::vector<std::string>>::iterator it
, std::vector<int> *f
, int *prev_position
, bool *deg
, bool *reporter
//, std::ofstream *logfile
);
void computeBorderTable(
std::string *X
, std::vector<int> *B
);
void preKMP(
std::string *pattern
, std::vector<int> *f
);
void KMP(
std::string *needle
, std::string *haystack
, std::vector<int> *f
, bool *deg
, std::vector<int> *report
, int *prev_position
, bool *reporter
//, std::ofstream *logfile
);
void computeBps(
std::vector<int> *Li
, std::vector<int> *B
, std::vector<int> *Bprime
, std::string *P
);
sdsl::csa_bitcompressed<> computeSuffixArray(
std::string s
);
void updateBitVector(
std::vector<bool> *BV
, std::vector<int> *Li_1
, int m
);
void reporting(
std::vector<int> *vector
, std::ofstream *outfile
);
void extend(
std::string *sj
, sdsl::cst_sada<> *stp
, sdsl::csa_bitcompressed<> *sap
, std::vector<int> *Li
, std::vector<bool> *E
, std::vector<bool> *PREV
, int *len
);
void extendToEnd(
std::string *sj
, sdsl::cst_sada<> *stp
, sdsl::csa_bitcompressed<> *sap
, std::vector<int> *Li
, std::vector<bool> *R
, int *len
, std::vector<int> *report
, int *prev_position
, std::string *P
, bool *deg
, bool *reporter
//, std::ofstream *logfile
);

// testing purposes only /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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

/* Begin logging */
/*
std::ofstream logfile; 
logfile.open("gism.log");
logfile << "log" << std::endl;
*/
/* Declare & Assign Seq Variables */

// pattern
std::string P;
std::ofstream tempfile; //used to construct suffix tree of P
std::string patfile = argv[2];

// text
std::list<std::vector<std::string>> T; //stores all of T
std::string textfile = argv[1];

// parse input
parseInput(&P, &T, textfile, patfile, &tempfile);

/* Begin Clock */

std::clock_t start_time = clock();

/* Pre-process Pattern */

// for KMP
std::vector<int> f(P.length(),0);
preKMP(&P, &f);

// suffix tree of P
sdsl::cst_sada<> stp; //documentation says construction is slow but fast operations, compared to other = vice versa
construct(stp, "temporary.gism", 1); //build suffix tree using temp pattern file
std::remove("temporary.gism");

// suffix array of pattern P
sdsl::csa_bitcompressed<> sap = computeSuffixArray(P);

/* GISM */

std::vector<int> report; //stores all pos in T where occ of P ends
std::vector<int> Li; //stores ending positions of prefixes of P found in T[i]
std::vector<int> Li_1; //kept in case L_{i-1} needs to be copied to L_i if epsilon in T[i]
int prev_position = -1; //used to calculate current pos (different to current segment T[i])

/* Loop Through Each T[i] */
for (std::list<std::vector<std::string>>::iterator it=T.begin(); it!=T.end(); it++){

	/* declare segment number */
	int i = std::distance(T.begin(),it);

	/* is this segment a degenerate pos? */
	bool deg = false;
	if ( (*it).size() > 1) deg = true;

	/* if yes, has it already been reported? */
	bool reporter = false;

	/* prepare all S_j in T[i] */
	std::vector<int> Bprime; //B'[j] = i s.t. i is ending pos of S_j in X
	bool epsilon = false; //flag empty string in T[i]
	std::vector<int> B; //border table
	std::string X; //as defined in paper








	//X = "ACGATC$A$ACGATCAAA$AAAACGATC$ACG$ATC";
	//Bprime = {7,17,27,31,35};
	prepareX(&Bprime, &epsilon, &P, &report, &i, &X, it, &f, &prev_position, &deg, &reporter /*, &logfile */ );
	//std::cout << X << std::endl;

for (std::vector<std::string>::iterator Ti = (*it).begin(); Ti != (*it).end(); Ti++){
	int len = (*Ti).length();
	if (len >= P.length() && (*Ti)!="E"){
	KMP( &P , &(*Ti) , &f , &deg , &report , &prev_position , &reporter /* , &(*logfile) */ );
	}
}








///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/***
if (i < 800){
std::cout << prev_position+1 << " " << "epsilon?" << epsilon << " " << X << std::endl;
}
***/

	if (reporter!=true){

	/* begin GISM algorithm */
	if (it==T.begin()) { //only for T[0] do:

		/* STEP 1: FIND PREFIXES OF P */
		computeBorderTable(&X, &B);
		computeBps(&Li, &B, &Bprime, &P);

	} else {

		/* Reset L_i and L_{i-1} */
		Li_1 = Li;
		Li.clear();

		/* STEP 1: FIND PREFIXES OF P */
		computeBorderTable(&X, &B);
		computeBps(&Li, &B, &Bprime, &P);

		/* initialise bitvector (E)xtend to aid extension of prefixes of P */
		std::vector<bool> E(P.length(),false);
		updateBitVector(&E, &Li_1, P.length());

		/* initialise bitvector (PREV)iously extended to aid extension of prefixes of P */
		std::vector<bool> PREV(P.length(),true);

		/* STEP 2: EXTEND PREFIXES OF P */

		for (std::vector<std::string>::iterator Ti = (*it).begin(); Ti != (*it).end(); Ti++){

			int len = (*Ti).length();

			if (len < P.length()) { //if S_j could occur in P

				extend(&(*Ti), &stp, &sap, &Li, &E, &PREV, &len); //extension function

			} //END_IF(S_j could occur in P)

		} //END_FOR(all S_j in T[i])

		/* STEP 3: EXTEND PREFIXES OF P TO END OF P*/

		std::vector<bool> R(P.length(),false); //bitvector (R)eport to aid reporting of occurences of P in T
		updateBitVector(&R, &Li_1, P.length());


		std::vector<std::string>::iterator Ti = (*it).begin();
		do{
			int len = (*Ti).length();
			extendToEnd(&(*Ti), &stp, &sap, &Li, &R, &len, &report, &prev_position, &P, &deg, &reporter /* , &logfile */ );
			Ti++;
		} while ( reporter!=true && Ti != (*it).end() );


//		for (std::vector<std::string>::iterator Ti = (*it).begin(); Ti != (*it).end(); Ti++){
//
//			int len = (*Ti).length();
//			extendToEnd(&(*Ti), &stp, &sap, &Li, &R, &len, &report, &prev_position, &P, &deg, &reporter /* , &logfile */ );
//
//		} //END_FOR(all S_j in T[i])


		if (epsilon==true) Li.insert(Li.end(), Li_1.begin(), Li_1.end());

	} //END_IF(T[0])

	} //END_IF(reporter!=true)

	/* update position in T */
	if ((*it).size() > 1) {
		prev_position++;
	} else {
		prev_position += (*it)[0].size();
	}

} //END_FOR(GISM)


/* Stop Clock */

std::clock_t end_time = clock();
double elapsed_secs = double(end_time - start_time) / CLOCKS_PER_SEC;

/* Report */

// create output file
std::ofstream outfile; 
std::string output_file = argv[3];
outfile.open(output_file);

// record time
outfile << elapsed_secs << std::endl;

// record number of ending positions reported
outfile << "#" << report.size() << std::endl;

// record all ending positions
reporting(&report, &outfile);

//close files
outfile << std::endl;
outfile.close();

//end logging
//logfile.close();

return 0;

} //END_MAIN

/**********************************************************************************************/
/************************************ FUNCTION DEFINITIONS ************************************/
/**********************************************************************************************/
//
/************************************************************************************************
function: compute suffix array of input string
time complexity:
space complexity:
*/
sdsl::csa_bitcompressed<> computeSuffixArray
(
std::string s
)
{
sdsl::csa_bitcompressed<> SA;
construct_im(SA, s, 1);
return SA;
}
/************************************************************************************************
function: extend prefixes of P to end of P (suffixes of P) using letters in T (prefixes of S_j)
time complexity:
space complexity:
*/
void extendToEnd
( //PARAMS
  std::string *sj
, sdsl::cst_sada<> *stp
, sdsl::csa_bitcompressed<> *sap
, std::vector<int> *Li
, std::vector<bool> *R
, int *len
, std::vector<int> *report
, int *prev_position
, std::string *P
, bool *deg
, bool *reporter
//, std::ofstream *logfile
) //END_PARAMS
{ //FUNCTION
std::string s = (*sj);
int occ = -1;
sdsl::cst_sada<>::size_type lb;
sdsl::cst_sada<>::size_type rb;
sdsl::cst_sada<>::node_type v = (*stp).root();
auto it = s.begin();
for (uint64_t char_pos = 0; it != s.end(); ++it){ //for each letter in S_j???????
	if ( forward_search( (*stp), v, it-s.begin(), *it, char_pos) > 0 ){ //if letter occurs in P
		occ = it-s.begin();
		lb = (*stp).lb(v);
		rb = (*stp).rb(v);
		if (occ > -1){ //prefix of S_j found in P, but where in P? Must be suffix!
			for (int i = lb; i < rb+1; i++){
				int startpos = (*sap)[i];
				int endpos = startpos+occ; //can be extended to P[endpos]
				if ( (*R)[startpos] && endpos == ( (*P).length()-1 ) ){
					if ( (*deg)==true  && (*reporter)==false  ) {
						(*report).push_back((*prev_position)+1);
						(*reporter)==true;
					} else {
						(*report).push_back((*prev_position)+occ+1);
					} //END_IF
				} //END_IF
			} //END_FOR
		} //END_IF
	} else { //if letter does not occur in P
		break;
	} //END_IF
} //END_FOR
} //END_FUNCTION(extendToEnd)
/************************************************************************************************
function:
time complexity:
space complexity:
*/
void extend
( //PARAMS
  std::string *sj
, sdsl::cst_sada<> *stp
, sdsl::csa_bitcompressed<> *sap
, std::vector<int> *Li
, std::vector<bool> *E
, std::vector<bool> *PREV
, int *len
) //END_PARAMS
{ //FUNCTION
std::string s = (*sj);
int occ = -1;
sdsl::cst_sada<>::size_type lb;
sdsl::cst_sada<>::size_type rb;
sdsl::cst_sada<>::node_type v = (*stp).root();
auto it = s.begin();
for (uint64_t char_pos = 0; it != s.end(); ++it){
	if ( forward_search( (*stp), v, it-s.begin(), *it, char_pos) > 0 ){ //if S_j found in P
		occ = it-s.begin();
		lb = (*stp).lb(v);
		rb = (*stp).rb(v);
	} else {
		break;
	} //END_IF
} //END_FOR
if (occ == (*len)-1){ // whole sj found in P, but where?
	for (int i = lb; i < rb+1; i++){
		int startpos = (*sap)[i];
		int endpos = startpos+occ;
		if ((*E)[startpos] && (*PREV)[endpos]){
			(*PREV)[endpos]=false;
			(*Li).push_back(endpos);
		} //END_IF
	} //END_FOR
} //END_IF
} //END_FUNCTION(extend)
/************************************************************************************************
function: print unique values in vector report
time complexity: linear in size of vector report
space complexity: none extra
*/
void reporting
( //PARAMS
  std::vector<int> *vector
, std::ofstream *outfile
) //END_PARAMS
{ //FUNCTION
if ( (*vector).size()==0 ){
	std::cout << "Nothing to report." << std::endl;
} else {
	(*outfile) << (*vector)[0];
	for (int i=1; i<(*vector).size(); i++){
		if((*vector)[i] != (*vector)[i-1]) (*outfile) << std::endl << (*vector)[i];
	} //END_FOR
} //END_IF
} //END_FUNCTION
// FOR COMPARSON WITH RITU ONLY:
/*
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
*/
/************************************************************************************************
function to
time complexity:
space complexity:
*/
void prepareX
( //PARAMS
  std::vector<int> *Bprime
, bool *epsilon
, std::string *P
, std::vector<int> *report
, int *i
, std::string *X
, std::list<std::vector<std::string>>::iterator it
, std::vector<int> *f
, int *prev_position
, bool *deg
, bool *reporter
//, std::ofstream *logfile
) //END_PARAMS
{ //FUNCTION

std::stringstream x;
x << (*P) << "$";
for (std::vector<std::string>::iterator j=(*it).begin(); j!=(*it).end(); j++){ //for each S_j in T[i]
	if ((*j) == "E"){
		(*epsilon) = true;
	} else {
		x << (*j) << "$";
		(*Bprime).push_back(x.str().length()-2);
	}
}
(*X) = x.str(); 
(*X).pop_back();





/*
std::stringstream x;
x << (*P) << "$"; //concatenate unique symbol to P to initiaise string X
for (std::vector<std::string>::iterator j=(*it).begin(); j!=(*it).end(); j++){ //for each S_j in T[i]
	if ((*j) == "E"){
		(*epsilon) = true; //if S_j is empty string, set flag to true 
	} else {
		x << (*j) << "$"; //concatenate S_j and unique letter to string X
		if ( (*j).length() >= (*P).length() ){ //if P could occur in S_j
			//KMP( &(*P) , &(*j) , f , &(*deg) , &(*report) , &(*prev_position) , &(*reporter) , &(*logfile) );
			KMP( &(*P) , &(*j) , f , &(*deg) , &(*report) , &(*prev_position) , &(*reporter) );
		}
		(*Bprime).push_back(x.str().length()-2); //in B': store ending pos of S_j in X
	} //END_IF
} //END_FOR
(*X) = x.str(); //concatenation of P and all S_j, separated by unique chars
(*X).pop_back(); //remove unecessary unique letter at end pos of X
*/




} //END_FUNCTION
/************************************************************************************************
function:
time complexity:
space complexity:
*/
void parseInput
( //PARAMS
std::string *P,
std::list<std::vector<std::string>> *T,
std::string textfile,
std::string patfile,
std::ofstream *tempfile
) //END_PARAMS
{ //FUNCTION
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
	} //END_IF
} //END_FOR
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
} //END_FUNCTION
/************************************************************************************************
function to
*/
void updateBitVector
( //PARAMS
  std::vector<bool> *BV
, std::vector<int> *Li_1
, int m
) //END_PARAMS
{ //FUNCTION
for (std::vector<int>::iterator p = (*Li_1).begin(); p != (*Li_1).end(); p++){
	if ((*p)!=-1000){
		(*BV)[((*p)+1)] = true;
	}
}
} //END_FUNCTION
/************************************************************************************************
function to
*/
void printSeqs
( //PARAMS
  std::list<std::vector<std::string>> *T
, std::string *P
) //END_PARAMS
{ //FUNCTION
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
} //END_FUNCTION
/************************************************************************************************
function to
*/
void computeBps
( //PARAMS
  std::vector<int> *Li //stores all B_p,s (for all pos in T)
, std::vector<int> *B //border table
, std::vector<int> *Bprime //B'[j] = i s.t. i is ending pos of S_j in X
, std::string *P //pattern string
) //END_PARAMS
{ //FUNCTION
for (int i = 0; i != (*Bprime).size(); i++)
{
	int Bi = (*Bprime)[i];
	if ((*B)[Bi] != 0)
	{
		for(int x = (*B)[Bi]; x != 0; x = (*B)[x-1]){
			(*Li).push_back(x-1);
		}
	} else {
		(*Li).push_back(-1000);
	} //END_IF
} //END_FOR
} //END_FUNCTION(computeBps)
/************************************************************************************************
function to
*/
void computeBorderTable
( //PARAMS
  std::string *X
, std::vector<int> *B
) //END_PARAMS
{ //FUNCTION
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
		} //END_IF
	} //END_WHILE
	p++;
} //END_FOR
(*B)[m-1] = p;
} //END_FUNCTION(computeBorderTable)
/************************************************************************************************
function to
*/
void preKMP
( //PARAMS
  std::string *pattern
, std::vector<int> *f
) //PARAMS
{ //FUNCTION
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
		} //END_IF
	} //END_IF
} //END_WHILE
} //END_FUNCTION(preKMP)
/************************************************************************************************
function:
time complexity:
space complexity:
source: https://en.wikibooks.org/wiki/Algorithm_Implementation/String_searching/Knuth-Morris-Pratt_pattern_matcher#Python
*/



void KMP
( //PARAMS
  std::string *pattern
, std::string *text
, std::vector<int> *f
, bool *deg
, std::vector<int> *report
, int *prev_position
, bool *reporter
//, std::ofstream *logfile
) //END_PARAMS
{ //FUNCTION
int m = (*pattern).length();
int n = (*text).length();
int i = 0;
int j = 0;
while (i<n)
{
	if ((*text)[i] == (*pattern)[j])
	{
		i++;
		j++;
	} //END_IF
	if (j==m)
	{
		if ((*deg)==true) {
			(*report).push_back((*prev_position)+1);
			(*reporter) = true;
			break;
		} else {
			(*report).push_back((*prev_position) + i-1 + 1 );
		} //END_IF
		j = (*f)[j-1];
	} //END_IF
	if  ( (*text)[i] != (*pattern)[j] )
	{
		if (j != 0){
			j = (*f)[j-1];
		} else {
			i++;
		} //END_IF
	} //END_IF
} //END_WHILE
} //END_FUNCTION(kmp)











//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
// NOTES
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


// FOR TESTING PURPOSES ONLY


/*********************          Print any int vector         **************************/
void printVector(std::vector<int> *vector)
{
for (std::vector<int>::iterator it = (*vector).begin(); it != (*vector).end(); it ++){
	std::cout << *it << " "; 
}
std::cout << std::endl << std::endl;
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



/*
while (i < n)
{
if (j == -1){
	i++;
	j = 0;
} else if ( (*text)[i] == (*pattern)[j] ) {
	j++;
		if ((*deg)==true) {
			(*report).push_back((*prev_position)+1);
			(*reporter) = true;
			break;
		} else {
			(*report).push_back((*prev_position) + i-1 + 1 );
		} //END_IF
	i++;
} else {
	j = (*f)[j];
}
}
*/

/*
(*f)[0] = -1;
for (int i = 1; i < m; i++){
	j = (*f)[i-1];
	while ( j >= 0 ){
		if ( (*pattern)[j] == (*pattern)[i-1] ){
			break;
		} else {
			j = (*f)[j];
		}
	}
	(*f)[i] = j+1;
}
*/



















