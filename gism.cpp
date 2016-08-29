#include <iostream> //cout,endl
#include <string> //string,
#include <vector> //vector,push_back
#include <fstream> //ifstream,is_open,good,getline,close
#include <sstream> //
#include "sdsl/suffix_trees.hpp"
#include "sdsl/util.hpp"
#include <iterator>

std::vector<int> computeBorderTable(std::string X, std::vector<int> B);
std::vector<int> computeBorder(std::string temp, std::vector<int> B);

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
/// if len(s) >= m ???
std::stringstream x;
std::string X;
std::string alpha = "mnopqrsvwxyz";
int unique = 0;
std::vector<int> B;

for (std::list<std::vector<std::string>>::iterator i=T.begin(); i!=T.end(); i++){
	tempVector = *i;
	for (std::vector<std::string>::iterator j=tempVector.begin(); j!=tempVector.end(); j++){
		x << *j << alpha[unique];
		unique++;
	}
	X = x.str();
	X.pop_back();
	//string X is ready
	std::cout << X << " : " << std::endl;
	computeBorderTable(X, B);
	//clean up
	x.str("");
	x.clear();
	unique = 0;
	std::cout << std::endl;
}

return 0;
}

std::vector<int> computeBorderTable(std::string X, std::vector<int> B)
{
	std::string temp;
	B.push_back(0);
	for (int i = 2; i < X.length(); i++){
		temp = X.substr(0,i);
		std::cout << temp << std::endl;
		//compute border length of string temp
		computeBorder(temp, B);

	}
	temp = X;
	std::cout << temp << std::endl;
	computeBorder(temp, B);

	//print border table
	for (std::vector<int>::iterator it = B.begin(); it != B.end(); it++){
		std::cout << *it << " ";
	}
	std::cout << std::endl << std::endl;

	return B;
}

std::vector<int> computeBorder(std::string temp, std::vector<int> B)
{
	int len = temp.length();
	std::cout << "Length of temp is " << len << std::endl;
	int borderlen = 0;
	std::cout << "Borderlen initialised to " << borderlen << std::endl;
	for (int j = 0; j < len; j++){
		std::cout << "Looping through temp string" << std::endl;
		char abc = temp[j];
		std::cout << "First letter in temp is " << abc << " at pos " << j << std::endl;
		char xyz = temp[len-j-1];
		std::cout << "Last letter in temp is " << xyz << " at pos " << len-j-1 << std::endl;
		if (abc==xyz){
			std::cout << "First and last letter match" << std::endl;
			borderlen++;
		} else if (xyz != 'A' | xyz != 'G' | xyz != 'C' | xyz != 'T'){
			std::cout << "Last letter is not ACGT" << std::endl;
			break;
		} else {
			std::cout << "Letters do not match" << std::endl;
			break;
		}
	}
	std::cout << "Borderlen = " << borderlen << std::endl <<std::endl;
	B.push_back(borderlen);
	
	return B;
}
