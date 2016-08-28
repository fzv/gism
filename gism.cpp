#include <iostream> //cout,endl
#include <string> //string,
#include <vector> //vector,push_back
#include <fstream> //ifstream,is_open,good,getline,close
#include <sstream> //
#include "sdsl/suffix_trees.hpp"

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

//Construct Suffix Tree of pattern P
std::string file = "pattern";
sdsl::cst_sada<> cst;
construct(cst, file, 1);
std::cout << "size " << cst.size() << std::endl;
for (sdsl::cst_sada<>::const_iterator it = cst.begin(); it!=cst.end(); it++){
	std::cout << cst.depth(*it) << "-[" << cst.lb(*it) << "," << cst.rb(*it) << "]" << std::endl;
}


return 0;
}
