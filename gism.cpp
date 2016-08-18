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
			//std::cout << line << std::endl;
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
	std::cout << "just read " << t[i] << std::endl; //print current symbol
	if (t[i]=='{'){ //if new pos in T
		tempVector.clear(); //create vector to hold letters in current pos in T
		tempString.str("");
		tempString.clear();
		std::cout << "printing tempVector\n";
		for (auto i=0; i< tempVector.size(); i++) std::cout << tempVector[i] << std::endl;
		std::cout << "\nprinting tempString\n";
		std::cout << tempString.str() << std::endl;
	} else if (t[i]=='}'){ //if reach end of current pos in T
		tempVector.push_back(tempString.str()); //add current s to tempVector
		T.push_back(tempVector); //fill current pos in T with tempVector
		std::cout << "printing tempVector\n";
		for (auto i=0; i< tempVector.size(); i++) std::cout << tempVector[i] << std::endl;
		std::cout << "\nprinting tempString\n";
		std::cout << tempString.str() << std::endl;
	} else if (t[i]==','){
		tempVector.push_back(tempString.str()); //add current s to tempVector
		tempString.str("");
		tempString.clear();
		std::cout << "printing tempString\n";
		std::cout << tempString.str() << std::endl;
		std::cout << "printing tempVector\n";
		for (auto i=0; i< tempVector.size(); i++) std::cout << tempVector[i] << std::endl;
		std::cout << std::endl;
	} else {
		std::cout << "printing tempString\n";
		tempString << t[i];
		std::cout << tempString.str() << std::endl;
	}
	std::cout << std::endl << std::endl;
}
tempVector.clear();
tempString.str("");
tempString.clear();

std::cout << std::endl << std::endl << "trying to print T properly" << std::endl << std::endl;
for (auto i=0; i<T.size(); i++){
	tempVector = T[i];
}



//Construct Suffix Tree of pattern P


return 0;
}
