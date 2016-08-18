#include <iostream> //cout,endl
#include <string> //string,
#include <vector> //vector,push_back
#include <fstream> //ifstream,is_open,good,getline,close
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

for (int i; i<t.length(); i++){
	std::cout << t[i] << " ";
	if(t[i]=="{"){
		tempVector.clear();
	}
}


//Construct Suffix Tree of pattern P


return 0;
}
