#ifndef PEDIGREE_H
#define PEDIGREE_H

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <utility>
#include <algorithm>
#include "CommonAlgorithms.cpp"

struct Individual
{
	std::string FamilyID;
	std::string IndividualID;
	std::string FatherID;
	std::string MotherID;
	int Sex;
	int AffectionStatus;
};

class Pedigree
{
	public:
	void PopulatePedigree(std::string&);
	std::map<std::string,std::vector<Individual> > Families;
	std::map<std::string,unsigned> NumberOfUnaffected;
	std::map<std::string,unsigned> NumberOfAffected;
};

void Pedigree::PopulatePedigree(std::string& PedigreeFilePath)
{
	std::ifstream PedigreeFile;
	
	if (PedigreeFilePath.empty()){ std::cerr << "No pedigree specified, exiting\n"; exit(0);}
	
	PedigreeFile.open(PedigreeFilePath.c_str());
	if (!PedigreeFile.is_open()) { std::cerr << "Pedigree file could not be opened, exiting\n"; exit(0);}
	
	std::string PedigreeLine;
	
	std::clog << "/\\ Designating pedigree...\n";
	
	while(getline(PedigreeFile,PedigreeLine))
	{
		std::vector<std::string> SplitPedigreeLine;
		tokenize(PedigreeLine,SplitPedigreeLine);
		Individual CurrentIndividual;
		CurrentIndividual.FamilyID = SplitPedigreeLine[0];
		CurrentIndividual.IndividualID = SplitPedigreeLine[1];
		CurrentIndividual.FatherID = SplitPedigreeLine[2];
		CurrentIndividual.MotherID = SplitPedigreeLine[3];
		CurrentIndividual.Sex = str2int(SplitPedigreeLine[4]);
		CurrentIndividual.AffectionStatus = str2int(SplitPedigreeLine[5]);
		
		if (CurrentIndividual.AffectionStatus == 1) NumberOfUnaffected[CurrentIndividual.FamilyID]++;
		if (CurrentIndividual.AffectionStatus == 2) NumberOfAffected[CurrentIndividual.FamilyID]++;
		
		Families[CurrentIndividual.FamilyID].push_back(CurrentIndividual);
	}
	PedigreeFile.close();
	for (std::map<std::string,std::vector<Individual> >::iterator CurrentFamily = Families.begin(); CurrentFamily != Families.end(); ++CurrentFamily)
	{
		for (std::vector<Individual>::iterator CurrentIndividual = CurrentFamily->second.begin(); CurrentIndividual != CurrentFamily->second.end(); ++CurrentIndividual)
		{
			std::clog << "/\\\tAdded individual " << CurrentIndividual->IndividualID << " (";
			if (CurrentIndividual->AffectionStatus == 2) std::clog << "affected) ";
			else std::clog << "unaffected) ";
			std::clog << "of family " << CurrentIndividual->FamilyID << "\n";
		}
	}
	std::clog << "/\\ Pedigree designation complete...\n";
}

#endif