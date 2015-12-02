#ifndef PARSEBED_H
#define PARSEBED_H

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
#include "Pedigree.cpp"

Variant ParseBED(std::string& BEDLine)
{
	Variant CurrentVariant;
	
	std::vector<std::string> SplitBEDLine;
	tokenize(BEDLine, SplitBEDLine);
	
	CurrentVariant.Chromosome = ProcessChromosomeString(SplitBEDLine[0]);
	
	CurrentVariant.Position = str2int(SplitBEDLine[1]);
	CurrentVariant.EndPosition = str2int(SplitBEDLine[2]);
	
	if (SplitBEDLine[3] == "NA") CurrentVariant.Novel = true;
	else CurrentVariant.Novel = false;
	
	CurrentVariant.PassedFilters = true;
	CurrentVariant.PassedFrequencyCutoffs = true;
	CurrentVariant.IsSecondary = true;
	
	return(CurrentVariant);
}

#endif