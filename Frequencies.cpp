#ifndef FREQUENCIES_H
#define FREQUENCIES_H

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
#include "Settings.cpp"

class AlleleFrequency
{
	public:
	AlleleFrequency(std::vector<std::string>&);
	AlleleFrequency();
	float Frequency;
	std::string Reference;
	std::string Alternate;
};

AlleleFrequency::AlleleFrequency(std::vector<std::string>& Data)
{
	Reference = Data[3];
	Alternate = Data[4];
	Frequency = ::atof(Data[5].c_str());
}

AlleleFrequency::AlleleFrequency()
{
	Reference = ".";
	Alternate = ".";
	Frequency = 0;
}

struct ReferenceAlleleFrequency
{
    std::string Reference;
    std::string Alternate;
	double Frequency;
    ReferenceAlleleFrequency(): Reference(""), Alternate(""), Frequency(0)
    { }
    ReferenceAlleleFrequency(std::string Reference_, std::string Alternate_, double Frequency_)
    : Reference(Reference_)
    , Alternate(Alternate_)
    , Frequency(Frequency_)
    { }
};

void FilterVariant(Variant& CurrentVariant, ProgramOptions& Settings)
{
	// PASS in the VCF filter column?
	
	CurrentVariant.PassedFilters = true;
	for (auto CurrentFilter = CurrentVariant.Filter.begin(); CurrentFilter != CurrentVariant.Filter.end(); ++CurrentFilter)
	{
		if (*CurrentFilter == "PASS" || *CurrentFilter == ".")
		{}
		else
		{
//			CurrentVariant.PassedFilters = false;
		}
	}
	
	// Variant frequency above threshold?
	
//	if (CurrentVariant.Quality < Settings.QualityCutoff) CurrentVariant.PassedFrequencyCutoffs = false;
	
	if (CurrentVariant.Annotations.find("ArtifactAC") != CurrentVariant.Annotations.end())
	{
		std::vector<std::string> Frequencies;
		tokenize(CurrentVariant.Annotations["ArtifactAC"],Frequencies,",");
		for (auto CurrentFrequency = Frequencies.begin(); CurrentFrequency != Frequencies.end(); ++CurrentFrequency)
		{
			if (*CurrentFrequency == ".")
			{
				continue;
			}
//			std::cerr <<*CurrentFrequency << "\n";
			if (std::stoi(*CurrentFrequency) > Settings.ArtifactCutoff)
			{
				CurrentVariant.PassedFrequencyCutoffs = false;
			}
		}
	}
	
	if (CurrentVariant.Annotations.find("ExAC_AF") != CurrentVariant.Annotations.end())
	{
		std::vector<std::string> Frequencies;
		tokenize(CurrentVariant.Annotations["ExAC_AF"],Frequencies,",");
		
		for (auto CurrentFrequency = Frequencies.begin(); CurrentFrequency != Frequencies.end(); ++CurrentFrequency)
		{
			if (CurrentVariant.PassedFrequencyCutoffs == false) continue;
			
			if (*CurrentFrequency == ".")
			{
				CurrentVariant.PassedFrequencyCutoffs = true;
				break;
			}
			if (std::stof(*CurrentFrequency) > Settings.FrequenciesCutoff)
			{
				CurrentVariant.PassedFrequencyCutoffs = false;
			}
		}
	}
	
	
	// Any other filtering to do?
	
}
	



#endif