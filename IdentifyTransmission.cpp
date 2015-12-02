#ifndef IDENTIFYTRANSMISSION_H
#define IDENTIFYTRANSMISSION_H

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

void IdentifyTransmission(TranscriptRNA& CurrentTranscript,Pedigree& CurrentPedigree)
{
	for (std::vector<Variant>::iterator CurrentVariant = CurrentTranscript.TranscriptVariants.begin(); CurrentVariant != CurrentTranscript.TranscriptVariants.end(); ++CurrentVariant)
	{
		for (std::map<std::string,std::vector<Individual> >::iterator CurrentFamily = CurrentPedigree.Families.begin(); CurrentFamily != CurrentPedigree.Families.end();++CurrentFamily)
		{
			for (std::vector<Individual>::iterator IndividualIterator = CurrentFamily->second.begin(); IndividualIterator != CurrentFamily->second.end();++IndividualIterator)
			{
				unsigned short IndividualAlleleCount = 0;
				if (CurrentVariant->Genotypes[IndividualIterator->IndividualID].first != CurrentVariant->Reference && CurrentVariant->Genotypes[IndividualIterator->IndividualID].first != ".") IndividualAlleleCount++;
				if (CurrentVariant->Genotypes[IndividualIterator->IndividualID].second != CurrentVariant->Reference && CurrentVariant->Genotypes[IndividualIterator->IndividualID].second != ".") IndividualAlleleCount++;
				
				if (IndividualAlleleCount == 0) continue;
					
				if (CurrentVariant->Genotypes.find(IndividualIterator->MotherID) == CurrentVariant->Genotypes.end() ) CurrentVariant->Genotypes[IndividualIterator->MotherID] = std::pair<std::string,std::string>(".",".");
				if (CurrentVariant->Genotypes.find(IndividualIterator->FatherID) == CurrentVariant->Genotypes.end() ) CurrentVariant->Genotypes[IndividualIterator->FatherID] = std::pair<std::string,std::string>(".",".");

				bool MotherHasAllele = false; bool FatherHasAllele = false;
				
				if (CurrentVariant->Genotypes[IndividualIterator->MotherID].first != CurrentVariant->Reference && CurrentVariant->Genotypes[IndividualIterator->MotherID].first != ".") MotherHasAllele = true;
				if (CurrentVariant->Genotypes[IndividualIterator->MotherID].second != CurrentVariant->Reference && CurrentVariant->Genotypes[IndividualIterator->MotherID].second != ".") MotherHasAllele = true;
				if (CurrentVariant->Genotypes[IndividualIterator->FatherID].first != CurrentVariant->Reference && CurrentVariant->Genotypes[IndividualIterator->FatherID].first != ".") FatherHasAllele = true;
				if (CurrentVariant->Genotypes[IndividualIterator->FatherID].second != CurrentVariant->Reference && CurrentVariant->Genotypes[IndividualIterator->FatherID].second != ".") FatherHasAllele = true;
				
				if (CurrentVariant->Chromosome != "X" && CurrentVariant->Chromosome != "Y")
				{
					if (IndividualAlleleCount == 1)
					{
						if (CurrentVariant->Genotypes[IndividualIterator->MotherID].first == "." || CurrentVariant->Genotypes[IndividualIterator->MotherID].second == "." || CurrentVariant->Genotypes[IndividualIterator->FatherID].first == "." || CurrentVariant->Genotypes[IndividualIterator->FatherID].second == ".")
						{
							CurrentVariant->InheritancePatterns[IndividualIterator->IndividualID].push_back(".");
							continue;
						}
						
						if (MotherHasAllele && FatherHasAllele) CurrentVariant->InheritancePatterns[IndividualIterator->IndividualID].push_back("EITHER");
						if (MotherHasAllele && !FatherHasAllele) CurrentVariant->InheritancePatterns[IndividualIterator->IndividualID].push_back(IndividualIterator->MotherID);
						if (!MotherHasAllele && FatherHasAllele) CurrentVariant->InheritancePatterns[IndividualIterator->IndividualID].push_back(IndividualIterator->FatherID);
						if (!MotherHasAllele && !FatherHasAllele) CurrentVariant->InheritancePatterns[IndividualIterator->IndividualID].push_back("DeNovo");
						continue;
					}
					if (IndividualAlleleCount == 2)
					{
						CurrentVariant->InheritancePatterns[IndividualIterator->IndividualID].push_back(IndividualIterator->MotherID);
						CurrentVariant->InheritancePatterns[IndividualIterator->IndividualID].push_back(IndividualIterator->FatherID);
					}
				}
				else
				{
					if (CurrentVariant->Chromosome == "Y")
					{
						if (IndividualIterator->Sex == 2) continue;
						if (CurrentVariant->Genotypes[IndividualIterator->FatherID].first == "." || CurrentVariant->Genotypes[IndividualIterator->FatherID].second == "." || FatherHasAllele)
						{
							CurrentVariant->InheritancePatterns[IndividualIterator->IndividualID].push_back(IndividualIterator->FatherID);
							continue;
						}
						else CurrentVariant->InheritancePatterns[IndividualIterator->IndividualID].push_back("DeNovo");
					}
					else
					{
						if (IndividualIterator->Sex == 1 || (IndividualAlleleCount == 1 && IndividualIterator->Sex == 2) )
						{
							if ((MotherHasAllele && FatherHasAllele) || CurrentVariant->Genotypes[IndividualIterator->MotherID].first == "." || CurrentVariant->Genotypes[IndividualIterator->MotherID].second == "." || CurrentVariant->Genotypes[IndividualIterator->FatherID].first == "." || CurrentVariant->Genotypes[IndividualIterator->FatherID].second == "." )
							{
								CurrentVariant->InheritancePatterns[IndividualIterator->IndividualID].push_back("EITHER");
								continue;
							}
							if (MotherHasAllele && !FatherHasAllele)
							{
								CurrentVariant->InheritancePatterns[IndividualIterator->IndividualID].push_back(IndividualIterator->MotherID);
								continue;
							}
							if (!MotherHasAllele && FatherHasAllele)
							{
								CurrentVariant->InheritancePatterns[IndividualIterator->IndividualID].push_back(IndividualIterator->FatherID);
								continue;
							}
							if (!MotherHasAllele && !FatherHasAllele && CurrentVariant->Genotypes[IndividualIterator->MotherID].first != "." && CurrentVariant->Genotypes[IndividualIterator->MotherID].second != "." && CurrentVariant->Genotypes[IndividualIterator->FatherID].first != "." && CurrentVariant->Genotypes[IndividualIterator->FatherID].second != "." )
							{
								 CurrentVariant->InheritancePatterns[IndividualIterator->IndividualID].push_back("DeNovo");
								 continue;
							}
						}
						else
						{
							if (IndividualAlleleCount == 2 && IndividualIterator->Sex == 2)
							{
								CurrentVariant->InheritancePatterns[IndividualIterator->IndividualID].push_back(IndividualIterator->MotherID);
								CurrentVariant->InheritancePatterns[IndividualIterator->IndividualID].push_back(IndividualIterator->FatherID);
							}
						}
					}
				}
			}
		}
	}
}

#endif