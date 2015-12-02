#ifndef IMPRINTTEST_H
#define IMPRINTTEST_H

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

void TestGeneWithImprintedGeneModel(TranscriptRNA& CurrentTranscript, Pedigree& CurrentPedigree, std::vector<std::string>& AnnotationsToInclude)
{
	std::vector<std::string> ImprintModels;
	ImprintModels.push_back("Paternal");
	ImprintModels.push_back("Materal");
	
	for (std::map<std::string,std::vector<Individual> >::iterator CurrentFamily = CurrentPedigree.Families.begin(); CurrentFamily != CurrentPedigree.Families.end();++CurrentFamily)
	{
		for (std::vector<Variant>::iterator CurrentVariant = CurrentTranscript.TranscriptVariants.begin(); CurrentVariant != CurrentTranscript.TranscriptVariants.end(); ++CurrentVariant)
		{
			// Get annotations
			bool SufficientSeverity = false;
			
			for (std::vector<SnpEffDB>::iterator CurrentAnnotationSeverity = CurrentVariant->AnnotationDB.begin(); CurrentAnnotationSeverity != CurrentVariant->AnnotationDB.end(); ++CurrentAnnotationSeverity)
			{
				if (( CurrentAnnotationSeverity->EffectImpact == "HIGH" || CurrentAnnotationSeverity->EffectImpact == "MODERATE" || CurrentAnnotationSeverity->EffectImpact == "LOW") && CurrentTranscript.Name == CurrentAnnotationSeverity->TranscriptID)
				{
					SufficientSeverity = true;
					continue;
				}
			}
			if (CurrentVariant->Annotations["dbNSFP_clinvar_clnsig"] != "." && !CurrentVariant->Annotations["dbNSFP_clinvar_clnsig"].empty()) SufficientSeverity = true;
			
			// Discard variants that are filtered or not impactful
			if ((!CurrentVariant->PassedFilters || !CurrentVariant->PassedFrequencyCutoffs || !SufficientSeverity) && !CurrentVariant->IsSecondary) continue;
			
			bool VariationFitsModel = true;
			std::string ModelType;
			
			
			for ( auto CurrentImprintModel = ImprintModels.begin(); CurrentImprintModel != ImprintModels.end(); ++CurrentImprintModel)
			{
				short TotalAlleles = 0;
				
				for (std::vector<Individual>::iterator IndividualIterator = CurrentFamily->second.begin(); IndividualIterator != CurrentFamily->second.end();++IndividualIterator)
				{
					if (CurrentVariant->Genotypes.find(IndividualIterator->IndividualID) == CurrentVariant->Genotypes.end()) CurrentVariant->Genotypes[IndividualIterator->IndividualID] = std::pair<std::string,std::string>(".",".");
					
				}
				for (std::vector<Individual>::iterator IndividualIterator = CurrentFamily->second.begin(); IndividualIterator != CurrentFamily->second.end();++IndividualIterator)
				{
					//Ignore comparisons with missing data
					if (CurrentVariant->Genotypes[IndividualIterator->IndividualID].first == "." || CurrentVariant->Genotypes[IndividualIterator->IndividualID].second == "." || CurrentVariant->IsSecondary) continue;
					
					short AlleleCount = 0;
					if (CurrentVariant->Genotypes[IndividualIterator->IndividualID].first != CurrentVariant->Reference && CurrentVariant->Genotypes[IndividualIterator->IndividualID].first != ".") AlleleCount++;
					if (CurrentVariant->Genotypes[IndividualIterator->IndividualID].second != CurrentVariant->Reference && CurrentVariant->Genotypes[IndividualIterator->IndividualID].second != ".") AlleleCount++;
					TotalAlleles = TotalAlleles + AlleleCount;
					
					if (IndividualIterator->AffectionStatus == 2)
					{
						if (AlleleCount == 0)
						{
							VariationFitsModel = false;
							break;
						}
						
						if (AlleleCount == 1)
						{
							for (std::vector<std::string>::iterator CurrentInheritance = CurrentVariant->InheritancePatterns[IndividualIterator->IndividualID].begin(); CurrentInheritance != CurrentVariant->InheritancePatterns[IndividualIterator->IndividualID].end(); ++CurrentInheritance)
							{
								if (*CurrentInheritance == "DeNovo") ModelType = "DeNovo";
							}
							continue;
						}
						else
						{
							VariationFitsModel = false;
							break;
						}
					}
					
					if (IndividualIterator->AffectionStatus == 1)
					{
						if (AlleleCount != 0) VariationFitsModel = false;
						else continue;
					}
				}
				
				// Check if all data is missing...
				if (TotalAlleles == 0) VariationFitsModel = false;
				
				
				if (VariationFitsModel || CurrentVariant->IsSecondary)
				{
					std::cout << CurrentTranscript.GeneName << "\t" << CurrentTranscript.Name << "\t" << CurrentFamily->first << "\t" << ModelType << "\t" << CurrentVariant->Chromosome << ":" << CurrentVariant->Position << "\t";
					
					for (std::vector<SnpEffDB>::iterator CurrentAnnotation = CurrentVariant->AnnotationDB.begin(); CurrentAnnotation != CurrentVariant->AnnotationDB.end(); ++CurrentAnnotation)
					{
						if (CurrentTranscript.Name != CurrentAnnotation->TranscriptID || CurrentAnnotation->EffectImpact == "MODIFIER") continue;
						std::cout << CurrentAnnotation->FunctionalClass << "," << CurrentAnnotation->EffectImpact << "," << CurrentAnnotation->Effect;
					}
					
					std::cout << "\t";
					
					for (std::vector<Individual>::iterator CurrentIndividual = CurrentPedigree.Families[CurrentFamily->first].begin(); CurrentIndividual != CurrentPedigree.Families[CurrentFamily->first].end(); ++CurrentIndividual)
					{
						if (CurrentIndividual->AffectionStatus == 2)
						{
							for (std::vector<std::string>::iterator Inheritance = CurrentVariant->InheritancePatterns[CurrentIndividual->IndividualID].begin(); Inheritance != CurrentVariant->InheritancePatterns[CurrentIndividual->IndividualID].end(); ++Inheritance)
							{
								std::cout << *Inheritance << "->" <<  CurrentIndividual->IndividualID;
								std::cout << ",";
							}
						}
					}
					
					
					if (CurrentVariant->Novel && !CurrentVariant->IsSecondary)
					{
						std::cout << "\tdbSNP-novel";
					}
					else
					{
						for (std::vector<std::string>::iterator CurrentID = CurrentVariant->ID.begin(); CurrentID != CurrentVariant->ID.end(); ++CurrentID)
						{
							std::cout << *CurrentID;
						}
					}
					
					for (auto CurrentAnnotation = AnnotationsToInclude.begin(); CurrentAnnotation != AnnotationsToInclude.end(); ++CurrentAnnotation)
					{
						std::cout << "\t" << CurrentVariant->Annotations[*CurrentAnnotation];
					}
					
					//				if (CurrentVariant->IsSecondary) std::cout << "(secondary)";
					
					std::cout << "\n";
				}
			}
		}
	}
	return;
}

#endif