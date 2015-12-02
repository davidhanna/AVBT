#ifndef GENETEST_H
#define GENETEST_H

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

class TwoHitModel
{
	public:
	TwoHitModel();
	void OutputTestResults(TranscriptRNA&, Pedigree&, std::map<std::string,std::vector<Individual> >::iterator, std::vector<std::string>& );
	void SexLinkedTest(TranscriptRNA&, Pedigree&, std::map<std::string,std::vector<Individual> >::iterator);
	void AutosomalTest(TranscriptRNA&, Pedigree&, std::map<std::string,std::vector<Individual> >::iterator);
	bool VariationFitsModel;
	
	std::string ModelType;
	std::vector<Variant> PrimaryVariants;
	std::vector<Variant> SecondaryVariants;
};

TwoHitModel::TwoHitModel()
{
	VariationFitsModel = true;
}

// MODIFY FOR SECONDARY VARIANTS
void TwoHitModel::OutputTestResults(TranscriptRNA& CurrentTranscript,Pedigree& CurrentPedigree, std::map<std::string,std::vector<Individual> >::iterator CurrentFamily, std::vector<std::string>& AnnotationsToInclude)
{

	for (std::vector<Variant>::iterator CurrentVariant = PrimaryVariants.begin(); CurrentVariant != PrimaryVariants.end(); ++CurrentVariant)
	{
		std::cout << CurrentTranscript.GeneName << "\t" << CurrentTranscript.Name << "\t" << CurrentFamily->first << "\t" << ModelType << "\t";
		std::cout << CurrentVariant->Chromosome << ":" << CurrentVariant->Position << "\t";
		
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
						
						std::vector<std::string>::iterator NextInheritance = Inheritance;
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
				std::cout << "\t" << *CurrentID;
			}
			
		}
		for (auto CurrentAnnotation = AnnotationsToInclude.begin(); CurrentAnnotation != AnnotationsToInclude.end(); ++CurrentAnnotation)
		{
			std::cout << "\t" << CurrentVariant->Annotations[*CurrentAnnotation];
		}
		std::cout << "\n";
	}
//	for (std::vector<Variant>::iterator CurrentVariant = SecondaryVariants.begin(); CurrentVariant != SecondaryVariants.end(); ++CurrentVariant)
//	{
//		std::cout << CurrentVariant->Chromosome << ":" << CurrentVariant->Position << "|UNKNOWN|(secondary)\t";
//	}
//	std::cout << "\n";
}

void TwoHitModel::SexLinkedTest(TranscriptRNA& CurrentTranscript, Pedigree& CurrentPedigree, std::map<std::string,std::vector<Individual> >::iterator CurrentFamily)
{
	for (std::vector<Individual>::iterator IndividualIterator = CurrentFamily->second.begin(); IndividualIterator != CurrentFamily->second.end();++IndividualIterator)
	{
		if (IndividualIterator->AffectionStatus == 2)
		{
			bool HasMotherAllele = false; bool HasFatherAllele = false; bool HasEitherAllele = false; unsigned short AlleleCount = 0;
			bool IsFounderOrMissing = false;
			
			for (std::vector<Variant>::iterator CurrentVariant = PrimaryVariants.begin(); CurrentVariant != PrimaryVariants.end(); ++CurrentVariant)
			{
				for (std::vector<std::string>::iterator RecordIterator = CurrentVariant->InheritancePatterns[IndividualIterator->IndividualID].begin(); RecordIterator != CurrentVariant->InheritancePatterns[IndividualIterator->IndividualID].end(); ++RecordIterator)
				{
					++AlleleCount;
					if (*RecordIterator == ".") IsFounderOrMissing = true;
					if (*RecordIterator == IndividualIterator->MotherID) HasMotherAllele = true;
					if (*RecordIterator == IndividualIterator->FatherID) HasFatherAllele = true;
					if (*RecordIterator == "EITHER") HasEitherAllele = true;
				}
			}
			if (IsFounderOrMissing)
			{
				if ( (AlleleCount < 2 && IndividualIterator->Sex == 2) || (AlleleCount < 1 && IndividualIterator->Sex == 1) )
				{
					VariationFitsModel = false;
					break;
				}
			}
			else if ( ( ( (!HasMotherAllele && !HasEitherAllele) || (!HasFatherAllele && !HasEitherAllele) ) && IndividualIterator->Sex == 2 ) || (AlleleCount == 0 && IndividualIterator->Sex == 1))
			{
				VariationFitsModel = false;
				break;
			}
		}
		
		ModelType = "X-Linked-Recessive";
		
//		for (std::vector<Individual>::iterator IndividualPrimeIterator = CurrentFamily->second.begin(); IndividualPrimeIterator != CurrentFamily->second.end();++IndividualPrimeIterator)
//		{
//			if (IndividualIterator==IndividualPrimeIterator) continue;
//
//			bool SameParents = false; unsigned short Matches = 0; unsigned short PhenotypeComparison = IndividualIterator->AffectionStatus + IndividualPrimeIterator->AffectionStatus;
//			
//			if (PhenotypeComparison == 2) continue;
//			
//			if (IndividualIterator->MotherID == IndividualPrimeIterator->MotherID && IndividualIterator->FatherID == IndividualPrimeIterator->FatherID && IndividualIterator->FatherID != ".") SameParents = true;
//
//			bool CompleteMatch = true;
//			for (std::vector<Variant>::iterator CurrentVariant = PrimaryVariants.begin(); CurrentVariant != PrimaryVariants.end(); ++CurrentVariant)
//			{
//				if (CurrentVariant->IsSecondary) continue;
//				if (CurrentVariant->Genotypes[IndividualIterator->IndividualID] != CurrentVariant->Genotypes[IndividualPrimeIterator->IndividualID]) CompleteMatch = false;
//				
//				if (CurrentVariant->Genotypes[IndividualIterator->IndividualID].first != CurrentVariant->Reference && CurrentVariant->Genotypes[IndividualIterator->IndividualID].first != "." && CurrentVariant->Genotypes[IndividualIterator->IndividualID].first == CurrentVariant->Genotypes[IndividualPrimeIterator->IndividualID].first) Matches++;
//				if (CurrentVariant->Genotypes[IndividualIterator->IndividualID].second != CurrentVariant->Reference && CurrentVariant->Genotypes[IndividualIterator->IndividualID].second != "." && CurrentVariant->Genotypes[IndividualIterator->IndividualID].second == CurrentVariant->Genotypes[IndividualPrimeIterator->IndividualID].second) Matches++;
//			}
//			if ( (SameParents && PhenotypeComparison == 4 && !CompleteMatch) || PhenotypeComparison == 3 && CompleteMatch || (!SameParents && PhenotypeComparison == 4 && Matches == 0) )
//			{
//				VariationFitsModel = false;
//				break;
//			}
//		}
	}
}

void TwoHitModel::AutosomalTest(TranscriptRNA& CurrentTranscript, Pedigree& CurrentPedigree, std::map<std::string,std::vector<Individual> >::iterator CurrentFamily)
{	
	for (std::vector<Individual>::iterator IndividualIterator = CurrentFamily->second.begin(); IndividualIterator != CurrentFamily->second.end();++IndividualIterator)
	{
		if (IndividualIterator->AffectionStatus == 2)
		{
			bool HasMotherAllele = false; bool HasFatherAllele = false; bool HasEitherAllele = false; unsigned short AlleleCount = 0; unsigned short BigAlleleCount = 0;
			bool IsFounderOrMissing = false;
			
			for (std::vector<Variant>::iterator CurrentVariant = PrimaryVariants.begin(); CurrentVariant != PrimaryVariants.end(); ++CurrentVariant)
			{				 
				++BigAlleleCount;
				for (std::vector<std::string>::iterator RecordIterator = CurrentVariant->InheritancePatterns[IndividualIterator->IndividualID].begin(); RecordIterator != CurrentVariant->InheritancePatterns[IndividualIterator->IndividualID].end(); ++RecordIterator)
				{
					++AlleleCount;
					if (*RecordIterator == ".") IsFounderOrMissing = true;
					if (*RecordIterator == IndividualIterator->MotherID) HasMotherAllele = true;
					if (*RecordIterator == IndividualIterator->FatherID) HasFatherAllele = true;
					if (*RecordIterator == "EITHER") HasEitherAllele = true;
				}
			}
			if (IsFounderOrMissing)
			{
				if (AlleleCount < 2)
				{
					VariationFitsModel = false;
					break;
				}
			}
			else if ( (!HasMotherAllele && !HasEitherAllele) || (!HasFatherAllele && !HasEitherAllele) )
			{
				VariationFitsModel = false;
				break;
			}
			
			if (BigAlleleCount == 1) ModelType = "Autosomal-Homozygous-Recessive";
			else ModelType = "Autosomal-Compound-Heterozygous-Recessive";
		}

		for (std::vector<Individual>::iterator IndividualPrimeIterator = CurrentFamily->second.begin(); IndividualPrimeIterator != CurrentFamily->second.end();++IndividualPrimeIterator)
		{
			if (IndividualIterator==IndividualPrimeIterator) continue;

			bool SameParents = false; unsigned short Matches = 0; unsigned short PhenotypeComparison = IndividualIterator->AffectionStatus + IndividualPrimeIterator->AffectionStatus;
			
			if (PhenotypeComparison == 2) continue;
			
			if (IndividualIterator->MotherID == IndividualPrimeIterator->MotherID && IndividualIterator->FatherID == IndividualPrimeIterator->FatherID && IndividualIterator->FatherID != ".") SameParents = true;

			bool CompleteMatch = true;
			for (std::vector<Variant>::iterator CurrentVariant = PrimaryVariants.begin(); CurrentVariant != PrimaryVariants.end(); ++CurrentVariant)
			{
				if (CurrentVariant->IsSecondary) continue;
				if (CurrentVariant->Genotypes[IndividualIterator->IndividualID] != CurrentVariant->Genotypes[IndividualPrimeIterator->IndividualID]) CompleteMatch = false;
				
				if (CurrentVariant->Genotypes[IndividualIterator->IndividualID].first != CurrentVariant->Reference && CurrentVariant->Genotypes[IndividualIterator->IndividualID].first != "." && CurrentVariant->Genotypes[IndividualIterator->IndividualID].first == CurrentVariant->Genotypes[IndividualPrimeIterator->IndividualID].first) Matches++;
				if (CurrentVariant->Genotypes[IndividualIterator->IndividualID].second != CurrentVariant->Reference && CurrentVariant->Genotypes[IndividualIterator->IndividualID].second != "." && CurrentVariant->Genotypes[IndividualIterator->IndividualID].second == CurrentVariant->Genotypes[IndividualPrimeIterator->IndividualID].second) Matches++;
			}
			if ( (SameParents && PhenotypeComparison == 4 && !CompleteMatch) || (PhenotypeComparison == 3 && CompleteMatch)  || (!SameParents && PhenotypeComparison == 4 && Matches == 0) )
			{
				VariationFitsModel = false;
				break;
			}
		}
	}
	
}

void TestGeneWithRecessiveModel(TranscriptRNA& CurrentTranscript, Pedigree& CurrentPedigree, std::vector<std::string>& AnnotationsToInclude)
{
	if (CurrentTranscript.TranscriptVariants.size() == 0 || CurrentTranscript.Chromosome == "Y") return;
		
	for (std::map<std::string,std::vector<Individual> >::iterator CurrentFamily = CurrentPedigree.Families.begin(); CurrentFamily != CurrentPedigree.Families.end();++CurrentFamily)
	{
		TwoHitModel RecessiveGeneStatistic;

		// Find suitable variants and pass them into sample structures for pairwise comparison
		
		for (std::vector<Variant>::iterator CurrentVariant = CurrentTranscript.TranscriptVariants.begin(); CurrentVariant != CurrentTranscript.TranscriptVariants.end(); ++CurrentVariant)
		{
			if (CurrentVariant->IsSecondary)
			{
				RecessiveGeneStatistic.SecondaryVariants.push_back(*CurrentVariant);
				continue;
			}
			
			bool SufficientSeverity = false;
			for (std::vector<SnpEffDB>::iterator CurrentAnnotationSeverity = CurrentVariant->AnnotationDB.begin(); CurrentAnnotationSeverity != CurrentVariant->AnnotationDB.end(); ++CurrentAnnotationSeverity)
			{
				if (CurrentAnnotationSeverity->EffectImpact == "HIGH" || CurrentAnnotationSeverity->EffectImpact == "MODERATE" || CurrentAnnotationSeverity->EffectImpact == "LOW"  && CurrentTranscript.Name == CurrentAnnotationSeverity->TranscriptID)
				{
					SufficientSeverity = true;
					break;
				}
			}
			if (CurrentVariant->Annotations["dbNSFP_clinvar_clnsig"] != "." && !CurrentVariant->Annotations["dbNSFP_clinvar_clnsig"].empty()) SufficientSeverity = true;
			
			if (!CurrentVariant->PassedFilters || !CurrentVariant->PassedFrequencyCutoffs || !SufficientSeverity) continue;
			bool HomozygousInUnaffected = false; bool InAffected = false;
			for (std::vector<Individual>::iterator IndividualIterator = CurrentFamily->second.begin(); IndividualIterator != CurrentFamily->second.end();++IndividualIterator)
			{
				if (IndividualIterator->AffectionStatus == 1)
				{
					if (CurrentVariant->Genotypes[IndividualIterator->IndividualID].first != CurrentVariant->Reference && CurrentVariant->Genotypes[IndividualIterator->IndividualID].first != "." && CurrentVariant->Genotypes[IndividualIterator->IndividualID].second != CurrentVariant->Reference && CurrentVariant->Genotypes[IndividualIterator->IndividualID].second != ".") HomozygousInUnaffected = true;
				}
				if (IndividualIterator->AffectionStatus == 2)
				{
					if ( (CurrentVariant->Genotypes[IndividualIterator->IndividualID].first != CurrentVariant->Reference && CurrentVariant->Genotypes[IndividualIterator->IndividualID].first != "." ) || (CurrentVariant->Genotypes[IndividualIterator->IndividualID].second != CurrentVariant->Reference && CurrentVariant->Genotypes[IndividualIterator->IndividualID].second != ".")) InAffected = true;
				}
			}
			if (InAffected && !HomozygousInUnaffected) RecessiveGeneStatistic.PrimaryVariants.push_back(*CurrentVariant);
		}
		
		// All variants should be collected for this family, now collected in gene statistic module
		if (RecessiveGeneStatistic.PrimaryVariants.empty() && RecessiveGeneStatistic.SecondaryVariants.empty())
		{
			RecessiveGeneStatistic.VariationFitsModel = false;
			continue;
		}
		
		// Begin tests for model
		if (CurrentTranscript.Chromosome == "X") RecessiveGeneStatistic.SexLinkedTest(CurrentTranscript, CurrentPedigree, CurrentFamily);
		if (CurrentTranscript.Chromosome != "X" && CurrentTranscript.Chromosome != "Y") RecessiveGeneStatistic.AutosomalTest(CurrentTranscript, CurrentPedigree, CurrentFamily);
		
		// Output results
		if (RecessiveGeneStatistic.VariationFitsModel || !RecessiveGeneStatistic.SecondaryVariants.empty()) RecessiveGeneStatistic.OutputTestResults(CurrentTranscript,CurrentPedigree, CurrentFamily, AnnotationsToInclude);
	}
	return;
}


void TestGeneWithSingleHitModel(TranscriptRNA& CurrentTranscript, Pedigree& CurrentPedigree, std::vector<std::string>& AnnotationsToInclude)
{
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
			
			if (CurrentVariant->Chromosome == "X" || CurrentVariant->Chromosome == "Y") ModelType = CurrentVariant->Chromosome + "-Linked-Dominant";
			else ModelType = "Autosomal-Dominant";

			short TotalAlleles = 0;
			
			for (std::vector<Individual>::iterator IndividualIterator = CurrentFamily->second.begin(); IndividualIterator != CurrentFamily->second.end();++IndividualIterator)
			{
				if (CurrentVariant->Genotypes.find(IndividualIterator->IndividualID) == CurrentVariant->Genotypes.end()) CurrentVariant->Genotypes[IndividualIterator->IndividualID] = std::pair<std::string,std::string>(".",".");
								
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

					if ((AlleleCount == 1 && ModelType == "Autosomal-Dominant") || ( AlleleCount == 2 && ModelType != "Autosomal-Dominant" && IndividualIterator->Sex == 1))
					{
						for (std::vector<std::string>::iterator CurrentInheritance = CurrentVariant->InheritancePatterns[IndividualIterator->IndividualID].begin(); CurrentInheritance != CurrentVariant->InheritancePatterns[IndividualIterator->IndividualID].end(); ++CurrentInheritance) if (*CurrentInheritance == "DeNovo") ModelType = "DeNovo";
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
	return;
}
#endif