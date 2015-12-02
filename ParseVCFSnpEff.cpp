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



struct VCFHeaderPrototype
{
	std::vector<std::string> HeaderLines;
	std::vector<std::string> SampleIDs;
	std::map<std::string,std::vector<std::string> > Metadata;
};

// EFF= Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_Change| Amino_Acid_Length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank  | Genotype_Number [ | ERRORS | WARNINGS ] )
// stop_lost(HIGH|MISSENSE|tGa/tTa|p.Ter444Leuext*?/c.1331G>T|477|LOC101929983|protein_coding|CODING|NM_001300891.2|3|1|WARNING_TRANSCRIPT_INCOMPLETE)
// missense_variant(MODERATE|MISSENSE|Ctg/Atg|p.Leu311Met/c.931C>A|444|ANGPT2|protein_coding|CODING|NM_001118888|6|1)


struct SnpEffDB
{
	std::string Effect;
	std::string EffectImpact;
	std::string FunctionalClass;
	std::string CodonChange;
	std::string AminoAcidChange;
	std::string AminoAcidLength;
	std::string GeneName;
	std::string TranscriptBiotype;
	std::string GeneCoding;
	std::string TranscriptID;
	std::string ExonRank;
	std::string GenotypeNumber;
};

struct Variant
{
	std::string Chromosome;
	unsigned long int Position;
	unsigned long int EndPosition;
	std::vector<std::string> ID;
	std::string Reference;
	std::vector<std::string> Alleles;
	double Quality;
	std::vector<std::string> Filter;
	std::map<std::string,std::vector<std::string> > Info;
	std::vector<std::string> Format;
	
	std::map<std::string,std::pair<std::string,std::string> > Genotypes;
	std::map<std::string,std::vector<std::string> > InheritancePatterns;
	std::vector<SnpEffDB> AnnotationDB;
	std::map<std::string,std::string> Annotations;

	bool IsSecondary;
	bool Novel;
	bool PassedFilters;
	bool PassedFrequencyCutoffs;
};

void ParseSnpEff(std::string& SnpEffLine, Variant& CurrentVariant)
{
	std::vector<std::string> SplitSnpEffLine;
	
	tokenize(SnpEffLine,SplitSnpEffLine,",");
	
	for (std::vector<std::string>::iterator CurrentAnnotation = SplitSnpEffLine.begin(); CurrentAnnotation != SplitSnpEffLine.end(); ++CurrentAnnotation)
	{
		SnpEffDB SnpEffAnnotation;
		
		std::vector<std::string> ParsedAnnotations;
		SnpEffAnnotation.Effect = CurrentAnnotation->substr(0,CurrentAnnotation->find_first_of("("));
		std::string TrimmedAnnotation = CurrentAnnotation->substr(CurrentAnnotation->find_first_of("(")+1);
		TrimmedAnnotation.pop_back();
		tokenize(TrimmedAnnotation,ParsedAnnotations,"|");
		
		SnpEffAnnotation.EffectImpact = ParsedAnnotations[0];
		SnpEffAnnotation.FunctionalClass = ParsedAnnotations[1];
		SnpEffAnnotation.CodonChange = ParsedAnnotations[2];
		SnpEffAnnotation.AminoAcidChange = ParsedAnnotations[3];
		SnpEffAnnotation.AminoAcidLength = ParsedAnnotations[4];
		SnpEffAnnotation.GeneName = ParsedAnnotations[5];
		SnpEffAnnotation.TranscriptBiotype = ParsedAnnotations[6];
		SnpEffAnnotation.GeneCoding = ParsedAnnotations[7];
		SnpEffAnnotation.TranscriptID = ParsedAnnotations[8];
		SnpEffAnnotation.ExonRank = ParsedAnnotations[9];
		SnpEffAnnotation.GenotypeNumber = ParsedAnnotations[10];
		
		CurrentVariant.AnnotationDB.push_back(SnpEffAnnotation);
	}
}

Variant ParseVCF(std::string& VCFLine,VCFHeaderPrototype& Header)
{ 
	Variant CurrentVariant;
	
	std::vector<std::string> SplitVCFLine;
	tokenize(VCFLine, SplitVCFLine);
	
	CurrentVariant.Chromosome = ProcessChromosomeString(SplitVCFLine[0]);
	
	CurrentVariant.Position = str2int(SplitVCFLine[1]);
	
	tokenize(SplitVCFLine[2],CurrentVariant.ID,";");
	if (std::find(CurrentVariant.ID.begin(), CurrentVariant.ID.end(), ".") != CurrentVariant.ID.end()) CurrentVariant.Novel = true;
	else CurrentVariant.Novel = false;
	
	CurrentVariant.Reference = SplitVCFLine[3];
	CurrentVariant.Alleles.push_back(CurrentVariant.Reference);
	
	tokenize(SplitVCFLine[4],CurrentVariant.Alleles,",");
	
	CurrentVariant.Quality = atof(SplitVCFLine[5].c_str());
	
	tokenize(SplitVCFLine[6],CurrentVariant.Filter,";");
	
	std::vector<std::string> SplitInfoLine;
	tokenize(SplitVCFLine[7],SplitInfoLine,";");
	
	std::vector<std::string> SplitAnnotations;
	for (std::vector<std::string>::iterator InfoLineIterator = SplitInfoLine.begin(); InfoLineIterator != SplitInfoLine.end(); ++InfoLineIterator)
	{
		std::vector<std::string> SplitInfoBit;
		tokenize(*InfoLineIterator,SplitInfoBit,"=");
		if (SplitInfoBit[0] == "EFF" || SplitInfoBit[0] == "SPL")
		{
			ParseSnpEff(SplitInfoBit[1],CurrentVariant);
			continue;
		}
		if (SplitInfoBit.size() > 1) CurrentVariant.Annotations[SplitInfoBit[0]] = SplitInfoBit[1];
	}
	
	//Rectify annotation bin
	
	tokenize(SplitVCFLine[8],CurrentVariant.Format,":");
	
	std::vector<std::string>::iterator SampleNameIterator = Header.SampleIDs.begin();
	for (std::vector<std::string>::iterator SampleIterator = SplitVCFLine.begin()+9; SampleIterator != SplitVCFLine.end(); ++SampleIterator)
	{
		if (*SampleIterator == "./." || *SampleIterator == ".")
		{
			CurrentVariant.Genotypes[*SampleNameIterator]=std::pair<std::string,std::string>(".",".");
			++SampleNameIterator;
			continue;
		}

		std::vector<std::string> SplitFormatLine,SplitGenotypeLine;
		tokenize(*SampleIterator,SplitFormatLine,":");
		tokenize(SplitFormatLine[0],SplitGenotypeLine,"/");
		
		if (SplitGenotypeLine[0] == "." )
		{
			CurrentVariant.Genotypes[*SampleNameIterator]=std::pair<std::string,std::string>(".",".");
			++SampleNameIterator;
			continue;
		}
		if (SplitGenotypeLine[0] != "." || SplitGenotypeLine[1] != ".")
		{
			CurrentVariant.Genotypes[*SampleNameIterator]=std::pair<std::string,std::string>(CurrentVariant.Alleles[str2int(SplitGenotypeLine[0])],CurrentVariant.Alleles[str2int(SplitGenotypeLine[1])]);
		}
		else
		{
			CurrentVariant.Genotypes[*SampleNameIterator]=std::pair<std::string,std::string>(".",".");
		}
		++SampleNameIterator;
	}
	
	CurrentVariant.PassedFrequencyCutoffs = true;
	return(CurrentVariant);
}

VCFHeaderPrototype ParseVCFHeader(std::ifstream& InputVCF,Pedigree& FamilyStructure)
{
	VCFHeaderPrototype VCFHeader;
	while (true)
	{
		std::string HeaderLine;
		getline(InputVCF, HeaderLine);
		if (InputVCF.peek()!='#')
		{
			std::vector<std::string> SplitHeaderLine;
			tokenize(HeaderLine,SplitHeaderLine);
			for (std::vector<std::string>::iterator HeaderIterator = SplitHeaderLine.begin()+9;HeaderIterator != SplitHeaderLine.end(); ++HeaderIterator)
			{
				VCFHeader.SampleIDs.push_back(*HeaderIterator);
			}
			
			for (std::map<std::string,std::vector<Individual> >::iterator CurrentFamily = FamilyStructure.Families.begin();CurrentFamily != FamilyStructure.Families.end(); ++CurrentFamily)
			{
				for (std::vector<Individual>::iterator CurrentIndividual = CurrentFamily->second.begin(); CurrentIndividual != CurrentFamily->second.end(); ++CurrentIndividual)
				{
					if (std::find(VCFHeader.SampleIDs.begin(), VCFHeader.SampleIDs.end(), CurrentIndividual->IndividualID) == VCFHeader.SampleIDs.end())
					{
						std::clog << "WARNING: Unable to find " << CurrentIndividual->IndividualID << " in VCF file and will use null genotypes in all tests.\n";
					}
				}
			}
			break;
		}
		std::vector<std::string> SplitHeaderLine;
		tokenize(HeaderLine,SplitHeaderLine,"=");
		VCFHeader.Metadata[SplitHeaderLine[0]].push_back(SplitHeaderLine[1]);
	}
	return(VCFHeader);
}
