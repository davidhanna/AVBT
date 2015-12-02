
// This verison is seeking to add tests for Mitochrondrial sequences, mosaicism, imprinted genes, repeat triplet disorders and other non-mendelian inheritance

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
#include "ParseVCFSnpEff.cpp"
#include "ParseRefSeq.cpp"
#include "CommonAlgorithms.cpp"
#include "GeneTestsIII.cpp"
#include "Pedigree.cpp"
#include "Frequencies.cpp"
#include "ParseBED.cpp"
#include "IdentifyTransmission.cpp"
#include "Settings.cpp"
#include "Imprinted.cpp"

int main(int argc, char** argv)
{
	ProgramOptions Settings;
	for (int i=1;i<argc;i++)
	{
		if ( (std::string)argv[i] == "-RefSeq" ) Settings.RefSeqPath = (std::string) argv[i+1];
		if ( (std::string)argv[i] == "-CutoffFrequency" ) Settings.FrequenciesCutoff = ::atof(argv[i+1]);
		if ( (std::string)argv[i] == "-VCF" ) Settings.VCFFilePaths.push_back(std::make_pair((std::string) argv[i+1],false));
		if ( (std::string)argv[i] == "-SupplementalVCF" ) Settings.VCFFilePaths.push_back(std::make_pair((std::string) argv[i+1],true));
		if ( (std::string)argv[i] == "-PED" ) Settings.PedigreePath = (std::string) argv[i+1];
		if ( (std::string)argv[i] == "-BED" ) Settings.BEDPath = (std::string) argv[i+1];
		if ( (std::string)argv[i] == "-QualityCutoff" ) Settings.QualityCutoff = ::atof(argv[i+1]);
		if ( (std::string)argv[i] == "-ArtifactCutoff" ) Settings.ArtifactCutoff = ::atoi(argv[i+1]);
		
	}
	
	std::vector<std::string> AnnotationsToInclude;
	{
	AnnotationsToInclude.push_back("ExAC_AF");
	AnnotationsToInclude.push_back("dbNSFP_Interpro_domain");
	AnnotationsToInclude.push_back("dbNSFP_SLR_test_statistic");
	AnnotationsToInclude.push_back("dbNSFP_Ancestral_allele");
	AnnotationsToInclude.push_back("dbNSFP_Ensembl_transcriptid");
	AnnotationsToInclude.push_back("dbNSFP_SIFT_score");
	AnnotationsToInclude.push_back("dbNSFP_SIFT_converted_rankscore");
	AnnotationsToInclude.push_back("dbNSFP_SIFT_pred");
	AnnotationsToInclude.push_back("dbNSFP_Polyphen2_HDIV_score");
	AnnotationsToInclude.push_back("dbNSFP_Polyphen2_HDIV_rankscore");
	AnnotationsToInclude.push_back("dbNSFP_Polyphen2_HDIV_pred");
	AnnotationsToInclude.push_back("dbNSFP_Polyphen2_HVAR_score");
	AnnotationsToInclude.push_back("dbNSFP_Polyphen2_HVAR_rankscore");
	AnnotationsToInclude.push_back("dbNSFP_Polyphen2_HVAR_pred");
	AnnotationsToInclude.push_back("dbNSFP_LRT_score");
	AnnotationsToInclude.push_back("dbNSFP_LRT_converted_rankscore");
	AnnotationsToInclude.push_back("dbNSFP_LRT_pred");
	AnnotationsToInclude.push_back("dbNSFP_MutationTaster_score");
	AnnotationsToInclude.push_back("dbNSFP_MutationTaster_converted_rankscore");
	AnnotationsToInclude.push_back("dbNSFP_MutationTaster_pred");
	AnnotationsToInclude.push_back("dbNSFP_MutationAssessor_score");
	AnnotationsToInclude.push_back("dbNSFP_MutationAssessor_rankscore");
	AnnotationsToInclude.push_back("dbNSFP_MutationAssessor_pred");
	AnnotationsToInclude.push_back("dbNSFP_FATHMM_score");
	AnnotationsToInclude.push_back("dbNSFP_FATHMM_rankscore");
	AnnotationsToInclude.push_back("dbNSFP_FATHMM_pred");
	AnnotationsToInclude.push_back("dbNSFP_RadialSVM_score");
	AnnotationsToInclude.push_back("dbNSFP_RadialSVM_rankscore");
	AnnotationsToInclude.push_back("dbNSFP_RadialSVM_pred");
	AnnotationsToInclude.push_back("dbNSFP_LR_score");
	AnnotationsToInclude.push_back("dbNSFP_LR_rankscore");
	AnnotationsToInclude.push_back("dbNSFP_LR_pred");
	AnnotationsToInclude.push_back("dbNSFP_Reliability_index");
	AnnotationsToInclude.push_back("dbNSFP_VEST3_score");
	AnnotationsToInclude.push_back("dbNSFP_VEST3_rankscore");
	AnnotationsToInclude.push_back("dbNSFP_CADD_raw");
	AnnotationsToInclude.push_back("dbNSFP_CADD_raw_rankscore");
	AnnotationsToInclude.push_back("dbNSFP_CADD_phred");
	AnnotationsToInclude.push_back("dbNSFP_GERP++_NR");
	AnnotationsToInclude.push_back("dbNSFP_GERP++_RS");
	AnnotationsToInclude.push_back("dbNSFP_GERP++_RS_rankscore");
	AnnotationsToInclude.push_back("dbNSFP_phyloP46way_primate");
	AnnotationsToInclude.push_back("dbNSFP_phyloP46way_primate_rankscore");
	AnnotationsToInclude.push_back("dbNSFP_phyloP46way_placental");
	AnnotationsToInclude.push_back("dbNSFP_phyloP46way_placental_rankscore");
	AnnotationsToInclude.push_back("dbNSFP_phyloP100way_vertebrate");
	AnnotationsToInclude.push_back("dbNSFP_phyloP100way_vertebrate_rankscore");
	AnnotationsToInclude.push_back("dbNSFP_phastCons46way_primate");
	AnnotationsToInclude.push_back("dbNSFP_phastCons46way_primate_rankscore");
	AnnotationsToInclude.push_back("dbNSFP_phastCons46way_placental");
	AnnotationsToInclude.push_back("dbNSFP_phastCons46way_placental_rankscore");
	AnnotationsToInclude.push_back("dbNSFP_phastCons100way_vertebrate");
	AnnotationsToInclude.push_back("dbNSFP_phastCons100way_vertebrate_rankscore");
	AnnotationsToInclude.push_back("dbNSFP_SiPhy_29way_pi");
	AnnotationsToInclude.push_back("dbNSFP_SiPhy_29way_logOdds");
	AnnotationsToInclude.push_back("dbNSFP_SiPhy_29way_logOdds_rankscore");
	AnnotationsToInclude.push_back("dbNSFP_LRT_Omega");
	AnnotationsToInclude.push_back("dbNSFP_clinvar_rs");
	AnnotationsToInclude.push_back("dbNSFP_clinvar_clnsig");
	AnnotationsToInclude.push_back("dbNSFP_clinvar_trait");
	AnnotationsToInclude.push_back("ArtifactAC");
}
	
	
	std::vector<std::string> ImprintedGenes;
	std::ifstream ImprintedGeneFile;
	ImprintedGeneFile.open("/data/vol2/home/dshanna/Imprinting/GeneList.simple.txt");
	std::string CurrentImprintedGene;
	while(getline(ImprintedGeneFile,CurrentImprintedGene))
	{
		ImprintedGenes.push_back(CurrentImprintedGene);
	}
	
	std::clog << "\n\t\t      __      ______ _______ \n\t\t     /\\ \\    / /  _ \\__   __|\n\t\t    /  \\ \\  / /| |_) | | |   \n\t\t   / /\\ \\ \\/ / |  _ <  | |   \n\t\t  / ____ \\  /  | |_) | | |   \n\t\t /_/    \\_\\/   |____/  |_|   \n\t\t                             \n\t\t\t\t     D.S. Hanna, University of Washington\n\t\t\t\t     V2.0\n\n\n";
	
	std::string Usage = "/\\  Standard usage for this program:\n\t-RefSeq <path> (File containing gene definitions)\n\t-FrequencyFile <path> (File containing population allele frequencies)\n\t-CutoffFrequency <float> (Threshold for inclusion)\n\t-VCF <path> (Path to snpEff annotated VCF file, can be specified multiple times)\n\t-SupplementalVCF <path> (Path to VCF file with untrusted genotypes, can be specified multiple times)\n\t-PED <path> (Standard 6 column pedigree file)\n\t-BED <path> (BED format file for CNV calls)\n";
	if (argc<6){ std::cout << Usage; return 0;}
	
	std::map<std::string, std::vector<TranscriptRNA> > RefSeqTranscripts;
	LoadRefSeq(RefSeqTranscripts,Settings.RefSeqPath);
	
	Pedigree FamilyStructure;
	FamilyStructure.PopulatePedigree(Settings.PedigreePath);
	
	// Process all VCF files specified
	for (std::vector<std::pair<std::string,bool> >::iterator CurrentVCFPath = Settings.VCFFilePaths.begin(); CurrentVCFPath != Settings.VCFFilePaths.end(); ++CurrentVCFPath)
	{
		std::ifstream VCFFile;
		VCFFile.open(CurrentVCFPath->first.c_str());
		
		if (!VCFFile.is_open()) { std::cerr << "Cannot open VCF, exiting\n"; exit(0);}
		VCFHeaderPrototype VCFHeader = ParseVCFHeader(VCFFile,FamilyStructure);
		
		std::clog << "/\\ Currently processing VCF file " << CurrentVCFPath->first << "\n";

		Variant CurrentVariant;
		std::string VCFLine;	

		unsigned int LineCount = 0;
		unsigned int Count = 0;
		
		while(getline(VCFFile,VCFLine))
		{
			LineCount++;
			if (LineCount == 10000)
			{
				Count++; LineCount = 0;
				std::clog << "/\\ " << Count*10000 << " variants processed\n";
			}
			
			CurrentVariant = ParseVCF(VCFLine,VCFHeader);
			
			FilterVariant(CurrentVariant,Settings);
			
			// Find if variant is secondary
			if (CurrentVCFPath->second) CurrentVariant.IsSecondary = true;
			else CurrentVariant.IsSecondary = false;
			
			// Place variants
			for (std::vector<TranscriptRNA>::iterator TranscriptIterator = RefSeqTranscripts[CurrentVariant.Chromosome].begin(); TranscriptIterator != RefSeqTranscripts[CurrentVariant.Chromosome].end(); ++TranscriptIterator)
			{
				if (TranscriptIterator->InRange(CurrentVariant))
				{
					if (CurrentVariant.IsSecondary)
					{
						bool FoundByPrimary = false;
						for (std::vector<Variant>::iterator PreviousVariants = TranscriptIterator->TranscriptVariants.begin();PreviousVariants != TranscriptIterator->TranscriptVariants.end();++PreviousVariants)
						{
							if (PreviousVariants->Position == CurrentVariant.Position && CurrentVariant.Alleles == PreviousVariants->Alleles) FoundByPrimary = true;
						}
						if (FoundByPrimary == false) TranscriptIterator->TranscriptVariants.push_back(CurrentVariant);
					}
					else
					{
						 TranscriptIterator->TranscriptVariants.push_back(CurrentVariant);
					}
				}
			}
		}
		VCFFile.close();
	}
	
	// Process BED files
	if (!Settings.BEDPath.empty())
	{
		std::ifstream InputBED;
		std::string BEDLine;
		InputBED.open(Settings.BEDPath.c_str());
		while(getline(InputBED,BEDLine))
		{
			Variant CurrentVariant = ParseBED(BEDLine);
			for (std::vector<TranscriptRNA>::iterator TranscriptIterator = RefSeqTranscripts[CurrentVariant.Chromosome].begin(); TranscriptIterator != RefSeqTranscripts[CurrentVariant.Chromosome].end(); ++TranscriptIterator)
			{
				if (TranscriptIterator->InRangeBED(CurrentVariant))
				{
					TranscriptIterator->TranscriptVariants.push_back(CurrentVariant);
				}
			}
		}
		InputBED.close();
	}
	
	// Begin analyses
	
	// Header
	
	std::cout << "Gene\tTranscript\tFamilyID\tModel\tPosition\tAnnotations\tVariantTransmissions\tdbSNP";
	for (auto CurrentAnnotation = AnnotationsToInclude.begin(); CurrentAnnotation != AnnotationsToInclude.end(); ++CurrentAnnotation)
	{
		std::cout << "\t" << *CurrentAnnotation;
	}
	std::cout << "\n";
	
	
	for (std::map<std::string,std::vector<TranscriptRNA> >::iterator ChromosomeTranscriptIterator = RefSeqTranscripts.begin(); ChromosomeTranscriptIterator != RefSeqTranscripts.end(); ++ChromosomeTranscriptIterator)
	{
		for (std::vector<TranscriptRNA>::iterator TranscriptIterator = ChromosomeTranscriptIterator->second.begin(); TranscriptIterator != ChromosomeTranscriptIterator->second.end(); ++TranscriptIterator)
		{
			IdentifyTransmission(*TranscriptIterator,FamilyStructure);
			TestGeneWithRecessiveModel(*TranscriptIterator,FamilyStructure,AnnotationsToInclude);
			TestGeneWithSingleHitModel(*TranscriptIterator,FamilyStructure,AnnotationsToInclude);
			TestGeneWithImprintedGeneModel(*TranscriptIterator,FamilyStructure,AnnotationsToInclude);
		}
	}
	
	std::clog << "/\\ Analysis complete.\n";
	
	return 0;
}
