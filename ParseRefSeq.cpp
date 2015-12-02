#ifndef PARSEREFSEQ_H
#define PARSEREFSEQ_H

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

class TranscriptRNA
{
	public:
	bool InRange(Variant&);
	bool InRangeBED(Variant&);
	bool InExon(Variant&);
	std::string Name;
	std::string Chromosome;
	bool Strand;
	unsigned long int TxStart;
	unsigned long int TxEnd;
	unsigned long int CDSStart;
	unsigned long int CDSEnd;
	unsigned long int ExonCounts;
	std::vector<unsigned long int> ExonStarts;
	std::vector<unsigned long int> ExonEnds;
	std::string GeneName;
	std::string CDSStartStat;
	std::string CDSEndStat;
	std::vector<int> ExonFrames;
	int Score;

	std::vector<Variant> TranscriptVariants;
};

bool TranscriptRNA::InRange(Variant& CurrentVariant)
{
	return ((TxStart > CurrentVariant.Position ? TxStart <= (CurrentVariant.Position + CurrentVariant.Reference.size()) : (CurrentVariant.Position + CurrentVariant.Reference.size()) <= TxEnd));
};

bool TranscriptRNA::InRangeBED(Variant& CurrentVariant)
{
	return (TxStart > CurrentVariant.Position ? TxStart <= CurrentVariant.EndPosition : (CurrentVariant.EndPosition <= TxEnd));
};


bool TranscriptRNA::InExon(Variant& CurrentVariant)
{
	std::vector<unsigned long int>::iterator Starts = ExonStarts.begin();
	std::vector<int>::iterator ExonFrame = ExonFrames.begin();
	for (std::vector<unsigned long int>::iterator Ends = ExonEnds.begin(); Ends != ExonEnds.end(); ++Ends)
	{
		if ( *Starts > CurrentVariant.Position ? *Starts <= (CurrentVariant.Position + CurrentVariant.Reference.size()) : (CurrentVariant.Position + CurrentVariant.Reference.size()) <= *Ends)
		{
			if ( CDSStart > CurrentVariant.Position ? CDSStart <= (CurrentVariant.Position + CurrentVariant.Reference.size()) : (CurrentVariant.Position + CurrentVariant.Reference.size()) <= CDSEnd)

			return true;
		}
	
		Starts++;
		ExonFrame++;
	}
	return false;
};

void LoadRefSeq(std::map<std::string, std::vector<TranscriptRNA> >& RefSeqTranscripts,std::string& RefSeqPath)
{
	std::ifstream RefSeqFile;
	std::string RefSeqLine;
	std::string Header;
	
	if (RefSeqPath.empty()) { std::cerr << "No RefSeq file specified, exiting\n"; exit(0); }
	
	RefSeqFile.open(RefSeqPath.c_str());
	
	if (!RefSeqFile.is_open()) { std::cerr << "RefSeq file could not be opened, exiting\n"; exit(0); }
	
	getline(RefSeqFile,Header);
	
	while(getline(RefSeqFile,RefSeqLine))
	{
		std::vector<std::string> SplitRefSeqLine;
		tokenize(RefSeqLine, SplitRefSeqLine);
		
		TranscriptRNA CurrentRNA;
		
		CurrentRNA.Name = SplitRefSeqLine[1];
		CurrentRNA.Chromosome = ProcessChromosomeString(SplitRefSeqLine[2]);
		if (SplitRefSeqLine[3][0] == '+') CurrentRNA.Strand = true;
		else CurrentRNA.Strand = false;
		
		CurrentRNA.TxStart = str2int(SplitRefSeqLine[4]);
		CurrentRNA.TxEnd = str2int(SplitRefSeqLine[5]);
		CurrentRNA.CDSStart = str2int(SplitRefSeqLine[6]);
		CurrentRNA.CDSEnd = str2int(SplitRefSeqLine[7]);
		CurrentRNA.ExonCounts = str2int(SplitRefSeqLine[8]);
		
		std::vector<std::string> SplitExonStarts;
		tokenize(SplitRefSeqLine[9], SplitExonStarts,",");
		for (std::vector<std::string>::iterator it = SplitExonStarts.begin(); it != SplitExonStarts.end();++it)
		{
			CurrentRNA.ExonStarts.push_back(str2int(*it));
		}

		std::vector<std::string> SplitExonEnds;
		tokenize(SplitRefSeqLine[10], SplitExonEnds,",");
		for (std::vector<std::string>::iterator it = SplitExonEnds.begin(); it != SplitExonEnds.end();++it)
		{
			CurrentRNA.ExonEnds.push_back(str2int(*it));
		}
		
		CurrentRNA.Score = str2int(SplitRefSeqLine[11]);
		CurrentRNA.GeneName = SplitRefSeqLine[12];
		CurrentRNA.CDSStartStat = SplitRefSeqLine[13];
		CurrentRNA.CDSEndStat = SplitRefSeqLine[14];
		
		std::vector<std::string> SplitExonFrames;
		tokenize(SplitRefSeqLine[15],SplitExonFrames,",");
		for (std::vector<std::string>::iterator it = SplitExonFrames.begin(); it != SplitExonFrames.end();++it)
		{
			if (*it == "-1") CurrentRNA.ExonFrames.push_back(-1);
			else CurrentRNA.ExonFrames.push_back(str2int(*it));
		}
		
		RefSeqTranscripts[ProcessChromosomeString(SplitRefSeqLine[2])].push_back(CurrentRNA);
	}
	std::clog << "/\\ Loaded RefSeqDB...\n";
}

#endif