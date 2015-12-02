#ifndef SETTINGS_H
#define SETTINGS_H

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

struct ProgramOptions
{
	ProgramOptions();
	std::string RefSeqPath;
	std::vector<std::pair<std::string,bool> > VCFFilePaths;
	std::string PedigreePath;
	std::string BEDPath;
	float FrequenciesCutoff;
	float QualityCutoff;
	int ArtifactCutoff;
};
ProgramOptions::ProgramOptions()
{
	FrequenciesCutoff = 1.0;
	QualityCutoff = 1.0;
	ArtifactCutoff = 0;
}


#endif