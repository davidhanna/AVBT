#ifndef COMMONALGORITHMS_H
#define COMMONALGORITHMS_H

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

inline long int str2int( const std::string &str )
{
	register const char *p = str.data(), *pEnd = (str.data() + str.size());
	register long int i = 0;
	bool negative = false;
	while( p != pEnd )
	{
		if( *p == '-' ) negative = true;
		if( isdigit(*p) ) i = i * 10 + (*p++ - '0');
		else p++;
	}
	if (negative) return (i*(-1));
	else return i;
}

inline void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters = "\t", const bool trimEmpty = false)  
{
   std::string::size_type pos, lastPos = 0;
   while(true)
   {
      pos = str.find_first_of(delimiters, lastPos);
      if(pos == std::string::npos)
      {
         pos = str.length();

         if(pos != lastPos || !trimEmpty)
            tokens.push_back(std::vector<std::string>::value_type(str.data()+lastPos,
                  (std::vector<std::string>::value_type::size_type)pos-lastPos ));

         break;
      }
      else
      {
         if(pos != lastPos || !trimEmpty)
            tokens.push_back(std::vector<std::string>::value_type(str.data()+lastPos,
                  (std::vector<std::string>::value_type::size_type)pos-lastPos ));
      }

      lastPos = pos + 1;
   }
};

std::string ProcessChromosomeString(std::string Chromosome)
{
	if (Chromosome.find("chr") != std::string::npos) Chromosome.erase(0,3);
/*
	if (Chromosome.find('_') != std::string::npos)
	{
		std::vector<std::string> SplitChromosome;
		tokenize(Chromosome,SplitChromosome,"_");
		if (Chromosome.find("gl") != std::string::npos)
		{
			Chromosome = "GL" + SplitChromosome[1].substr(2) + ".1";
		}
		else
		{
			Chromosome = SplitChromosome[0];
		}
	}
*/
	return Chromosome;
	std::cout << Chromosome << "\n";
}

#endif