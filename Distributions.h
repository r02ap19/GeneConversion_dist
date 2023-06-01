#pragma once

#define GenomeDK 0

#include <stdio.h>
#include <stdlib.h>
#if GenomeDK 
#include <unistd.h>
#else
#include <tchar.h> 
#include <direct.h>
#include <io.h>
#endif
#include <vector>
#include <set>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <numeric>
#include <time.h>
#include <random>
#include <iterator>

#include "Parameters.h"
#include "Haplotype.h"

using namespace std;

class Data {
public:
	Data();
	~Data();
	int holder;
	vector<Haplotype> haps;
	vector<Haplotype> tmp_haps;
	vector<Haplotype> converted_haps;
	vector<int> polymorphic_positions;

private:
};

clock_t extime;

Parameters para;
Data dataset;
string dir, dirOut;
ofstream seq, conv;

std::random_device rd;
std::mt19937 rdgen(rd());

std::geometric_distribution<> geo(para.succes_p);
std::normal_distribution<> normal(para.mean_normal, para.var_normal);
std::uniform_int_distribution<> uni(para.min, para.max);
//std::_Beta_distribution<> beta(para.alpha, para.beta);
std::uniform_int_distribution<> WatsonCrick_bases(0, 3);


// Functions declaration
const string Int2Str(const int x);
void Simulate(void);
void initialize(void);
void add_SNPs(void);
void add_GeneConversions(void);
void outSeq_header(void);
void outConv_header(void);
void outSeq(int h_number, char base, std::ofstream* out);
void outConv(int pos, int tract_length, int ID, std::ofstream* out);
