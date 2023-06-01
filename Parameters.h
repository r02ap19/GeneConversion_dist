#pragma once

#include <fstream>

using namespace std;

class Parameters
{
public:
	Parameters();
	~Parameters();
	//Simulation
	int SimNr;

	//Genome
	int bp; //Number of base pairs in haplotype
	double SNP_density; //Density of SNPs
	double rate; //rate parameter for gamma distribution of inter SNP distances (see Field et al. 2016).
	double shape; //shape parameter for gamma distribution of inter SNP distances (see Field et al. 2016).
	bool coalescent_dist; // if true, distances between SNPs are gamma distributed, approximating the expectation from coalescent theory (see Field et al. 2016).
	int realized_genome_size; //size of the genome when genome is is determined by inter-SNP distances i.e. when coalescent_dist == true. Else, realized_genome_size = bp.
	bool binom_dist; //if true, each base has a specific probability of being a SNP (p_success) such that the distribution of SNPs (and inter-SNP distances) follow a binomial dist.
	double p_success; //per base probability of being a SNP

	//Gene conversions (GC)
	int GC_events; // number of gene conversion events (fixed for now, but could be sampled).
	string distribution; // distribution of GC locations: geometric, normal, uniform, beta or exponential
	double succes_p; // probability of success, geometric distribution
	double mean_normal; // mean of normal distribution
	double var_normal; // variance, normal distribution
	int min; //lower bound for uniformly distributed GC tract lengths
	int max; //upper bound for uniformaly distributed GC tract lengths
	double alpha; //alpha parameter in beta distribution for beta distributed GC tract lengths
	double beta; //beta parameter in beta distribution for beta distributed GC tract lengths

	void outPara(string dir); //Parameter output
};


