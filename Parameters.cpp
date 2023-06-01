#include "Parameters.h"

Parameters::Parameters() {

	//Simulation
	SimNr = 65;

	//Genome
	bp = 10000000;
	SNP_density = 0.0000415; // given parameters below, acheived approximate SNP density of Zhao et al. 2003
	shape = 0.25567091;
	rate = 237.3453398;// scale parameter equivalent to the rate parameter 0.00421327 inferred by Peter from the PacBio data//OLD:scale parameter equivalent to the rate parameter 1/0.00166 yielding approx 0.00083 SNPs/base and 10Mbp genome
	coalescent_dist = true;
	binom_dist = false;
	p_success = 0.00083;

	//Gene conversions (GC)
	GC_events = 100;
	distribution = "geometric"; // distribution of GC lengths: geometric, normal, uniform beta or exponential
	succes_p = 0.01;
	mean_normal = 100.0;
	var_normal = 10;
	min = 1;
	max = 10;
	alpha = 0.5;
	beta = 0.5;

}

Parameters::~Parameters() {
}

void Parameters::outPara(string name) {

	ofstream out;
	out.open(name.c_str());

	//Simulation
	out << "SimNr\t" << SimNr << endl;

	//Genome
	out << "bp\t" << bp << endl;
	out << "SNP_density\t" << SNP_density << endl;
	out << "rate\t" << rate << endl;
	out << "shape\t" << shape << endl;
	out << "realized_genome_size\t" << realized_genome_size << endl;
	out << "coalescent_dist\t" << coalescent_dist << endl;
	out << "binom_dist\t" << binom_dist << endl;
	out << "p_success\t" << p_success << endl;

	//Gene conversions
	out << "GC_events\t" << GC_events << endl;
	out << "distribution\t" << distribution << endl;
	out << "succes_p\t" << succes_p << endl;
	out << "mean_normal\t" << mean_normal << endl;
	out << "var_normal\t" << var_normal << endl;
	out << "min\t" << min << endl;
	out << "max\t" << max << endl;
	out << "alpha\t" << alpha << endl;
	out << "beta\t" << beta << endl;

	out.close();
}

