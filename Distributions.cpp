#include "Distributions.h"

#if GenomeDK
int main(int argc, char* argv[])
{
	// Get the current directory.
	char* buffer = getcwd(NULL, 0);
	dir = buffer;
	free(buffer);
	dir = dir + "/"; //Current directory path
	dirOut = dir + "Outputs/"; //Outpus folder path

	para.SimNr = std::atoi(argv[1]);
	para.rate = std::atof(argv[2]);

	Simulate();
	cout << "Simulation completed" << endl;
	return 0;
}
#else
int _tmain(int argc, _TCHAR* argv[])
{
	// Get the current directory.
	char* buffer = _getcwd(NULL, 0);
	dir = buffer;
	free(buffer);
	dir = dir + "\\"; //Current directory path
	dirOut = dir + "Outputs\\"; //Outpus folder path
	extime = clock();

	Simulate();
	std::cout << "Simulation completed" << endl;
	return 0;
}
#endif
const string Int2Str(const int x)
{
	ostringstream o;
	if (!(o << x)) return "ERROR";
	return o.str();
}

void Simulate(void) {
	std::cout << "Simulation nr. " << para.SimNr << endl;
	string name = dirOut + "Sim" + Int2Str(para.SimNr) + "_Para.txt";

	outSeq_header();
	outConv_header();
	cout << "Initializing...";
	initialize();
	cout << "done" <<  endl;

	extime = clock() - extime;
	std::cout << "time = " << (float)extime / CLOCKS_PER_SEC << " sec" << endl;
	extime = clock();
	cout << "\n" << endl;

	cout << "Adding markers...";
	add_SNPs();
	cout << "done" << endl;

	extime = clock() - extime;
	std::cout << "time = " << (float)extime / CLOCKS_PER_SEC << " sec" << endl;
	extime = clock();
	cout << "\n" << endl;

	cout << "Conducting gene conversion events...";
	add_GeneConversions();
	cout << "done" << endl;

	extime = clock() - extime;
	std::cout << "time = " << (float)extime / CLOCKS_PER_SEC << " sec" << endl;
	extime = clock();
	cout << "\n" << endl;

	para.outPara(name);
	if (seq.is_open()) seq.close();
	if (conv.is_open()) conv.close();
}

void initialize(void) {
	Haplotype H1, H2;
	int base_value, base_counter, distance;
	std::map<int, char>::iterator iter;
	base_counter = 0;
	std::gamma_distribution<> dists(para.shape, para.rate);

	if (para.coalescent_dist == true) {
		int N_SNPs = round(para.SNP_density * para.bp);
		for (int x = 0; x < N_SNPs; x++) {
			distance = round(dists(rdgen));
			if (distance < 1) {
				distance = 1;
			}
			cout << distance << endl;
			base_counter = base_counter + distance;
			dataset.polymorphic_positions.push_back(base_counter);
		}
		para.realized_genome_size = base_counter;
	} else {
		para.realized_genome_size = para.bp;
	}

	for (int x = 0; x < para.realized_genome_size; x++) {
		base_value = WatsonCrick_bases(rdgen);
		if (base_value == 0)(H1.sequence[x] = 'A');
		if (base_value == 1)(H1.sequence[x] = 'T');
		if (base_value == 2)(H1.sequence[x] = 'C');
		if (base_value == 3)(H1.sequence[x] = 'G');
	}

	H2 = H1;
	H2.number = 2;

	dataset.tmp_haps.push_back(H1);
	dataset.tmp_haps.push_back(H2);
}

void add_SNPs(void) {
	int N_SNPs = round(para.SNP_density * para.bp);
	int pos;
	double prob;
	std::map<int, char>::iterator iter;
	std::vector<Haplotype>::iterator iter2;

	std::uniform_int_distribution<> genome_pos(0, para.realized_genome_size);
	std::uniform_real_distribution<> SNP_prob(0.0, 1.0);
	if (para.coalescent_dist == true) {
		for (int x = 0; x < N_SNPs; x++) {
			dataset.tmp_haps[0].sequence[dataset.polymorphic_positions[x]] = 'X';
			dataset.tmp_haps[0].SNPs[dataset.polymorphic_positions[x]] = 'X';
		}
	}
	if (para.binom_dist == true) {
		for (int x = 0; x < para.realized_genome_size; x++) {
			prob = SNP_prob(rdgen);
			if (prob < para.p_success) {
				dataset.tmp_haps[0].sequence[x] = 'X';
				dataset.tmp_haps[0].SNPs[x] = 'X';
			}
		}
	} else {
		pos = genome_pos(rdgen);

		for (int x = 0; x < N_SNPs; x++) {
			pos = genome_pos(rdgen);
			while (dataset.tmp_haps[0].sequence[pos] == 'X') {
				pos = genome_pos(rdgen);
			}
			dataset.tmp_haps[0].sequence[pos] = 'X';
			dataset.tmp_haps[0].SNPs[pos] = 'X';
		}
	}

	dataset.haps = dataset.tmp_haps;
}

//GC function
void add_GeneConversions() {
	int pos, GC_tractLength;
	std::map<int, char>::iterator iter;
	std::map<int, tract>::iterator iter2;
	//std::map<int, tract>::iterator iter3;
	tract event;
	std::uniform_int_distribution<> genome_pos(0, para.realized_genome_size);

	for (int x = 0; x < para.GC_events; x++) {
		if (para.distribution == "geometric") {
			GC_tractLength = round(geo(rdgen));
		}
		if (para.distribution == "normal") {
			GC_tractLength = round(normal(rdgen));
			cout << "test" << endl;
		}
		if (para.distribution == "uniform") {
			GC_tractLength = round(uni(rdgen));
		}
		//if (para.distribution == "beta") {
		//	GC_tractLength = round(beta(rdgen));
		//}
		pos = genome_pos(rdgen);
		event.start = pos;
		event.length = GC_tractLength;
		dataset.haps[1].GC_events[x] = event;

		for (int y = 0; y < GC_tractLength; y++) {
			dataset.haps[1].sequence[pos + y] = dataset.haps[0].sequence[pos + y];
		}
	}

	cout << "done" << endl;

	extime = clock() - extime;
	std::cout << "time = " << (float)extime / CLOCKS_PER_SEC << " sec" << endl;
	extime = clock();
	cout << "\n" << endl;

	cout << "Writing output files...";

	dataset.converted_haps = dataset.haps;

	for (iter2 = dataset.converted_haps[1].GC_events.begin(); iter2 != dataset.converted_haps[1].GC_events.end(); iter2++) {
		outConv(iter2->second.start, iter2->second.length, iter2->first, &conv);
	}
	for (iter = dataset.converted_haps[0].sequence.begin(); iter != dataset.converted_haps[0].sequence.end(); iter++) {
		outSeq(dataset.converted_haps[0].number, iter->second, &seq);
	}
	for (iter = dataset.converted_haps[1].sequence.begin(); iter != dataset.converted_haps[1].sequence.end(); iter++) {
		outSeq(dataset.converted_haps[1].number, iter->second, &seq);
	}
}

Data::Data() {
	holder = 0;
}
Data::~Data() {
}

void outSeq_header(void) {
	string name;
	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_outSeq.txt";
	seq.open(name.c_str());
	seq << "Haplotype\tseq" << endl;
}

void outConv_header(void) {
	string name;
	name = dirOut + "Sim" + Int2Str(para.SimNr) + "_outConv.txt";
	conv.open(name.c_str());
	conv << "pos\ttract_length\tID" << endl;
}

void outSeq(int h_number, char base, std::ofstream* out) {
	*out << h_number << "\t" << base << "\t";
	*out << endl;
}

void outConv(int pos, int tract_length, int ID, std::ofstream* out) {
	*out << pos << "\t" << tract_length << "\t" << ID << "\t";
	*out << endl;
}

