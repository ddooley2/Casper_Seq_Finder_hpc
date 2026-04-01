#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "CrisprGroup.h"
#include "Read.h"
#include "pameval.h"
#include "Write.h"
#include "sqlite3.h"

using namespace std;

int main(int argc, char* argv[])
//int main()
{
	clock_t tStart = clock();
	//vector<string> argv = { "Executable", "saCas9", "NNGRRT", "TRUE","FALSE","TRUE","4","16","0","bsu","C:/Users/Tfry/Desktop/","C:/Users/Tfry/Desktop/CASPERapp/CASPERinfo","C:/Users/Tfry/Desktop/Recommended_CSPR_Files/FNA/bacillus_subtilis_168.fna", "bsu", "notes_go_here", "DATA:CRISPRSCAN" };
	//vector<string> argv = { "Executable", "asCas12", "TTTV", "TRUE","TRUE","TRUE","0","16","8","sce","C:/Users/Tfry/Desktop/","C:/Users/Tfry/Desktop/CASPERapp/CASPERinfo","C:/Users/Tfry/Desktop/Recommended_CSPR_Files/FNA/sce.fna", "sce", "notes_go_here", "DATA:CRISPRSCAN" };
	//vector<string> argv = { "Executable", "saCas9", "NNGRRT", "TRUE","FALSE","TRUE","4","16","0","sce","C:/Users/Tfry/Desktop/","C:/Users/Tfry/Desktop/CASPERapp/CASPERinfo","C:/Users/Tfry/Desktop/Recommended_CSPR_Files/FNA/sce.fna", "sce", "notes_go_here", "DATA:CRISPRSCAN" };
	//vector<string> argv = { "Executable", "spCas9", "NGG", "TRUE","FALSE","TRUE","4","16","0","ac","C:/Users/Tfry/Desktop/","C:/Users/Tfry/Desktop/CASPERapp/CASPERinfo","C:/Users/Tfry/Desktop/Recommended_CSPR_Files/FNA/acinetobacter_baumanii.fna", "ac", "notes_go_here", "DATA:CRISPRSCAN" };
	//vector<string> argv = { "Executable", "spCas9", "NGG", "TRUE","FALSE","TRUE","4","16","0","gallus","C:/Users/Tfry/Desktop/","C:/Users/Tfry/Desktop/CASPERapp/CASPERinfo","C:/Users/Tfry/Desktop/Recommended_CSPR_Files/FNA/gallus_gallus.fna", "gallus", "notes_go_here", "DATA:CRISPRSCAN" };
	//vector<string> argv = { "Executable", "saCas9", "NNGRRT", "TRUE","FALSE","TRUE","4","16","0","gallus","C:/Users/Tfry/Desktop/","C:/Users/Tfry/Desktop/CASPERapp/CASPERinfo","C:/Users/Tfry/Desktop/Recommended_CSPR_Files/FNA/gallus_gallus.fna", "gallus", "notes_go_here", "DATA:CRISPRSCAN" };
	//vector<string> argv = { "Executable", "asCas12", "TTTV", "TRUE","TRUE","TRUE","0","20","4","ac","C:/Users/Tfry/Desktop/","C:/Users/Tfry/Desktop/CASPERapp/CASPERinfo","C:/Users/Tfry/Desktop/Recommended_CSPR_Files/FNA/acinetobacter_baumanii.fna", "ac", "notes_go_here", "DATA:CRISPRSCAN" };
	//vector<string> argv = { "Executable", "asCas12", "TTTV", "TRUE","TRUE","TRUE","0","20","4","gallus","C:/Users/Tfry/Desktop/","C:/Users/Tfry/Desktop/CASPERapp/CASPERinfo","C:/Users/Tfry/Desktop/Recommended_CSPR_Files/FNA/gallus_gallus.fna", "gallus", "notes_go_here", "DATA:CRISPRSCAN" };
	//vector<string> argv = { "Executable", "saCas9", "NNGRRT", "TRUE","FALSE","TRUE","4","16","0","ac","C:/Users/Tfry/Desktop/","C:/Users/Tfry/Desktop/CASPERapp/CASPERinfo","C:/Users/Tfry/Desktop/Recommended_CSPR_Files/FNA/acinetobacter_baumanii.fna", "ac", "notes_go_here", "DATA:CRISPRSCAN" };

	//vector<string> argv = { "Executable", "spCas9", "NGG", "TRUE","FALSE","TRUE","4","16","0","ac","C:/Users/Tfry/Desktop/","C:/Users/Tfry/Desktop/CASPERapp/CASPERinfo","C:/Users/Tfry/Desktop/Recommended_CSPR_Files/FNA/acinetobacter_baumanii.fna", "ac", "notes_go_here", "DATA:CRISPRSCAN" };
	//vector<string> argv = { "Executable", "saCas9", "NNGRRT", "TRUE","FALSE","TRUE","4","16","0","ac","C:/Users/Tfry/Desktop/","C:/Users/Tfry/Desktop/CASPERapp/CASPERinfo","C:/Users/Tfry/Desktop/Recommended_CSPR_Files/FNA/acinetobacter_baumanii.fna", "ac", "notes_go_here", "DATA:CRISPRSCAN" };
	//vector<string> argv = { "Executable", "asCas12", "TTTV", "TRUE","TRUE","TRUE","0","20","4","ac","C:/Users/Tfry/Desktop/","C:/Users/Tfry/Desktop/CASPERapp/CASPERinfo","C:/Users/Tfry/Desktop/Recommended_CSPR_Files/FNA/acinetobacter_baumanii.fna", "ac", "notes_go_here", "DATA:CRISPRSCAN" };


	//variables
	CrisprGroup genome;
	Read read;
	Write write;
	string endo, pam, org_code, output_path, score_file, input_file, org_name, notes, cspr_filename, db_filename, pam_regex, on_target_data;
	int five_length, seed_length, three_length, seq_length, pam_length, total = 0;
	bool directionality, reps;
	bool strand = true;
	bool mt = true;
	vector<string> chroms;
	vector<string> sequences;
	vector<unsigned long> compressed_seeds, kstats;
	vector<int> on_target_scores;
	vector<long> seed_locs;
	vector<int> seed_cnts;
	vector<int> uniques, repeats;

	//read in argv
	endo = string(argv[1]);
	pam = string(argv[2]);
	if (string(argv[3]) == "TRUE")
	{
		mt = true;
	}
	else
	{
		mt = false;
	}

	if (string(argv[4]) == "TRUE")
	{
		directionality = true;
	}
	else
	{
		directionality = false;
	}
	if (string(argv[5]) == "TRUE")
	{
		reps = true;
	}
	else
	{
		reps = false;
	}
	five_length = stoi(string(argv[6]));
	seed_length = stoi(string(argv[7]));
	three_length = stoi(string(argv[8]));
	org_code = string(argv[9]);
	output_path = string(argv[10]);
	score_file = string(argv[11]);
	input_file = string(argv[12]);
	org_name = string(argv[13]);
	notes = string(argv[14]);
	on_target_data = string(argv[15]);
	
	// Parse optional thread count parameter
	int num_threads = std::thread::hardware_concurrency();
	if (argc > 16) {
		num_threads = stoi(string(argv[16]));
	}
	if (num_threads <= 0) num_threads = 4;
	cout << "Using " << num_threads << " threads." << endl;
	
	seq_length = five_length + seed_length + three_length;
	pam_length = pam.size();
	total = seq_length + pam_length;

	//generate list of valid_pams
	pameval PamEval(pam);
	PamEval.generatePamsWrapper();

	//filenames
	cspr_filename = output_path + org_code + "_" + endo + ".cspr";
	db_filename = output_path + org_code + "_" + endo + "_repeats.db";

	//read input file
	cout << "Reading input file." << endl;
	read.read_file(input_file, sequences, chroms, kstats);
	seed_cnts.resize(sequences.size());

	cout << "Number of Chromosomes/Scaffolds: " << sequences.size() << endl;

	//find pams
	cout << "Finding Targets." << endl;
	genome.findPAMs(directionality, mt, sequences, compressed_seeds, seed_locs, on_target_scores, seed_cnts, strand, PamEval, pam_length, seq_length, five_length, seed_length, score_file, on_target_data, endo, directionality, pam);

	//process targets
	cout << "Processing Targets." << endl;
	genome.process_targets(uniques, repeats, compressed_seeds, seed_locs);

	//writing data
	if (directionality)
	{
		cout << "Extracting Uniques..." << endl;
		thread t1([&write, &uniques, &sequences, &seed_locs, &on_target_scores, &seed_cnts, &kstats, &org_name, &cspr_filename, &score_file, &chroms, &notes, &pam_length, &seq_length]() {write.write_uniques_dir(uniques, sequences, seed_locs, on_target_scores, seed_cnts, kstats, org_name, cspr_filename, score_file, chroms, notes, pam_length, seq_length); });

		cout << "Extracting Repeats..." << endl;
		thread t2([&write, &db_filename, &repeats, &sequences, &seed_locs, &compressed_seeds, &on_target_scores, &seed_cnts, &score_file, &five_length, &three_length, &seed_length, &pam_length, &seq_length]() {write.write_repeats_dir(db_filename, repeats, sequences, seed_locs, compressed_seeds, on_target_scores, seed_cnts, score_file, five_length, three_length, seed_length, pam_length, seq_length); });

		t1.join();
		cout << "Uniques extracted." << endl;
		t2.join();
		cout << "Repeats extracted." << endl;
	
	}
	else
	{
		cout << "Extracting Uniques..." << endl;
		thread t1([&write, &uniques, &sequences, &seed_locs, &on_target_scores, &seed_cnts, &kstats, &org_name, &cspr_filename, &score_file, &chroms, &notes, &pam_length, &seq_length]() {write.write_uniques(uniques, sequences, seed_locs, on_target_scores, seed_cnts, kstats, org_name, cspr_filename, score_file, chroms, notes, pam_length, seq_length); });
		
		cout << "Extracting Repeats..." << endl;
		thread t2([&write, &db_filename, &repeats, &sequences, &seed_locs, &on_target_scores, &compressed_seeds, &seed_cnts, &score_file, &five_length, &three_length, &seed_length, &pam_length, &seq_length]() {write.write_repeats(db_filename, repeats, sequences, seed_locs, on_target_scores, compressed_seeds, seed_cnts, score_file, five_length, three_length, seed_length, pam_length, seq_length); });
		
		t1.join();
		cout << "Uniques extracted." << endl;
		t2.join();
		cout << "Repeats extracted." << endl;
	}

	cout << "Finished." << endl;
	printf("Time taken: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
	//system("pause");
	return 0;
}