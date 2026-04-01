#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>
#include "Score.h"
#include "sqlite3.h"
#include <math.h>
#include "pameval.h"



using namespace std;

class Write
{
public:
	void write_uniques(vector<int> &uniques, vector<string> &sequences, vector<long> &seed_locs, vector<int> &on_target_scores, vector<int> &seed_cnts, vector<unsigned long> &kstats, string &org_name, string &filename, string &score_file, vector<string> &chroms, string &notes, int &pam_length, int &seq_length);
	void write_repeats(string &filename, vector<int> &repeats, vector<string> &sequences, vector<long> &seed_locs, vector<int> &on_target_scores, vector<unsigned long> &comp_seeds, vector<int> &seed_cnts, string &score_file, int &five_length, int &three_length, int &seed_length, int &pam_length, int &seq_length);
	void write_uniques_dir(vector<int> &uniques, vector<string> &sequences, vector<long> &seed_locs, vector<int> &on_target_scores, vector<int> &seed_cnts, vector<unsigned long> &kstats, string &org_name, string &filename, string &score_file, vector<string> &chroms, string &notes, int &pam_length, int &seq_length);
	void write_repeats_dir(string &filename, vector<int> &repeats, vector<string> &sequences, vector<long> &seed_locs, vector<unsigned long> &compressed_seeds, vector<int> &on_target_scores, vector<int> &seed_cnts, string &score_file, int &five_length, int &three_length, int &seed_length, int &pam_length, int &seq_length);
	int convertCharBase4(char &c);
	string reverseComplement(string &str);
	int findChromosome(int idx, vector<int> &boundaries);
};