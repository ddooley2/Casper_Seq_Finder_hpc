#include "Write.h"
#include <sstream>
using namespace std;

//abs comparison used for sorting vector of longs
bool abs_cmp(long i1, long i2)
{
	return (abs(i1) < abs(i2));
}


void Write::write_uniques(vector<int> &uniques, vector<string> &sequences, vector<long> &seed_locs, vector<int> &on_target_scores, vector<int> &seed_cnts, vector<unsigned long> &kstats, string &org_name, string &filename, string &score_file, vector<string> &chroms, string &notes, int &pam_length, int &seq_length)
{
	//variables
	string comp, seq, genome, kstat, misc, full_seq, curr_pam;
	long pos = 0;
	int i = 0;
	int j = 0;
	int running_cnt = 0;
	int leftover_padding = 35 - 6 - seq_length - pam_length;

	FILE* outputfile;
	outputfile = fopen(filename.c_str(),"w");
	
	char *file_buffer = new char[1024 * 1024];
	setvbuf(outputfile, file_buffer, _IOFBF, 1024 * 1024);

	//output first details - org name, kstats, misc
	genome = "GENOME: " + org_name;
	kstat = "KARYSTATS: ";
	for (int i = 0; i < kstats.size(); i++)
	{
		kstat += to_string(kstats[i]) + ",";
	}
	misc = "MISCELLANEOUS: " + notes;

	fprintf(outputfile, "%s\n", genome.c_str());
	fprintf(outputfile, "%s\n", kstat.c_str());
	fprintf(outputfile, "%s\n", misc.c_str());

	string ss_temp;
	stringstream ss(ss_temp);
	vector<int> temp_indexes;
	//loop through all chromosomes, get locations of pams, write out pams
	for (int curr_chrom = 0; curr_chrom < sequences.size(); curr_chrom++)
	{
		fprintf(outputfile, "%s\n", chroms[curr_chrom].c_str());
		
		//get locations pertaining to pams found in the current chromosome
		running_cnt += seed_cnts[curr_chrom];
		while (true)
		{
			if (j >= uniques.size())
			{
				break;
			}
			if (uniques[j] >= running_cnt)
			{
				break;
			}
			temp_indexes.push_back(uniques[j]);
			j++;
		}
		//sort locations based on absolute value
		sort(temp_indexes.begin(), temp_indexes.end(),
			[&seed_locs](int A, int B) -> bool {
			return abs(seed_locs[A]) < abs(seed_locs[B]);
		});

		//store reverse complement of current chromosome
		comp = reverseComplement(sequences[curr_chrom]);

		//loop through locations in current chromosomes, extract sequence, calculate score, write out
		for (int i = 0; i < temp_indexes.size(); i++)
		{
			if (seed_locs[temp_indexes[i]] > 0)
			{
				seq = sequences[curr_chrom].substr(seed_locs[temp_indexes[i]] - seq_length, seq_length + pam_length);
				curr_pam = seq.substr(seq_length, pam_length);
				seq = seq.substr(0, seq_length);
				fprintf(outputfile, "%ld,%s,%s,%i\n", seed_locs[temp_indexes[i]], seq.c_str(), curr_pam.c_str(), on_target_scores[temp_indexes[i]]);
			}
			else
			{
				pos = comp.size() + seed_locs[temp_indexes[i]] + 1;
				seq = comp.substr(pos - seq_length, seq_length + pam_length);
				curr_pam = seq.substr(seq_length, pam_length);
				seq = seq.substr(0, seq_length);
				fprintf(outputfile, "%ld,%s,%s,%i\n", seed_locs[temp_indexes[i]], seq.c_str(), curr_pam.c_str(), on_target_scores[temp_indexes[i]]);
			}
		}
		//clear temporary locations vector
		temp_indexes.clear();
	}
	//close file
	fclose(outputfile);
	delete[] file_buffer;

	//cleanup - clear out uniques vector
	uniques.clear();
	uniques.shrink_to_fit();
}


void Write::write_repeats(string& filename, vector<int> &repeats, vector<string> &sequences, vector<long> &seed_locs, vector<int> &on_target_scores, vector<unsigned long> &compressed_seeds, vector<int> &seed_cnts, string &score_file, int &five_length, int &three_length, int &seed_length, int &pam_length, int &seq_length)
{
	//variables
	sqlite3 *db;
	sqlite3_stmt *stmt;
	char *zErrMsg = 0;
	string seq, locs, seed, scores, fives, threes, pams, cs, full_seq;
	int leftover_padding = 35 - 6 - seq_length - pam_length;
	int rc;
	int i = 0;
	int curr_chrom = 0;
	int cnt = 0;
	long pos = 0;
	vector<string> comps(sequences.size());
	char num_buf[32];

	// Build chromosome boundaries for O(1) lookup
	vector<int> chrom_boundaries(seed_cnts.size() + 1);
	chrom_boundaries[0] = 0;
	for (size_t j = 0; j < seed_cnts.size(); j++)
	{
		chrom_boundaries[j + 1] = chrom_boundaries[j] + seed_cnts[j];
	}

	//open and setup DB with aggressive optimizations
	rc = sqlite3_open(filename.c_str(), &db);
	rc = sqlite3_exec(db, "PRAGMA journal_mode = OFF;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "PRAGMA synchronous = OFF;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "PRAGMA locking_mode = EXCLUSIVE;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "PRAGMA temp_store = MEMORY;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "PRAGMA cache_size = -512000;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "PRAGMA mmap_size = 1073741824;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "PRAGMA page_size = 65536;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "DROP TABLE IF EXISTS repeats;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "VACUUM;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "CREATE TABLE repeats (seed TEXT, chromosome TEXT, location TEXT, three TEXT, five TEXT, pam TEXT, score TEXT, count INT);", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "BEGIN TRANSACTION;", NULL, 0, &zErrMsg);
	
	//prepare statement for inserts
	rc = sqlite3_prepare_v2(db, "INSERT INTO repeats (seed, chromosome, location, three, five, pam, score, count) VALUES (?,?,?,?,?,?,?,?)", -1, &stmt, NULL);
	
	// Pre-allocate string buffers
	locs.reserve(10000);
	cs.reserve(5000);
	threes.reserve(5000);
	fives.reserve(5000);
	pams.reserve(5000);
	scores.reserve(10000);
	
	//get all reverse comps of chromosome sequences
	for (i = 0; i < sequences.size(); i++)
	{
		comps[i] = reverseComplement(sequences[i]);
	}
	i = 0;
	//loop through repeats indices (they should be grouped together, so loop until compressed seed value isnt equal anymore)
	while (i < repeats.size())
	{
		// Clear and reuse string buffers
		locs.clear();
		cs.clear();
		threes.clear();
		fives.clear();
		pams.clear();
		scores.clear();
		
		// Find chromosome using binary search on boundaries
		curr_chrom = findChromosome(repeats[i], chrom_boundaries);
		
		int len = snprintf(num_buf, sizeof(num_buf), "%ld", seed_locs[repeats[i]]);
		locs.append(num_buf, len);

		if (seed_locs[repeats[i]] > 0)
		{
			seq = sequences[curr_chrom].substr(seed_locs[repeats[i]] - seq_length, seq_length + pam_length);
			full_seq = sequences[curr_chrom].substr(seed_locs[repeats[i]] - seq_length - leftover_padding, 35);
		}
		else
		{
			pos = comps[curr_chrom].size() + seed_locs[repeats[i]] + 1;
			seq = comps[curr_chrom].substr(pos - seq_length, seq_length + pam_length);
			full_seq = comps[curr_chrom].substr(pos - seq_length - leftover_padding, 35);
		}
		
		len = snprintf(num_buf, sizeof(num_buf), "%d", curr_chrom + 1);
		cs.append(num_buf, len);
		threes = seq.substr(five_length + seed_length, three_length);
		fives = seq.substr(0, five_length);
		pams = seq.substr(seq_length, pam_length);
		len = snprintf(num_buf, sizeof(num_buf), "%d", on_target_scores[repeats[i]]);
		scores.append(num_buf, len);
		cnt = 1;
		seed = seq.substr(five_length, seed_length);

		while (true)
		{
			if (i >= repeats.size() - 1)
			{
				i++;
				break;
			}
			if (compressed_seeds[repeats[i]] != compressed_seeds[repeats[i + 1]])
			{
				i++;
				break;
			}
			cnt++;
			
			len = snprintf(num_buf, sizeof(num_buf), ",%ld", seed_locs[repeats[i + 1]]);
			locs.append(num_buf, len);
			
			curr_chrom = findChromosome(repeats[i + 1], chrom_boundaries);
			
			if (seed_locs[repeats[i + 1]] > 0)
			{
				seq = sequences[curr_chrom].substr(seed_locs[repeats[i + 1]] - seq_length, seq_length + pam_length);
				full_seq = sequences[curr_chrom].substr(seed_locs[repeats[i + 1]] - seq_length - leftover_padding, 35);
			}
			else
			{
				pos = comps[curr_chrom].size() + seed_locs[repeats[i + 1]] + 1;
				seq = comps[curr_chrom].substr(pos - seq_length, seq_length + pam_length);
				full_seq = comps[curr_chrom].substr(pos - seq_length - leftover_padding, 35);
			}
			
			len = snprintf(num_buf, sizeof(num_buf), ",%d", curr_chrom + 1);
			cs.append(num_buf, len);
			
			threes += ",";
			threes += seq.substr(five_length + seed_length, three_length);
			fives += ",";
			fives += seq.substr(0, five_length);
			pams += ",";
			pams += seq.substr(seq_length, pam_length);
			
			len = snprintf(num_buf, sizeof(num_buf), ",%d", on_target_scores[repeats[i + 1]]);
			scores.append(num_buf, len);
			i++;
		}

		//bind parameters with SQLITE_STATIC and execute prepared statement
		sqlite3_bind_text(stmt, 1, seed.c_str(), -1, SQLITE_STATIC);
		sqlite3_bind_text(stmt, 2, cs.c_str(), -1, SQLITE_STATIC);
		sqlite3_bind_text(stmt, 3, locs.c_str(), -1, SQLITE_STATIC);
		sqlite3_bind_text(stmt, 4, threes.c_str(), -1, SQLITE_STATIC);
		sqlite3_bind_text(stmt, 5, fives.c_str(), -1, SQLITE_STATIC);
		sqlite3_bind_text(stmt, 6, pams.c_str(), -1, SQLITE_STATIC);
		sqlite3_bind_text(stmt, 7, scores.c_str(), -1, SQLITE_STATIC);
		sqlite3_bind_int(stmt, 8, cnt);
		sqlite3_step(stmt);
		sqlite3_reset(stmt);
	}

	//finalize statement, end db transaction, close db file
	sqlite3_finalize(stmt);
	rc = sqlite3_exec(db, "END TRANSACTION;", NULL, 0, &zErrMsg);
	sqlite3_close(db);
}


void Write::write_uniques_dir(vector<int> &uniques, vector<string> &sequences, vector<long> &seed_locs, vector<int> &on_target_scores, vector<int> &seed_cnts, vector<unsigned long> &kstats, string &org_name, string &filename, string &score_file, vector<string> &chroms, string &notes, int &pam_length, int &seq_length)
{
	//variables
	string comp, seq, genome, kstat, misc, full_seq, curr_pam;
	long pos = 0;
	int i = 0;
	int j = 0;
	int running_cnt = 0;
	int leftover_padding = 35 - 6 - seq_length - pam_length;

	FILE* outputfile;
	outputfile = fopen(filename.c_str(), "w");
	
	char *file_buffer = new char[1024 * 1024];
	setvbuf(outputfile, file_buffer, _IOFBF, 1024 * 1024);

	//output first details - org name, kstats, misc
	genome = "GENOME: " + org_name;
	kstat = "KARYSTATS: ";
	for (int i = 0; i < kstats.size(); i++)
	{
		kstat += to_string(kstats[i]) + ",";
	}
	misc = "MISCELLANEOUS: " + notes;

	fprintf(outputfile, "%s\n", genome.c_str());
	fprintf(outputfile, "%s\n", kstat.c_str());
	fprintf(outputfile, "%s\n", misc.c_str());

	string ss_temp;
	stringstream ss(ss_temp);
	vector<int> temp_indexes;
	//loop through all chromosomes, get locations of pams, write out pams
	for (int curr_chrom = 0; curr_chrom < sequences.size(); curr_chrom++)
	{
		fprintf(outputfile, "%s\n", chroms[curr_chrom].c_str());

		//get locations pertaining to pams found in the current chromosome
		running_cnt += seed_cnts[curr_chrom];
		while (true)
		{
			if (j >= uniques.size())
			{
				break;
			}
			if (uniques[j] >= running_cnt)
			{
				break;
			}
			temp_indexes.push_back(uniques[j]);
			j++;
		}
		//sort locations based on absolute value
		sort(temp_indexes.begin(), temp_indexes.end(),
			[&seed_locs](int A, int B) -> bool {
			return abs(seed_locs[A]) < abs(seed_locs[B]);
		});

		//store reverse complement of current chromosome
		comp = reverseComplement(sequences[curr_chrom]);

		//loop through locations in current chromosomes, extract sequence, calculate score, write out
		for (int i = 0; i < temp_indexes.size(); i++)
		{
			if (seed_locs[temp_indexes[i]] > 0)
			{
				seq = sequences[curr_chrom].substr(seed_locs[temp_indexes[i]] - 1, seq_length + pam_length);
				curr_pam = seq.substr(0, pam_length);
				seq = seq.substr(pam_length, seq_length);
				fprintf(outputfile, "%ld,%s,%s,%i\n", seed_locs[temp_indexes[i]] + pam_length, seq.c_str(), curr_pam.c_str(), on_target_scores[temp_indexes[i]]);
			}
			else
			{
				pos = comp.size() + seed_locs[temp_indexes[i]];
				seq = comp.substr(pos, seq_length + pam_length);
				curr_pam = seq.substr(0, pam_length);
				seq = seq.substr(pam_length, seq_length);
				fprintf(outputfile, "%ld,%s,%s,%i\n", seed_locs[temp_indexes[i]] + pam_length, seq.c_str(), curr_pam.c_str(), on_target_scores[temp_indexes[i]]);
			}
		}
		//clear temporary locations vector
		temp_indexes.clear();
	}
	//close file
	fclose(outputfile);
	delete[] file_buffer;

	//cleanup - clear out uniques vector
	uniques.clear();
	uniques.shrink_to_fit();
}


void Write::write_repeats_dir(string &filename, vector<int> &repeats, vector<string> &sequences, vector<long> &seed_locs, vector<unsigned long> &compressed_seeds, vector<int> &on_target_scores, vector<int> &seed_cnts, string &score_file, int &five_length, int &three_length, int &seed_length, int &pam_length, int &seq_length)
{
	//variables
	sqlite3 *db;
	sqlite3_stmt *stmt;
	char *zErrMsg = 0;
	string seq, locs, seed, scores, fives, threes, pams, cs, full_seq;
	int leftover_padding = 35 - 6 - seq_length - pam_length;
	int rc;
	int i = 0;
	int curr_chrom = 0;
	long pos = 0;
	int cnt = 0;
	vector<string> comps(sequences.size());
	char num_buf[32];

	// Build chromosome boundaries for O(1) lookup
	vector<int> chrom_boundaries(seed_cnts.size() + 1);
	chrom_boundaries[0] = 0;
	for (size_t j = 0; j < seed_cnts.size(); j++)
	{
		chrom_boundaries[j + 1] = chrom_boundaries[j] + seed_cnts[j];
	}

	//open and setup DB with aggressive optimizations
	rc = sqlite3_open(filename.c_str(), &db);
	rc = sqlite3_exec(db, "PRAGMA journal_mode = OFF;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "PRAGMA synchronous = OFF;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "PRAGMA locking_mode = EXCLUSIVE;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "PRAGMA temp_store = MEMORY;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "PRAGMA cache_size = -512000;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "PRAGMA mmap_size = 1073741824;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "PRAGMA page_size = 65536;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "DROP TABLE IF EXISTS repeats;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "VACUUM;", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "CREATE TABLE repeats (seed TEXT, chromosome TEXT, location TEXT, three TEXT, five TEXT, pam TEXT, score TEXT, count INT);", NULL, 0, &zErrMsg);
	rc = sqlite3_exec(db, "BEGIN TRANSACTION;", NULL, 0, &zErrMsg);

	//prepare statement for inserts
	rc = sqlite3_prepare_v2(db, "INSERT INTO repeats (seed, chromosome, location, three, five, pam, score, count) VALUES (?,?,?,?,?,?,?,?)", -1, &stmt, NULL);

	// Pre-allocate string buffers
	locs.reserve(10000);
	cs.reserve(5000);
	threes.reserve(5000);
	fives.reserve(5000);
	pams.reserve(5000);
	scores.reserve(10000);

	//get all reverse comps of chromosome sequences
	for (i = 0; i < sequences.size(); i++)
	{
		comps[i] = reverseComplement(sequences[i]);
	}
	i = 0;

	//loop through repeats indices (they should be grouped together, so loop until compressed seed value isnt equal anymore)
	while (i < repeats.size())
	{
		// Clear and reuse string buffers
		locs.clear();
		cs.clear();
		threes.clear();
		fives.clear();
		pams.clear();
		scores.clear();
		
		// Find chromosome using binary search on boundaries
		curr_chrom = findChromosome(repeats[i], chrom_boundaries);
		
		int len = snprintf(num_buf, sizeof(num_buf), "%ld", seed_locs[repeats[i]] + pam_length);
		locs.append(num_buf, len);
		
		if (seed_locs[repeats[i]] > 0)
		{
			seq = sequences[curr_chrom].substr(seed_locs[repeats[i]] - 1, seq_length + pam_length);
			full_seq = sequences[curr_chrom].substr(seed_locs[repeats[i]] - 1 - 6, 35);
		}
		else
		{
			pos = comps[curr_chrom].size() + seed_locs[repeats[i]];
			seq = comps[curr_chrom].substr(pos, seq_length + pam_length);
			full_seq = comps[curr_chrom].substr(pos - 6, 35);
		}
		
		len = snprintf(num_buf, sizeof(num_buf), "%d", curr_chrom + 1);
		cs.append(num_buf, len);
		threes = seq.substr(pam_length + five_length + seed_length, three_length);
		fives = seq.substr(pam_length, five_length);
		pams = seq.substr(0, pam_length);
		len = snprintf(num_buf, sizeof(num_buf), "%d", on_target_scores[repeats[i]]);
		scores.append(num_buf, len);
		cnt = 1;
		
		while (true)
		{
			if (compressed_seeds[repeats[i]] != compressed_seeds[repeats[i + 1]] || i >= repeats.size())
			{
				i++;
				break;
			}
			cnt++;
			
			len = snprintf(num_buf, sizeof(num_buf), ",%ld", seed_locs[repeats[i + 1]] + pam_length);
			locs.append(num_buf, len);
			
			curr_chrom = findChromosome(repeats[i + 1], chrom_boundaries);
			
			if (seed_locs[repeats[i + 1]] > 0)
			{
				seq = sequences[curr_chrom].substr(seed_locs[repeats[i + 1]] - 1, seq_length + pam_length);
				full_seq = sequences[curr_chrom].substr(seed_locs[repeats[i + 1]] - 1 - 6, 35);
			}
			else
			{
				pos = comps[curr_chrom].size() + seed_locs[repeats[i + 1]];
				seq = comps[curr_chrom].substr(pos, seq_length + pam_length);
				full_seq = comps[curr_chrom].substr(pos - 6, 35);
			}
			
			len = snprintf(num_buf, sizeof(num_buf), ",%d", curr_chrom + 1);
			cs.append(num_buf, len);
			
			threes += ",";
			threes += seq.substr(pam_length + five_length + seed_length, three_length);
			fives += ",";
			fives += seq.substr(pam_length, five_length);
			pams += ",";
			pams += seq.substr(0, pam_length);
			
			len = snprintf(num_buf, sizeof(num_buf), ",%d", on_target_scores[repeats[i + 1]]);
			scores.append(num_buf, len);
			i++;
		}

		//bind parameters with SQLITE_STATIC and execute prepared statement
		seed = seq.substr(pam_length + five_length, seed_length);
		sqlite3_bind_text(stmt, 1, seed.c_str(), -1, SQLITE_STATIC);
		sqlite3_bind_text(stmt, 2, cs.c_str(), -1, SQLITE_STATIC);
		sqlite3_bind_text(stmt, 3, locs.c_str(), -1, SQLITE_STATIC);
		sqlite3_bind_text(stmt, 4, threes.c_str(), -1, SQLITE_STATIC);
		sqlite3_bind_text(stmt, 5, fives.c_str(), -1, SQLITE_STATIC);
		sqlite3_bind_text(stmt, 6, pams.c_str(), -1, SQLITE_STATIC);
		sqlite3_bind_text(stmt, 7, scores.c_str(), -1, SQLITE_STATIC);
		sqlite3_bind_int(stmt, 8, cnt);
		sqlite3_step(stmt);
		sqlite3_reset(stmt);
	}

	//finalize statement, end db transaction, close db file
	sqlite3_finalize(stmt);
	rc = sqlite3_exec(db, "END TRANSACTION;", NULL, 0, &zErrMsg);
	sqlite3_close(db);
}


int Write::convertCharBase4(char &c) {
	switch (c) {
	case 'A': return 0;
	case 'T': return 1;
	case 'C': return 2;
	case 'G': return 3;
	default: return 0;
	}
}

string Write::reverseComplement(string &str)
{
	string rc;
	rc.reserve(str.size());
	for (long i = str.size() - 1; i >= 0; i--)
	{
		char n = str[i];
		char reverse;
		switch (n) {
		case 'A': reverse = 'T'; break;
		case 'T': reverse = 'A'; break;
		case 'G': reverse = 'C'; break;
		case 'C': reverse = 'G'; break;
		default: reverse = 'N';
		}
		rc += reverse;
	}
	return rc;
}

int Write::findChromosome(int idx, vector<int> &boundaries)
{
	int lo = 0, hi = boundaries.size() - 2;
	while (lo < hi)
	{
		int mid = (lo + hi + 1) / 2;
		if (boundaries[mid] <= idx)
			lo = mid;
		else
			hi = mid - 1;
	}
	return lo;
}
