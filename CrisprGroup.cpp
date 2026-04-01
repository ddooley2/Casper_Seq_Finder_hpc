#include "CrisprGroup.h"

using namespace std;

/* Constructor: CrisprGroup
 * Usage: CrisprGroup crisprgroup;
 */

CrisprGroup::CrisprGroup()
{

}

/* Destructor: ~CrisprGroup
 * Usage: delete CrisprGroup crisprgroup;
 */

CrisprGroup::~CrisprGroup()
{

}

/* Function: findPAMs
 * -------------------------------------------------------------------------------------------------------
 * Usage: A function that goes through the entire sequence that was inputted by the user.  It
 * finds all instances of "GG" which is the signal for a PAM sequence.  It then takes that location and
 * makes a new instance of gRNA in which the sequence is placed into to fill the data of the object.
 */

void CrisprGroup::findPAMs(bool dir, bool mt, vector<string> &sequences, vector<unsigned long> &compressed_seeds, vector<long> &seed_locs, vector<int> &ot_scores, vector<int> &cnts, bool &strand, pameval &PamEval, int &pam_length, int &seq_length, int &five_length, int &seed_length, string &score_file, string &on_target_data, string &endo, bool &directionality, string &pam)
{
	Scoring score(score_file, on_target_data, directionality, endo, seq_length, pam);
	vector<vector<long> > locs(sequences.size());
	vector<vector<unsigned long> > comp_seeds(sequences.size());
	vector<vector<int> > on_target_scores(sequences.size());
	int i = 0;

	if (dir == true) //if considering directionality
	{
		if (mt) //if using multi-threading, use thread pool
		{
			size_t num_threads = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), sequences.size());
			if (num_threads == 0) num_threads = 4;
			ThreadPool pool(num_threads);
			std::vector<std::future<void>> futures;
			futures.reserve(sequences.size());
			
			for (size_t k = 0; k < sequences.size(); k++)
			{
				futures.push_back(pool.enqueue([this, &locs, &cnts, &comp_seeds, &on_target_scores, &sequences, &PamEval, k, &pam_length, &seq_length, &five_length, &seed_length, &score]() {
					find_seeds_dir(locs[k], cnts[k], comp_seeds[k], on_target_scores[k], sequences[k], PamEval, k, pam_length, seq_length, five_length, seed_length, score);
				}));
			}
			for (size_t j = 0; j < futures.size(); j++)
			{
				futures[j].get();
				if ((j + 1) % 10000 == 0 || j == futures.size() - 1) {
					cout << "Processed " << (j + 1) << " / " << futures.size() << " sequences" << endl;
				}
			}
		}
		else //if using single-threading, go one at a time
		{
			for (size_t k = 0; k < sequences.size(); k++)
			{
				find_seeds_dir(locs[k], cnts[k], comp_seeds[k], on_target_scores[k], sequences[k], PamEval, k, pam_length, seq_length, five_length, seed_length, score);
				cout << "Chromosome " << k + 1 << " complete." << endl;
			}
		}
	}
	else //if not considering directionality
	{
		if (mt)
		{
			size_t num_threads = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), sequences.size());
			if (num_threads == 0) num_threads = 4;
			ThreadPool pool(num_threads);
			std::vector<std::future<void>> futures;
			futures.reserve(sequences.size());
			
			for (size_t k = 0; k < sequences.size(); k++)
			{
				futures.push_back(pool.enqueue([this, &locs, &cnts, &comp_seeds, &on_target_scores, &sequences, &PamEval, k, &seq_length, &five_length, &seed_length, &score]() {
					find_seeds(locs[k], cnts[k], comp_seeds[k], on_target_scores[k], sequences[k], PamEval, k, seq_length, five_length, seed_length, score);
				}));
			}
			for (size_t j = 0; j < futures.size(); j++)
			{
				futures[j].get();
				if ((j + 1) % 10000 == 0 || j == futures.size() - 1) {
					cout << "Processed " << (j + 1) << " / " << futures.size() << " sequences" << endl;
				}
			}
		}
		else
		{
			for (size_t k = 0; k < sequences.size(); k++)
			{
				find_seeds(locs[k], cnts[k], comp_seeds[k], on_target_scores[k], sequences[k], PamEval, k, seq_length, five_length, seed_length, score);
				cout << "Chromosome " << k + 1 << " complete." << endl;
			}
		}
	}
	
	//merge data
	for (int i = 0; i < comp_seeds.size(); i++)
	{
		for (int j = 0; j < comp_seeds[i].size(); j++)
		{
			compressed_seeds.push_back(comp_seeds[i][j]);
			seed_locs.push_back(locs[i][j]);
			ot_scores.push_back(on_target_scores[i][j]);
		}
		comp_seeds[i].clear();
		locs[i].clear();
		on_target_scores[i].clear();
		comp_seeds[i].shrink_to_fit();
		locs[i].shrink_to_fit();
		on_target_scores[i].shrink_to_fit();
	}
}

//find seeds, non-directionality
void CrisprGroup::find_seeds(vector<long> &l, int &c, vector<unsigned long> &comp_seeds, vector<int> &on_target_scores, string &seq_pointer, pameval &PamEval, int chrom, int &seq_length, int &five_length, int &seed_length, Scoring &score)
{
	//check for 5 or more N's in seq before saving it
	string seq, seed, full_seq;
	long pos = 0;
	long loc = 0;
	int cnt = 0;
	int size = seq_pointer.size();
	string curr_pam = "";
	int leftover_padding = 35 - 6 - seq_length - PamEval.pam_list[0].size();
	//loop through each valid pam
	for (int j = 0; j < PamEval.pam_list.size(); j++)
	{
		curr_pam = PamEval.pam_list[j];
		pos = 0;
		while (pos != string::npos)
		{
			if (pos >= seq_length + leftover_padding && pos < seq_pointer.size() - 35)
			{
				seq = seq_pointer.substr(pos - seq_length, seq_length);
				full_seq = seq_pointer.substr(pos - seq_length - leftover_padding, 35);
				seed = seq.substr(five_length, seed_length);
				cnt = 0;

				for (int i = 0; i < seq.size(); i++)
				{
					if (seq[i] == 'N')
					{
						cnt++;
					}
					else
					{
						cnt = 0;
					}
					if (cnt >= 5)
					{
						break;
					}
				}
				if (cnt < 5)
				{
					comp_seeds.push_back(compressSeq(seed));
					l.push_back(pos);
					on_target_scores.push_back(round(score.scoreSequence(seq, full_seq, PamEval)));
					c++;
				}
			}
			pos = seq_pointer.find(curr_pam, pos + 1);
		}
	}
	
	reverseComplement(seq_pointer);
	
	for (int j = 0; j < PamEval.pam_list.size(); j++)
	{
		pos = 0;
		curr_pam = PamEval.pam_list[j];
		while (pos != string::npos)
		{
			if (pos >= seq_length + leftover_padding && pos < seq_pointer.size() - 35)
			{
				full_seq = seq_pointer.substr(pos - seq_length - leftover_padding, 35);
				seq = seq_pointer.substr(pos - seq_length, seq_length);
				seed = seq.substr(five_length, seed_length);
				cnt = 0;

				for (int i = 0; i < seq.size(); i++)
				{
					if (seq[i] == 'N')
					{
						cnt++;
					}
					else
					{
						cnt = 0;
					}
					if (cnt >= 5)
					{
						break;
					}
				}

				if (cnt < 5)
				{
					loc = size - pos + 1;
					loc *= -1;

					comp_seeds.push_back(compressSeq(seed));
					l.push_back(loc);
					on_target_scores.push_back(round(score.scoreSequence(seq, full_seq, PamEval)));
					c++;
				}
			}
			pos = seq_pointer.find(curr_pam, pos + 1);
		}
	}
	reverseComplement(seq_pointer);
}

//find seeds - directionality
void CrisprGroup::find_seeds_dir(vector<long> &l, int &c, vector<unsigned long> &comp_seeds, vector<int> &on_target_scores, string &seq_pointer, pameval &PamEval, int chrom, int &pam_length, int &seq_length, int &five_length, int &seed_length, Scoring &score)
{
	//check for 5 or more N's in seq before saving it
	string seq, seed, full_seq;
	long pos = 0;
	int cnt = 0;
	long loc = 0;
	string curr_pam = "";
	int size = seq_pointer.size();
	int leftover_padding = 35 - 6 - seq_length - PamEval.pam_list[0].size();

	for (int j = 0; j < PamEval.pam_list.size(); j++)
	{
		curr_pam = PamEval.pam_list[j];
		pos = 0;
		while (pos != string::npos)
		{
			if (pos >= seq_length + leftover_padding && pos < seq_pointer.size() - 35)
			{
				seq = seq_pointer.substr(pos + pam_length, seq_length);
				full_seq = seq_pointer.substr(pos - 6, 35);
				seed = seq.substr(five_length, seed_length);
				cnt = 0;
				for (int i = 0; i < seq.size(); i++)
				{
					if (seq[i] == 'N')
					{
						cnt++;
					}
					else
					{
						cnt = 0;
					}
					if (cnt >= 5)
					{
						break;
					}
				}
				if (cnt < 5)
				{
					comp_seeds.push_back(compressSeq(seed));
					l.push_back(pos + 1);
					on_target_scores.push_back(round(score.scoreSequence(seq, full_seq, PamEval)));
					c++;
				}
			}
				
			pos = seq_pointer.find(curr_pam, pos + 1);
		}
	}

	reverseComplement(seq_pointer);
	for (int j = 0; j < PamEval.pam_list.size(); j++)
	{
		curr_pam = PamEval.pam_list[j];
		pos = 0;
		while (pos != string::npos)
		{
			if (pos >= seq_length + leftover_padding && pos < seq_pointer.size() - 35)
			{
				seq = seq_pointer.substr(pos + pam_length, seq_length);
				full_seq = seq_pointer.substr(pos - 6, 35);
				seed = seq.substr(five_length, seed_length);
				cnt = 0;
				for (int i = 0; i < seq.size(); i++)
				{
					if (seq[i] == 'N')
					{
						cnt++;
					}
					else
					{
						cnt = 0;
					}
					if (cnt >= 5)
					{
						break;
					}
				}

				if (cnt < 5)
				{
					loc = size - pos;
					loc *= -1;
					comp_seeds.push_back(compressSeq(seed));
					l.push_back(loc);
					on_target_scores.push_back(round(score.scoreSequence(seq, full_seq, PamEval)));
					c++;
				}
			}
			pos = seq_pointer.find(curr_pam, pos + 1);
		}
	}
	reverseComplement(seq_pointer);
}

void CrisprGroup::process_targets(vector<int> &uniques, vector<int> &repeats, vector<unsigned long> &compressed_seeds, vector<long> &seed_locs)
{
	//sort seeds
	int i = 0;
	vector<int> indices(compressed_seeds.size());
	
	//sort a vector of indices that would follow a sort based on compressed seed vector
	iota(indices.begin(), indices.end(), 0);
	sort(indices.begin(), indices.end(),
		[&compressed_seeds](int A, int B) -> bool {
		return compressed_seeds[A] < compressed_seeds[B];
	});

	//now we know that if compressed seeds are sorted, a repeated seed will be next to its repeats in the vector
	//loop through compressed seeds, find out if that index leads to a unique or a repeated seed	
	if (compressed_seeds[indices[0]] != compressed_seeds[indices[1]])
	{
		uniques.push_back(indices[0]);
	}
	else
	{
		repeats.push_back(indices[0]);
	}
	for (i = 1; i < indices.size() - 2; i++)
	{
		if (compressed_seeds[indices[i - 1]] != compressed_seeds[indices[i]] && compressed_seeds[indices[i]] != compressed_seeds[indices[i + 1]])
		{
			uniques.push_back(indices[i]);
		}
		else
		{
			repeats.push_back(indices[i]);
		}
	}
	//handle last
	if (compressed_seeds[indices[i - 1]] != compressed_seeds[indices[i]])
	{
		uniques.push_back(indices[i]);
	}
	else
	{
		repeats.push_back(indices[i]);
	}

	sort(uniques.begin(), uniques.end());

	//cleanup - clear out indices vector
	indices.clear();
	indices.shrink_to_fit();
}

/* Function: reverseComplement
 * --------------------------------------------------------------------------------------------------------
 * Usage: a sequence in the form of a string is passed in by reference and the function returns the reverse
 * complement of the passed in sequence, inserting X's if there are any nucleotide discrepancies in the
 * original sequence.
*/

void CrisprGroup::reverseComplement(string &str)
{
	string rc = "";
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
	str = rc;
}

unsigned long CrisprGroup::compressSeq(string &s) {
	unsigned long compseq = 0;
	for (int i = 0; i < s.size(); i++)
	{
		compseq += convertCharBase4(s[i])*pow(4, i); //multiplying by power-4 converts to base10
	}
	return compseq; //base10 version of sequence string
}

int CrisprGroup::convertCharBase4(char &c) {
	switch (c) {
	case 'A': return 0;
	case 'T': return 1;
	case 'C': return 2;
	case 'G': return 3;
	default: return 0;
	}
}