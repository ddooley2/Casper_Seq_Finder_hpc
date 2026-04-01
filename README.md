# **CASPER SeqFinder Algorithm**

This algorithm searches through an input genomic sequence file (FASTA/FNA) and finds all unique and repeated guide-RNA (gRNA)
sequences for a given Cas endonuclease. It outputs two files:

1. CSPR file containing a list of all *unique* gRNAs
2. SQL DB file containing all *repeated* gRNAs

## Download And Compile Sqlite3
### Linux: 
1. Open terminal
2. Run the following command: `sudo apt-get install sqlite3`

### Mac and Linux (if you did not use apt-get):
1. Download sqlite3 source code
2. Open terminal
3. CD to sqlite3 source code folder
3. Build using the following command: `gcc -c sqlite3.c`
	* The build command should generate a .o file that will be used with compiling SeqFinder.
4. Copy the .o file to the SeqFinder source code folder

### Windows (Visual Studio 2017):
1. Download sqlite3 source code
2. Open "Developer Command Prompt for VS 2017"
3. CD to sqlite3 source code folder
4. Build the library files by running the following command: `lib /DEF:sqlite3.def /OUT:sqlite3.lib /MACHINE:x64`

## Download and Compile Seq-Finder
### Linux (if you used apt-get install sqlite3):
1. Download SeqFinder source code for Linux (in Repository)
2. Open terminal
3. CD to Seq-Finder source code folder
3. Run Command to compile SeqFinder: `g++ -std=c++11 -O2 *.cpp -pthread -l sqlite3 -o SeqFinder`
   * The `-O2` flag enables optimizations for better performance

### Mac and Linux (if you manually built sqlite3 .o file, make sure sqlite3 .o file is in same folder as SeqFinder source code):
1. Download SeqFinder source code for Mac or Linux (in Repository)
2. Open terminal
3. CD to SeqFinder source code folder
3. Run command to compile SeqFinder: `g++ -std=c++11 -O2 *.cpp -pthread sqlite3.o -o SeqFinder`
   * The `-O2` flag enables optimizations for better performance

### Windows (Visual Studio 2017):
1. Download SeqFinder source for Windows (in Repository)
2. Open Visual Studio 2017
3. Create a New Project
4. Import the SeqFinder source code files
5. Go to Project->Properties and do the following:
	* Set Configuration to "All Configurations"
	* Set Platform to "All Platforms"
	* In Configuration Properties->VC++ Directories add `C:\Path\To\sqlite3;` to "Include Directories"
	* In Configuration Properties->VC++ Directories add `C:\Path\To\sqlite;` to "Include Libraries"
	* In C/C++->General add `C:\Path\To\sqlite3;` to "Additional Include Directories" and "Additional #using Directories"
	* In Linker->General add `C:\Path\To\sqlite;` to "Additional Library Directories"
	* In Linker->Input add `C:\Path\To\sqlite\sqlite3.lib;` in "Additonaly Dependencies
	* Make sure when you add these paths that there are ';' seperating all paths/object in the line.
6. For debugging, make sure Debug and x64 are selected before running
	* If debugging, make sure you set the debugging command arguments. See "How to run SeqFinder" below.
7. For compiling, make sure Release and x64 are selected before running
	

## How to run SeqFinder
* CD to the directory containing the SeqFinder executable
	* The command line arguments for SeqFinder are as follows: `endonuclease pam multi_threading_bool 5_prime_pam_bool generate_repeats_file_bool
  5_prime_length seed_length 3_prime_length org_code output_folder_path CASPERinfo_file_path input_file_path org_name notes on_target_data [num_threads]`

* Example command: `./SeqFinder spCas9 NGG TRUE FALSE TRUE 4 16 0 org_code output_folder CASPERinfo input.fasta org_name notes DATA:CRISPRSCAN`
* Example with thread count: `./SeqFinder spCas9 NGG TRUE FALSE TRUE 4 16 0 org_code output_folder CASPERinfo input.fasta org_name notes DATA:CRISPRSCAN 16`

### Memory and Performance Recommendations
* **Small genomes** (<100 scaffolds, <500 Mb): Default settings work well
* **Medium genomes** (100-10,000 scaffolds, 500Mb-2Gb): Use 8-16 threads for optimal performance
* **Large fragmented genomes** (>10,000 scaffolds, >2Gb): 
  * Use 16-32 threads depending on available CPU cores
  * Ensure sufficient RAM: approximately 3-4x the genome size
  * Example: For a 1.9Gb genome with 583,000 scaffolds, allocate at least 8GB RAM
