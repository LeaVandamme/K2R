#include <iostream>
#include <stdio.h>
#include <chrono>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <typeinfo>
#include <getopt.h>
#include <omp.h>
#include "../headers/index_color.h"
#include "../headers/utils.h"
#include "../headers/MinimizerLister.h"
#include "../TurboPFor-Integer-Compression/include/ic.h"

using namespace std;
using namespace chrono;

uint16_t k(31);
uint16_t m(15);
uint16_t counting_bf_size(32);
uint16_t min_ab(2);
uint16_t max_ab(1000);
uint16_t num_thread(1);
string read_file(""), binary_prefix("binary_index");
bool max_ab_on(false);
bool homocompression(false);

void PrintHelp()
{
	cout <<
			"\n******************* K2R **************************************\n"


			"\n INDEX CONSTRUCTION\n"
            "-r                       :     Build index from file (FASTA allowed)\n"
			"-b                       :     Write index in binary files (default : binary_index) \n"
            "-t                       :     Number of threads used (default: 1, max: 5)\n"

			"\n TWEAK PARAMETERS\n"
			"-k                       :     k-mer size (default: " << intToString(k) << ")\n"
            "-m                       :     Minimizer size (default: " << intToString(m) << ")\n"
			"-s                       :     Counting bloom filter size (log(2), default " << intToString(counting_bf_size) << ")\n"
			"-h                       :     Homocompression of reads \n"
			"--min-ab                 :     Minimizers minimum abundance (default: << " << min_ab << ")\n"
			"--max-ab                 :     Minimizers maximum abundance (default: << " << max_ab << ", must be < 2^16)\n";


	exit(1);
}



void ProcessArgs(int argc, char** argv)
{
	const char* const short_opts = "r:b:t:hk:m:s:";
	const option long_opts[] =
	{
		{"min-ab", required_argument, nullptr, 'min'},
		{"max-ab", required_argument, nullptr, 'max'},
		{nullptr, no_argument, nullptr, 0}
	};
	while (true)
	{
		const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

		if (-1 == opt)
			break;

		switch (opt)
		{
            case 'r':
				read_file=optarg;
				break;
			case 'b':
				binary_prefix = optarg;
				break;
			case 'h':
				homocompression=true;
				break;
            case 't':
				num_thread=stoi(optarg);
				if(num_thread > 5){
					num_thread = 5;
				}
				break;
            case 'k':
				k=stoi(optarg);
				break;
            case 'm':
				m=stoi(optarg);
				break;
			case 's':
				counting_bf_size=stoi(optarg);
				break;
			case 'min':
				min_ab=stoi(optarg);
				break;
			case 'max':
				max_ab_on = true;
				max_ab=stoi(optarg);
				break;
			case '?':
				PrintHelp();
				break;
			default:
				PrintHelp();
				break;
		}
	}
}




int main(int argc, char *argv[]){

		putenv("OMP_STACKSIZE=128M");

		k=(31);
		m=(15);
		counting_bf_size=(32);
		min_ab=(2);
		max_ab=(1000);
		num_thread=(1);
		read_file=("");
		binary_prefix=("binary_index");
		max_ab_on=(false);
		homocompression=(false);



        ProcessArgs(argc, argv);
        if (argc < 2){
            PrintHelp();
        }

        Index_color index_color = Index_color(read_file, k, m, counting_bf_size, binary_prefix);

        string file_out, prefix_binary_file, binary_file_mmermap, binary_file_colormap;
		auto start_indexing = high_resolution_clock::now();
		cout << "INDEX CONSTRUCTION" << endl;
		cout << "=========================================================" << endl << endl;
		cout << "Start indexing the file : " + read_file << endl << endl;
		index_color.create_index_mmer_no_unique(read_file, k, m, min_ab, max_ab, counting_bf_size, homocompression, num_thread);
		auto end_indexing = high_resolution_clock::now();
		auto indexing = duration_cast<seconds>(end_indexing - start_indexing);

		cout << "Indexing takes " << indexing.count() << " seconds." << endl;
		uint64_t memory_used_index = getMemorySelfMaxUsed();
		cout << "Resource usage : " << intToString(memory_used_index) << " Ko" << endl;
		cout << endl << "SERIALIZATION" << endl;
		cout << "=========================================================" << endl;

		binary_file_colormap = binary_prefix + "_color.bin";
		binary_file_mmermap = binary_prefix + "_mmer.bin";
		auto start_serializing = high_resolution_clock::now();
		cout << endl << "Beginning of serialization in files : " + binary_file_mmermap + " & " + binary_file_colormap << endl;
		index_color.serialize_mmermap(binary_file_mmermap);
		index_color.serialize_colormap(binary_file_colormap);
		auto end_serializing = high_resolution_clock::now();
		auto serializing = duration_cast<seconds>(end_serializing - start_serializing);
		cout << "Serialization takes : " << serializing.count() << " seconds." << endl;

        return 0;
    }
