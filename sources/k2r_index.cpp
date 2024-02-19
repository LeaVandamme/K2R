#include <iostream>
#include <stdio.h>
#include <chrono>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <typeinfo>
#include <getopt.h>
#include <omp.h>
#include "../BitMagic/src/bm.h"
#include "../BitMagic/src/bmundef.h"
#include "../BitMagic/src/bmsparsevec.h"
#include "../include/fastDelta.h"
#include "../headers/index_color.h"
#include "../headers/utils.h"
#include "../headers/MinimizerLister.h"
#include "../TurboPFor-Integer-Compression/include/ic.h"

using namespace std;
using namespace chrono;

uint16_t k(31);
uint16_t m(15);
uint16_t counting_bf_size(32);
uint16_t min_occ(2);
uint16_t num_thread(1);
string read_file(""), binary_prefix("binary_index");
bool keep_unique(false);
bool homocompression(false);

void PrintHelp()
{
	cout <<
			"\n******************* K2R **************************************\n"


			"\n INDEX CONSTRUCTION\n"
            "--read-file              :     Build index from file (FASTA allowed)\n"
            "--keep-unique            :     Keep all minimizers in index, included uniques \n"
			"--binary-prefix          :     Write index in binary files (default : binary_index) \n"
            "--homocompression        :     Homocompression of reads \n"
            "--thread                 :     Number of threads used (default: 1)\n"

			"\n TWEAK PARAMETERS\n"
			"-k                       :     k-mer size (default: " << intToString(k) << ")\n"
            "-m                       :     Minimizer size (default: " << intToString(m) << ")\n"
			"-b                       :     Counting bloom filter size (log(2), default " << intToString(counting_bf_size) << ")\n"
            "-o                       :     Minimum occurence of minimizers (default: << " << min_occ << ")\n";

	exit(1);
}



void ProcessArgs(int argc, char** argv)
{
	const char* const short_opts = "k:m:t:b:o:";
	const option long_opts[] =
	{
		{"read-file", required_argument, nullptr, 'r'},
        {"keep-unique", no_argument, nullptr, 'u'},
		{"binary-prefix", required_argument, nullptr, 'p'},
        {"homocompression", no_argument, nullptr, 'h'},
        {"thread", required_argument, nullptr, 'th'},
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
            case 'u':
				keep_unique=true;
				break;
			case 'p':
				binary_prefix = optarg;
				break;
			case 'h':
				homocompression=true;
				break;
            case 'th':
				num_thread=stoi(optarg);
				break;
            case 'k':
				k=stoi(optarg);
				break;
            case 'm':
				m=stoi(optarg);
				break;
			case 'b':
				counting_bf_size=stoi(optarg);
				break;
			case 'o':
				min_occ=stoi(optarg);
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
        ProcessArgs(argc, argv);
        if (argc < 2){
            PrintHelp();
        }
        Index_color index_color = Index_color(read_file, k, m, counting_bf_size, binary_prefix);

        string file_out, prefix_binary_file, binary_file_mmermap, binary_file_colormap;

        if(keep_unique){
           cout<<"TODO"<<endl;
        }else{
            auto start_indexing = high_resolution_clock::now();
            index_color.create_index_mmer_no_unique(read_file, k, m, min_occ, counting_bf_size, homocompression, num_thread);
            auto end_indexing = high_resolution_clock::now();
            auto indexing = duration_cast<seconds>(end_indexing - start_indexing);

            cout << "\n";
            cout << "Indexing takes " << indexing.count() << " seconds." << endl;
            cout << "M-mer indexed : " << intToString(index_color.mmermap.size()) << " | " << endl;
            uint64_t memory_used_index = getMemorySelfMaxUsed();
            cout << "Resource usage (indexing) : " << intToString(memory_used_index) << " Ko" << endl;

            binary_file_colormap = binary_prefix + "_color.bin";
            binary_file_mmermap = binary_prefix + "_mmer.bin";
            auto start_serializing = high_resolution_clock::now();
            index_color.serialize_mmermap(binary_file_mmermap);
            index_color.serialize_colormap(binary_file_colormap);
            auto end_serializing = high_resolution_clock::now();
            auto serializing = duration_cast<seconds>(end_serializing - start_serializing);
            cout << "Serialization takes : " << serializing.count() << " seconds." << endl;

        }
        return 0;
    }
        

