#include <iostream>
#include <stdio.h>
#include <chrono>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <typeinfo>
#include <getopt.h>
#include "../headers/index_color.h"
#include "../headers/utils.h"
#include "../headers/MinimizerLister.h"
#include "../TurboPFor-Integer-Compression/include/ic.h"


using namespace std;
using namespace chrono;




string query_file(""), output_prefix("hits_"), binary_prefix("binary_index");
bool fof(false), fasta(false);
double threshold(0.4);
uint16_t num_thread(1);



void PrintHelp()
{
	cout <<
			"\n******************* K2R **************************************\n"


			"\n INDEX QUERY \n"
            "-s OR -f                 :     Sequence file (FASTA) if a unique sequence is queried (-s), file of file if several sequences are queried (-f) \n"
			"-o                       :     Write output reads in fasta file (default : query_output) \n"
            "-b                       :     Read index from binary files (default: binary_index) \n"

			"\n TWEAK PARAMETERS\n"
            "-r                       :     Rate of minimizer found in the read to keep it in results (between 0 and 1, default: " << intToString(threshold) << ")\n";

	exit(1);
}



void ProcessArgs(int argc, char** argv)
{
	const char* const short_opts = "s:f:o:b:t:r:";
	const option long_opts[] =
	{
	};
	while (true)
	{
		const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

		if (-1 == opt)
			break;

		switch (opt)
		{
            case 's':
				query_file=optarg;
                fasta=true;
				break;
            case 'f':
				query_file=optarg;
                fof=true;
				break;
            case 'o':
				output_prefix=optarg;
				break;
            case 'b':
				binary_prefix=optarg;
				break;
            case 't':
				num_thread=stoi(optarg);
				break;
            case 'r':
				threshold=stod(optarg);
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
    
    ProcessArgs(argc, argv);
    if (argc < 2){
        PrintHelp();
    }else{
        string bin_file_mmer = binary_prefix + "_mmer.bin";
        string bin_file_color = binary_prefix + "_color.bin";
        Index_color index_color = Index_color(bin_file_mmer, bin_file_color);
        index_color.read_stream=new ifstream(index_color.filename, ios::in);

        int opt;

        if(fof){
            auto start_querying_seq = high_resolution_clock::now();
            index_color.query_fof(query_file, output_prefix, threshold, num_thread);
            auto end_querying_seq = high_resolution_clock::now();
            auto querying_seq = duration_cast<nanoseconds>(end_querying_seq - start_querying_seq);
            cout << "Querying the file takes " << (float)querying_seq.count() << " ns." << endl;
        }
        else if(fasta){
            auto start_querying_seq = high_resolution_clock::now();
            index_color.query_fasta(query_file, output_prefix, threshold, num_thread);
            auto end_querying_seq = high_resolution_clock::now();
            auto querying_seq = duration_cast<nanoseconds>(end_querying_seq - start_querying_seq);
            cout << "Querying the sequence takes " << (float)querying_seq.count() << " ns." << endl;
        }
        else{
            PrintHelp();
        }
    }
    return 0;
}
