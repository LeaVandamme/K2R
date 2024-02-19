#include <iostream>
#include <stdio.h>
#include <chrono>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <typeinfo>
#include <getopt.h>
#include "../include/fastDelta.h"
#include "../headers/index_color.h"
#include "../headers/utils.h"
#include "../headers/MinimizerLister.h"
#include "../TurboPFor-Integer-Compression/include/ic.h"


using namespace std;
using namespace chrono;




string query_file(""), output_prefix("hits_"), binary_prefix("binary_index");
bool fof(false), fasta(false);
double threshold(1);
uint16_t num_thread(1);



void PrintHelp()
{
	cout <<
			"\n******************* K2R **************************************\n"


			"\n INDEX QUERY \n"
            "--query-file             :     Sequence(s) to query (FASTA allowed)\n"
			"--output-prefix          :     Write output reads in fasta file (default : query_output) \n"
            "--binary-prefix          :     Read index from binary files (default: binary_index) \n"
            "--thread                 :     Number of threads used (default: 1)\n"

			"\n TWEAK PARAMETERS\n"
            "-s                       :     When querying a unique sequence \n"
            "-t                       :     Rate of minimizer found in the read to keep it in results (between 0 and 1, default: " << intToString(threshold) << ")\n"
            "-f                       :     When querying several sequences \n";

	exit(1);
}



void ProcessArgs(int argc, char** argv)
{
	const char* const short_opts = "sft:";
	const option long_opts[] =
	{
		{"query-file", required_argument, nullptr, 'q'},
		{"output-prefix", required_argument, nullptr, 'p'},
        {"binary-prefix", required_argument, nullptr, 'b'},
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
            case 'q':
				query_file=optarg;
				break;
            case 'p':
				output_prefix=optarg;
				break;
            case 'b':
				binary_prefix=optarg;
				break;
            case 'th':
				num_thread=stoi(optarg);
				break;
            case 's':
				fasta=true;
				break;
            case 't':
				threshold=stod(optarg);
				break;
			case 'f':
				fof=true;
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
        string binary_file_mmermap, binary_file_colormap;


        binary_file_colormap = binary_prefix + "_color.bin";
        binary_file_mmermap = binary_prefix + "_mmer.bin";

        if(fof){
            auto start_querying_seq = high_resolution_clock::now();
            index_color.query_fof(query_file, output_prefix, threshold, num_thread);
            auto end_querying_seq = high_resolution_clock::now();
            auto querying_seq = duration_cast<nanoseconds>(end_querying_seq - start_querying_seq);
            cout << index_color.minimizer_match << endl;
            cout << "Querying the sequence takes " << (float)querying_seq.count() << " ns." << endl;
        }
        else if(fasta){
            string output_path = output_prefix+query_file;
            index_color.deserialize_mmermap(binary_file_mmermap);
            index_color.deserialize_colormap(binary_file_colormap);
            auto start_querying_seq = high_resolution_clock::now();
            index_color.query_fasta(query_file, output_path, threshold, num_thread);
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
