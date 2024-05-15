#include <iostream>
#include <stdio.h>
#include <chrono>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <typeinfo>
#include <getopt.h>
// #include "../BitMagic/src/bm.h"
// #include "../BitMagic/src/bmundef.h"
// #include "../BitMagic/src/bmsparsevec.h"
#include "../headers/index_color.h"
#include "../headers/utils.h"
#include "../headers/MinimizerLister.h"
#include "../TurboPFor-Integer-Compression/include/ic.h"
#include <map>

using namespace std;
using namespace chrono;




string query_file(""), output_prefix("query_output"), binary_prefix("binary_index");
bool sequence(false), fasta(false);



void PrintHelp()
{
	cout <<
			"\n******************* K2R **************************************\n"


			"\n INDEX STAT \n"
			"--output-prefix          :     Write output reads in fasta file (default : query_output) \n"
            "--binary-prefix          :     Read index from binary files (default: binary_index) \n";

	exit(1);
}




void ProcessArgs(int argc, char** argv)
{
	const char* const short_opts = "";
	const option long_opts[] =
	{
		{"output-prefix", required_argument, nullptr, 'p'},
        {"binary-prefix", required_argument, nullptr, 'b'},
		{nullptr, no_argument, nullptr, 0}
	};
	while (true)
	{
		const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);

		if (-1 == opt)
			break;

		switch (opt)
		{
            case 'p':
				output_prefix=optarg;
				break;
            case 'b':
				binary_prefix=optarg;
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

        //
        ifstream file_mmer(bin_file_mmer, ios::binary);
        file_mmer.seekg(0, ios::end);
        uint64_t file_mmer_size_bytes = file_mmer.tellg();
        ifstream file_color(bin_file_color, ios::binary);
        file_color.seekg(0, ios::end);
        uint64_t file_color_size_bytes = file_color.tellg();


        uint64_t nb_mmer = index_color.mmermap.size();
        uint64_t nb_color = 0;
        uint64_t norm_color_compress = 0;
        uint64_t norm_color_decompress = 0;
        uint64_t sum_diff_color = 0;
        uint64_t sum_diff_color_log = 0;
        uint32_t actual_size;
        std::map<uint64_t, uint64_t> ab{};
        uint64_t ab_le_1 = 0;
        uint64_t ab_le_2 = 0;
        uint64_t ab_le_4 = 0;
        uint64_t ab_le_8 = 0;
        uint64_t ab_le_16 = 0;
        uint64_t ab_le_32 = 0;
        uint64_t ab_le_64 = 0;
        uint64_t last = 0;
        for(uint i(0);i<1024;i++){
            color_map::iterator it = index_color.colormap[i].begin();
            nb_color += index_color.colormap[i].size();
            while (it != index_color.colormap[i].end()){
                norm_color_compress += it->second.compressed_array_size;
                ab[actual_size]++;
                if (actual_size <= 1) {
                    ab_le_1++;
                }
                if (actual_size <= 2) {
                    ab_le_2++;
                }
                if (actual_size <= 4) {
                    ab_le_4++;
                }
                if (actual_size <= 8) {
                    ab_le_8++;
                }
                if (actual_size <= 16) {
                    ab_le_16++;
                }
                if (actual_size <= 32) {
                    ab_le_32++;
                }
                if (actual_size <= 64) {
                    ab_le_64++;
                }
                vector<uint32_t> v = index_color.decompress_color(it->second);
                last = 0;
                for(uint32_t u : v) {
                    sum_diff_color += u - last;
                    sum_diff_color_log += log2(u - last + 1);
                    last = u;
                }

                it++;
            }
        }

        ofstream out_file;
        out_file.open (output_prefix+".csv");
        out_file << "filename (params)" << ",";
        out_file << "binary_prefix (params)" << ",";
        out_file << "k (params)" << ",";
        out_file << "m (params)" << ",";
        out_file << "size mmer file" << ",";
        out_file << "size color file" << ",";
        out_file << "nb mmer" << ",";
        out_file << "nb color" << ",";
        out_file << "sum size color compress" << ",";
        out_file << "sum size color decompress" << ",";
        out_file << "sum diff color" << ",";
        out_file << "sum diff color log" << ",";
        out_file << "ab lte 1" << ",";
        out_file << "ab lte 2" << ",";
        out_file << "ab lte 4" << ",";
        out_file << "ab lte 8" << ",";
        out_file << "ab lte 16" << ",";
        out_file << "ab lte 32" << ",";
        out_file << "ab lte 64" << endl;

        out_file << index_color.filename << ",";
        out_file << index_color.binary_prefix << ",";
        out_file << index_color.k << ",";
        out_file << index_color.m << ",";
        out_file << file_mmer_size_bytes << ",";
        out_file << file_color_size_bytes << ",";
        out_file << nb_mmer << ",";
        out_file << nb_color << ",";
        out_file << norm_color_compress << ",";
        out_file << norm_color_decompress << ",";
        out_file << sum_diff_color << ",";
        out_file << sum_diff_color_log << ",";
        out_file << ab_le_1 << ",";
        out_file << ab_le_2 << ",";
        out_file << ab_le_4 << ",";
        out_file << ab_le_8 << ",";
        out_file << ab_le_16 << ",";
        out_file << ab_le_32 << ",";
        out_file << ab_le_64 << endl;


        out_file.close();


        // for (const auto& [key, value] : ab) {
        //     cout << '[' << key << "] = " << value << "; ";
        // }
    }
    return 0;
}
