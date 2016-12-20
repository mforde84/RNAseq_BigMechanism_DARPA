#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <thread>
#include <sstream>
#include <fstream>
#include <math.h> 
#include <algorithm> 
#include "api/BamReader.h"

using namespace std;
using namespace BamTools;

struct bam_headers{
	
	string chr;
	int start;
	unsigned long lib_size;

};

struct gtf_headers{

	string chr;
	string source;
	string feature;
	string score;
	string strand;
	string frame;
	string annotation;
	int start;
	int end;
	short int bin;
	vector <int> counts;
	vector <float> rpkm;
	
};

vector <string> find_index(vector <vector <bam_headers> > &bams){
	
	//define vector for bam_index to chromosome
	
	vector <string> compute_holder;
	
	for (unsigned int bam_idx = 0; bam_idx < bams.size();bam_idx++){
		if (bams[bam_idx].size() == 0){
			compute_holder.push_back("EMPTY");
		}else{
			compute_holder.push_back(bams[bam_idx][0].chr);
		}
	}
	return compute_holder;
	
}

vector <vector <gtf_headers> > load_gtf(char* filename, int nthreads){
	
	vector<gtf_headers> push_matrix;
	gtf_headers holder;
	ifstream gtf_file(filename);
	string line;
	
	//iterate through gtf file, and load to memory
	//cout << "Loading GTF to memory" << "\n";
	if (gtf_file.is_open()){
		int sub_count = 0;
		string transfer_hold[8];
		while(getline(gtf_file,line)){
			istringstream iss(line);
			string token;
			while(getline(iss,token,'\t')){
				if (sub_count == 8){
					holder.chr = transfer_hold[0];
					holder.source = transfer_hold[1];
					holder.feature = transfer_hold[2];
					holder.start = atoi(transfer_hold[3].c_str());
					holder.end = atoi(transfer_hold[4].c_str());
					holder.score = transfer_hold[5];
					holder.strand = transfer_hold[6];
					holder.frame = transfer_hold[7];
					holder.annotation = token;
					push_matrix.push_back(holder);
					sub_count = 0;
				} else {
					transfer_hold[sub_count] = token;
					++sub_count;
				}
			}
		}
		//breakup gtf records to bins by num threads
		unsigned int divider_count = push_matrix.size() / nthreads;
		unsigned int remainder_count = push_matrix.size() % nthreads;
		vector <vector <gtf_headers> > thread_push;
		vector <gtf_headers> init_push;
		for (int i = 1; i < nthreads + 1; i++){
			thread_push.push_back(init_push);
			if (i < nthreads){
				for (unsigned int x = 1; x < divider_count + 1; x++){
					thread_push.at(i-1).push_back(push_matrix[(((i*divider_count)-divider_count)+x)-1]);
				}
			}else{
				for (unsigned int x = 1; x < divider_count + remainder_count + 1; x++){
					thread_push.at(i-1).push_back(push_matrix[(((i*divider_count)-divider_count)+x)-1]);
				}
			}
		}
		gtf_file.close();
		return(thread_push);
	}else{
		//cout << "GTF unsuccessfully loaded to memory. Check path to file, and annotation format. Exiting" << "\n";
		exit(-1);
	}
	
}

vector <vector <bam_headers> > load_bam(char* filename){
	
	unsigned long lib_count = 0;
	vector <vector <bam_headers> > push_matrix;
	vector <bam_headers> iter_chr;
	int iter_refid = -1;
	bam_headers bam_holder;
	BamReader reader;
	BamAlignment al;
	const vector<RefData>& references = reader.GetReferenceData();
	
	//load bam file to memory, bin by chromosome
	//bam file must be sorted by chr. otherwise the lookup will segfault	
	//cout << "Loading " << filename << " to memory" << "\n";
	if (reader.Open(filename)) {	
		while (reader.GetNextAlignmentCore(al)) {
			if (al.IsMapped()){
				lib_count++;
				if(al.RefID != iter_refid){
					iter_refid++;
					push_matrix.push_back(iter_chr);
				}else{
					bam_holder.chr = references[al.RefID].RefName;
					bam_holder.start = al.Position;
					push_matrix.at(iter_refid).push_back(bam_holder);
				}
			}
		}
		push_matrix[0][0].lib_size = lib_count;
		reader.Close();
		return(push_matrix);
	}else{
		//cout << "Could not open input BAM file. Exiting." << endl;
		exit(-1);
	}
	
}

short int find_bin(const string & gtf_chr, const vector <string> mapping){
	
	//f which chr. bin the gtf line is associated with 
	int bin_compare = -1;
	for (unsigned int i = 0; i < mapping.size(); i++){
		if(gtf_chr == mapping[i]){ 
			bin_compare = i;
		}
	}
	return(bin_compare);

}
	
vector<vector <int> > parse_chromosome ( vector <vector <bam_headers> > &bam_matrix){

	//cout << "Parsing chromosomes\n";
	vector <vector <int> > chr_matrix;
	vector <int> place_holder;
	
	//moves the bam chromosome bins to searchable 0 indexed vector
	for (unsigned int i = 0; i < bam_matrix.size(); i++){
		chr_matrix.push_back(place_holder);
		for (unsigned int x = 0; x < bam_matrix[i].size();x++){
				//cout << bam_matrix[i][x].start << "\n";
				chr_matrix.at(i).push_back(bam_matrix[i][x].start);
		}
	}
	//cout << "works" << "\n";
	return(chr_matrix);

}

int find_frame(vector <int> &chr_matrix, int start_gtf, int stop_gtf){
	
	//binary search for range of bam alignments satisfing gtf annotation range
	//dcount is equal to size of subrange			
	int counts;
	vector <int>::iterator low,up;
	low =std::lower_bound (chr_matrix.begin(), chr_matrix.end(), start_gtf); 
	up = std::upper_bound (chr_matrix.begin(), chr_matrix.end(), stop_gtf); 
	counts = (up - chr_matrix.begin()) - (low - chr_matrix.begin());
	return(counts);
			
} 

void write_counts(vector <vector <gtf_headers> > &gtf_matrix, char *filenames[], int args){

	//output headers
	cout << "chr" << "\t" << "source" << "\t" << "feature" << "\t" << "start" << "\t" << "end" << "\t" << "score" << "\t" << "strand" << "\t" << "frame" << "\t" << "annotation";
	for (int z = 3; z < args; z++){
		if (z < args - 1){
			cout << "\t" << filenames[z] << "\t" << filenames[z] << "_rpkm";
		}else{
			cout << "\t" << filenames[z] << "\t" << filenames[z] << "_rpkm" << "\n";
		}
	}
	//write counts
	for (unsigned int i = 0; i < gtf_matrix.size(); i++){		
		for (unsigned int x = 0; x < gtf_matrix[i].size();x++){
			cout << gtf_matrix[i][x].chr << "\t" << gtf_matrix[i][x].source << "\t" << gtf_matrix[i][x].feature << "\t" << gtf_matrix[i][x].start << "\t" << gtf_matrix[i][x].end << "\t" << gtf_matrix[i][x].score << "\t" << gtf_matrix[i][x].strand << "\t" << gtf_matrix[i][x].frame << "\t" << gtf_matrix[i][x].annotation << "\t";
			for (unsigned int y = 0; y < gtf_matrix[i][x].counts.size(); y++){
				if (y == gtf_matrix[i][x].counts.size() - 1){
					cout << gtf_matrix[i][x].counts[y] << "\t" << gtf_matrix[i][x].rpkm[y] << "\n";
				}else{
					cout << gtf_matrix[i][x].counts[y] << "\t" << gtf_matrix[i][x].rpkm[y] << "\t";
				}
			} 
		}
	}
		
}

void runtime_threads(vector <gtf_headers> &gtf_bin, vector <vector <int> > &chr_matrix, vector <string> &index_map, unsigned long lib_size){
	
	//thread process
	for (unsigned int i = 0; i < gtf_bin.size();i++){//gtf_bin.size()
		int holder_count = 0;
		gtf_bin[i].bin = find_bin(gtf_bin[i].chr,index_map);
		if(gtf_bin[i].bin != -1){
			holder_count = find_frame(std::ref(chr_matrix[gtf_bin[i].bin]), gtf_bin[i].start, gtf_bin[i].end);
			gtf_bin[i].counts.push_back(holder_count);
			gtf_bin[i].rpkm.push_back(((float)holder_count / ((float)lib_size / 1000000)) / (((float)gtf_bin[i].end - (float)gtf_bin[i].start) / 1000));
		}else{
			gtf_bin[i].counts.push_back(0);
			gtf_bin[i].rpkm.push_back(0);
		}
	}
	
}

void define_threads(vector <vector <int> > &chr_matrix, int nthreads, vector <vector <gtf_headers> > &gtf_matrix, vector <string> index_map, unsigned long lib_size){
	
	//instantiates then joins threads
	vector <thread> thread_holder;
	for (int i = 0; i < nthreads; i++){
		thread_holder.push_back(thread(runtime_threads,std::ref(gtf_matrix[i]),std::ref(chr_matrix),std::ref(index_map), lib_size));
	}
	for (int i = 0; i < thread_holder.size(); i++){
		thread_holder[i].join();
	}

}

int main (int argc, char *argv[]){

	// usage
	// ./fast_count_multi num_threads gtf_file bam_files
	
	//gtf annotations to memory
	vector <vector <gtf_headers> > gtf_matrix = load_gtf(argv[2], atoi(argv[1]));
	
	//iterate through bam files
	for(int i = 3;i < argc;i++){
		
		//bam alignments chr and bp position in memory
		vector <vector <bam_headers> > bam_matrix = load_bam(argv[i]);
		
		//parse chromosome start sites to searchable 2d vector
		vector <vector <int> > chr_matrix = parse_chromosome(std::ref(bam_matrix));
		
		//map chromosome to bam matrix index
		vector <string> index_mapping = find_index(std::ref(bam_matrix));
				
		//launch threads
		//cout << "Performing counts for " << argv[i] << "\n";
		define_threads(std::ref(chr_matrix), atoi(argv[1]), std::ref(gtf_matrix), index_mapping, bam_matrix[0][0].lib_size);
		
	}
	//write counts
	write_counts(std::ref(gtf_matrix),argv,argc);
	
}



