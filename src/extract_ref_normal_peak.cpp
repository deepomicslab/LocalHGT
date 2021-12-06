#include <iostream>
#include <fstream>
#include <cmath> 
#include <cstdlib> 
#include <string>
#include <vector>
#include <time.h>
#include <cstdlib>
#include <ctime>
#include <string.h>
#include <sys/stat.h> 
#include <pthread.h>
#include <thread>
#include <sstream>
#include <map>


using namespace std;

const short c = 300;
const char coder_num = 3;
const char least_depth = 3;
const int k = 32;
long long array_size = pow(2, k);
char *kmer_count_table = new char[array_size];
// char* kmer_count_table = (char*)malloc(array_size);

int MIN_KMER_NUM = 6; //6
int REF_NEAR = 200;
int DIFF = 3; //3
int PEAK_W = 5; //5
int NEAR = PEAK_W; //PEAK_W

class Split_reads{
    public:
        map<int, int> chr_kmer_count;
        map<int, int> chr_peak_index;
        int min_kmer_num = MIN_KMER_NUM;
        void add_peak(int peak_chr, int peak_pos, int peak_index);
        void check_split(int* peak_filter);
};

void Split_reads::add_peak(int peak_chr, int peak_pos, int peak_index){
    if (chr_kmer_count.find(peak_chr) == chr_kmer_count.end()){
        chr_kmer_count[peak_chr] = 1;
        chr_peak_index[peak_chr] = peak_index;
    }
    else{
        chr_kmer_count[peak_chr] += 1;
    }   
}

void Split_reads::check_split(int* peak_filter){
    map<int, int>::iterator iter;
    int chr, peak_index;
    iter = chr_kmer_count.begin();
    while(iter != chr_kmer_count.end()){
        // cout << iter->first << ":"<< iter->second <<endl;
        if (iter->second < min_kmer_num){
            chr_kmer_count.erase(iter);
        }
        iter ++ ;
    }
    // cout <<"********************"<<endl;
    if (chr_kmer_count.size() > 1){
        // map<int, int>::iterator iter;
        iter = chr_kmer_count.begin();
        while(iter != chr_kmer_count.end()){
            // cout << iter->first << ":"<< iter->second <<endl;
            chr = iter->first;
            peak_index = chr_peak_index[chr];
            peak_filter[peak_index] = 1;
            iter ++ ;
        }
        // cout << "------------------"<<endl;
    }
    
}

class Peaks{
    public:
        int my_peak_index = 0;
        int filter_peak_num = 0;
        int near = NEAR;
        int ref_gap = 500;
        int ref_near = REF_NEAR;
        int max_peak_num = 500000000;
        int *peak_loci = new int[max_peak_num*2];
        int *peak_filter = new int[max_peak_num];
        unsigned int *peak_kmer = new unsigned int[array_size];
        // memset(peak_kmer, 0, sizeof(unsigned int)*array_size);
        void add_peak(int ref_index, int pos, unsigned int* record_ref_index,int ref_len);
        void delete_array(void);
        void slide_reads(string fastq_file, bool* coder, int* base, 
                        char* comple, short *choose_coder, int down_sam_ratio, long start, long end);
        void count_filtered_peak(string interval_name);
};

void Peaks::add_peak(int ref_index, int pos, unsigned int* record_ref_index, int ref_len){
    peak_loci[2*my_peak_index] = ref_index;
    peak_loci[2*my_peak_index+1] = pos;
    
    int index;
    for (int near_pos = pos - near; near_pos < pos + 1; near_pos++){
        if (near_pos>0 & near_pos<ref_len){
            for (int p = 0; p < 3; p++){
                index = coder_num*near_pos+p;
                peak_kmer[record_ref_index[index]] = my_peak_index;
            }  
        }
    }  
    if (my_peak_index > max_peak_num){
        cout <<"too many peaks!"<<endl;
    }
    my_peak_index += 1;  
}

void Peaks::slide_reads(string fastq_file, bool* coder, int* base, 
                        char* comple, short *choose_coder, int down_sam_ratio, long start, long end){

    time_t t0 = time(0);
    ifstream fq_file; 
    fq_file.open(fastq_file);
    string reads_seq;
    int reads_int [150];
    int reads_comple_int [150];
    unsigned int lines = 0;
    int converted_reads [450];
    int complemented_reads [450];
    int m;
    int n;
    unsigned int kmer_index, comple_kmer_index, real_index, b;   
    int r ;
    short read_len = 0;
    unsigned seed;
    seed = time(0);
    srand(seed);
    int chr_index, peak_locus,peak_index;

    long pos = 0;
    for (long i = start; i>0; i--){
        fq_file.seekg(i, ios::beg);
        char j;
        fq_file.get(j);
        if (j == '@'){ //only read name has this symbol.
            pos = i;
            break;
        }       
    }
    fq_file.seekg(pos, ios::beg);
    long add_size = start;
    // cout<<start<<"\t"<<end<<"\t||||||||||\t"<<pos<<endl;
    while (fq_file >> reads_seq)
    {
        if (lines % 1000000 == 1000000-1){
            cout << "recheck reads\t"<<lines<<endl;
        }
        if (add_size>=end){
            break;
        }
        add_size += reads_seq.length();
        // if (lines == 1){
        //     cout << reads_seq << endl;
        // }

        if (lines % 4 == 1){
            time_t t1 = time(0);
            if (lines == 1){
                read_len = reads_seq.length();//cal read length
            }
            r = rand() % 100;
            if (r < down_sam_ratio){
                for (int j = 0; j < read_len; j++){
                    reads_int[j] = (int)reads_seq[j];
                    reads_comple_int[j] = comple[reads_int[j]];
                }
                int peak_chr = 0;
                Split_reads each_read;
                for (int j = 0; j < read_len-k+1; j++){
                    for (int i = 0; i < 3; i++){
                        kmer_index = 0;
                        comple_kmer_index = 0;
                        bool all_valid = true;
                        for (int z = 0; z<k; z++){
                            m = coder[c*choose_coder[z*3+i]+reads_int[j+z]]; // choose_coder[z*3+i] indicate which coder
                            n = coder[c*choose_coder[(k-1-z)*3+i]+reads_comple_int[j+z]];
                            if (m == 5){
                                all_valid = false;
                                break;
                            }
                            kmer_index += m*base[z]; 
                            comple_kmer_index += n*base[(k-1-z)];
                              
                        }
                        if (kmer_index > comple_kmer_index){ //use a smaller index
                            real_index = comple_kmer_index;
                        }   
                        else{
                            real_index = kmer_index;
                        }
                        
                        if (peak_kmer[real_index] != 0){
                            each_read.add_peak(peak_loci[2*peak_kmer[real_index]], peak_loci[2*peak_kmer[real_index]+1], peak_kmer[real_index]);
                            // cout << peak_loci[2*peak_kmer[real_index]]<<"\t"<<peak_loci[2*peak_kmer[real_index]+1]<<endl;
                            // chr_index = peak_loci[2*peak_kmer[real_index]];
                            // peak_locus = peak_loci[2*peak_kmer[real_index]+1];
                            // if (peak_chr != 0 & chr_index != peak_chr){
                            //     peak_filter[peak_index] = 1;
                            //     peak_filter[peak_kmer[real_index]] = 1;
                            //     // cout << peak_chr<<"\t"<< chr_index<< "\t"<<peak_locus<<endl;
                            // }
                            // peak_index = peak_kmer[real_index];
                            // peak_chr = chr_index;
                            
                        }
                        
                    }
                }
                // cout << "****************"<<endl;
                each_read.check_split(peak_filter);
            }         
        }
        lines++;
    }
    fq_file.close();
}

void Peaks::delete_array(void){
    delete [] peak_loci;
    delete [] peak_kmer;
    delete [] peak_filter;
}

void Peaks::count_filtered_peak(string interval_name){
    ofstream interval_file;
    interval_file.open(interval_name, ios::out | ios::trunc);
    int start = 1;
    int end = 1;
    int chr = 1;
    long final_ref_len = 0;
    for (int i = 0; i < my_peak_index; i++){
        if (peak_filter[i] == 1 ){
            filter_peak_num += 1;
            if (chr == peak_loci[2*i] & peak_loci[2*i+1]-ref_near - end < ref_gap){
                end = peak_loci[2*i+1]+ref_near;
            }
            else{
                interval_file << chr << "\t" << start << "\t"<< end << endl;
                final_ref_len += (end - start);
                chr = peak_loci[2*i];
                start = peak_loci[2*i+1]-ref_near;
                end = peak_loci[2*i+1]+ref_near;
            }
            
        }
    }
    interval_file << chr << "\t" << start << "\t"<< end << endl;
    final_ref_len += (end - start);
    interval_file.close();
    cout <<"##########final_ref_len"<<final_ref_len<<endl;
}

long slide_window(bool* record_ref_hit, int ref_len, int ref_index, long extract_ref_len,
                ofstream & interval_file, float hit_ratio, float perfect_hit_ratio,
                long & total_peak_num,Peaks & MyPeak,unsigned int* record_ref_index){ 
                //find density hits regions
    int frag_index = 0;
    int window = 500;
    int one_coder_bases = 0; // the number of bases supported by at least one coder in a window
    int three_coder_bases = 0; // the number of bases supported by three coder in a window
    int one_coder_min = window * hit_ratio;
    int three_coder_min = window * perfect_hit_ratio;    
    bool good_window;
    bool conti_flag = false;
    int start = 0 ;
    int end = 0;
    int* save_good_intervals = new int[2*ref_len/window];
    short *single_hit_num = new short[ref_len];
    short *trio_hit_num = new short[ref_len];
    bool *peak_hit = new bool[ref_len];
    int* save_peak_intervals = new int[2*ref_len];
    int peak_index = 0;
    
    for (int j = 0; j < ref_len; j++){
        short hit_coder_num = 0;
        for (int p = 0; p < 3; p++){
            // if (record_ref_hit[ref_len*p+j] == true){
            if (record_ref_hit[coder_num*j+p] == true){
                hit_coder_num += 1;
                // cout << j <<"hit"<<p << endl;
            }
        }
        
        if (hit_coder_num == 3){
            trio_hit_num[j] = 1;
        }
        else{
            trio_hit_num[j] = 0;
        }
        if (hit_coder_num > 0){
            single_hit_num[j] = 1;
        }
        else{
            single_hit_num[j] = 0;
        }
        
        if (j < window){
            if (hit_coder_num >= 1){
                one_coder_bases += 1;
            }
            if (hit_coder_num == 3){
                three_coder_bases += 1;
            }
        }
        else{
            one_coder_bases = one_coder_bases - single_hit_num[j-window] + single_hit_num[j];
            three_coder_bases = three_coder_bases - trio_hit_num[j-window] + trio_hit_num[j];
        }

        if (one_coder_bases >= one_coder_min & three_coder_bases >= three_coder_min){
            good_window = true;
        }
        else{
            good_window = false;
        }
        
        if (conti_flag == false & good_window == true){
            start = j - 2 * window;
            if (start < 1){
                start = 1;
            }
            conti_flag = true;
        }
        if (conti_flag == true & good_window == false){
            end = j + 2 * window;
            if (end > ref_len){
                end = ref_len;
            }
            if (frag_index > 0 & start - save_good_intervals[2*frag_index-1] < window ){
                save_good_intervals[2*frag_index-1] = end;
            }
            else{
                save_good_intervals[2 * frag_index] = start;
                save_good_intervals[2 * frag_index+1] = end;
                frag_index += 1;
            }
            conti_flag = false;
        }

        // find peak
        int w = PEAK_W;
        peak_hit[j] = false;
        if (j > 60){
            for (int m = 0; m < 5; m++){
                int diff = 0;
                for (int n = 0; n < w; n++){
                    diff += (single_hit_num[j-m-w-n] - single_hit_num[j-n]);
                }
                // original
                // if (abs(diff) >= 3){
                //     peak_hit[j] = true;
                //     break;
                // }
                if (diff >= DIFF){
                    peak_hit[j-m-w] = true;
                }
                if (diff <= -DIFF){
                    peak_hit[j] = true;
                }
            }
            
        }
        // find peak

    }
    if (conti_flag == true & good_window == true){
        end = ref_len;
        if (frag_index > 0 & start - save_good_intervals[2*frag_index-1] < window ){
            save_good_intervals[2*frag_index-1] = end;
        }
        else{
            save_good_intervals[2 * frag_index] = start;
            save_good_intervals[2 * frag_index+1] = end;
            frag_index += 1;
        }

    }

    for (int i = 0; i < frag_index; i++){
        // extract_ref_len += (save_good_intervals[2*i+1] - save_good_intervals[2*i]);
        for (int j =save_good_intervals[2*i]; j <save_good_intervals[2*i+1];j++ ){
            if (peak_hit[j] == true){
                total_peak_num += 1;
                MyPeak.add_peak(ref_index, j, record_ref_index, ref_len);
                start = j - 200;
                end = j + 200;
                if (end > ref_len){
                    end = ref_len;
                }
                if (start < 1){
                    start = 1;
                }
                if (peak_index > 0 & start - save_peak_intervals[2*peak_index-1] < window ){
                    save_peak_intervals[2*peak_index-1] = end;
                }
                else{
                    save_peak_intervals[2 * peak_index] = start;
                    save_peak_intervals[2 * peak_index+1] = end;
                    peak_index += 1;
                }

            }
        }
        
    }

    for (int i = 0; i < peak_index; i++){
        extract_ref_len += (save_peak_intervals[2*i+1] - save_peak_intervals[2*i]);
        // interval_file << ref_index << "\t" << save_peak_intervals[2*i] << "\t"<< save_peak_intervals[2*i+1] << endl;
    }
    // for (int i = 0; i < frag_index; i++){
    //     extract_ref_len += (save_good_intervals[2*i+1] - save_good_intervals[2*i]);
    //     interval_file << ref_index << "\t" << save_good_intervals[2*i] << "\t"<< save_good_intervals[2*i+1] << endl;
    // }
    // cout << frag_index << endl;
    delete [] save_good_intervals;
    delete [] single_hit_num;
    delete [] trio_hit_num;
    delete [] save_peak_intervals;
    delete [] peak_hit;
    // return save_good_intervals;
    return extract_ref_len;
}

float cal_tab_empty_rate(){
    double empty_num = 0;
    float empty_rate;
    float weak_rate = 0;
    double weak_num = 0;
    for (unsigned int j = 0; j < array_size; j++){  
        if ((int)kmer_count_table[j] == 0){
            empty_num += 1;
        }
        // else{
        //     cout << (int)kmer_count_table[j] <<empty_num<< endl;
        // }
        if ((int)kmer_count_table[j] != least_depth){
            weak_num += 1;
        }
        // cout << j << (int)kmer_count_table[j]<< empty_num <<endl;

        
    }
    empty_rate = empty_num/array_size;
    weak_rate = weak_num/array_size;
    cout << array_size << "\t" << weak_num<< "\t" << empty_num<<endl;
    cout << array_size << "\t" << weak_rate<< "\t" << empty_rate<<endl;
    return empty_rate;
}

void read_ref(string fasta_file, bool* coder, int* base, int k, char* comple,
             string index_name, short *choose_coder)
{
    ifstream fa_file;
    ofstream index_file;
    index_file.open(index_name, ios::out | ios::binary);
    fa_file.open(fasta_file, ios::in);
    string ref_seq, line_seq;
    ref_seq = "\0";
    int ref_len;
    char n;
    int m;
    int e;
    unsigned int kmer_index, comple_kmer_index, real_index;
    int ref_index = 0;
    long extract_ref_len = 0;
    long slide_ref_len = 0;
    string chr_name ;
    string pre_name = "start";
    time_t t0 = time(0);
    int covert_num, comple_num;
    short convert_ref[150];
    short complemented_ref[150];
    char support_coder;

    // save random coder
    for (int j = 0; j < 100; j++){
        index_file.write((char *)(&choose_coder[j]), sizeof(unsigned int));
    }
    cout <<"Start index ref..."<<endl;

    while (fa_file >> line_seq){ //remember the last chr
        if (line_seq[0] == '>'){
            chr_name = pre_name;
            pre_name = line_seq.substr(1);
            // if (ref_index > 10){
            //     break;
            // }
            ref_len= ref_seq.length();
            slide_ref_len += ref_len;
            // cout << chr_name <<"\t"<< ref_index <<"\t" << ref_len << endl;
            if (ref_len > k){   
                for (int i = 0; i < 3; i++){
                    kmer_index = 0; // start kmer
                    comple_kmer_index = 0;
                }
                int *ref_int = new int[ref_len];
                int *ref_comple_int = new int[ref_len];
                for (int j = 0; j < ref_len; j++){
                    ref_int[j] = (int)ref_seq[j];
                    ref_comple_int[j] = comple[ref_int[j]];
                }
                index_file.write((char *)(&ref_len), sizeof(unsigned int));
                for (int j = 0; j < ref_len-k+1; j++){
                    support_coder = 0;
                    for (int i = 0; i < 3; i++){
                        kmer_index = 0;
                        comple_kmer_index = 0;
                        bool all_valid = true;
                        for (int z = 0; z<k; z++){
                            m = coder[c*choose_coder[z*3+i]+ref_int[j+z]];
                            if (m == 5){
                                all_valid = false;
                                break;
                            }
                            kmer_index += m*base[z]; 
                            comple_kmer_index += coder[c*choose_coder[(k-1-z)*3+i]+ref_comple_int[j+z]]*base[(k-1-z)];  
                        }
                        // cout <<kmer_index<<"\t"<< j<<"\t"<<comple_kmer_index<<"\t"<<reads_seq[j] <<endl;
                        if (kmer_index > comple_kmer_index){ //use a smaller index
                            real_index = comple_kmer_index;
                        }   
                        else{
                            real_index = kmer_index;
                        }
                        if (all_valid == false){
                            real_index = 0;
                        }
                        index_file.write((char *)(&real_index), sizeof(real_index));
                        
                    }
                }
                
                // cout << chr_name<<"\t"<<endl;
                delete [] ref_int;
                delete [] ref_comple_int;
            }
            if (ref_index % 1000 == 0){
                time_t t1 = time(0);
                cout << chr_name<<"\t" << ref_index << "\t" <<ref_len <<" bp\t"<<slide_ref_len<<" bp\t" <<t1-t0<< endl;
            }  
            ref_index += 1;
            ref_seq = "\0";
        }
        else{
            ref_seq += line_seq;           
        }            
    }
    // for the last chr
    chr_name = pre_name;
    ref_len= ref_seq.length();
    slide_ref_len += ref_len;
    if (ref_len > k){   
        for (int i = 0; i < 3; i++){
            kmer_index = 0; // start kmer
            comple_kmer_index = 0;
        }
        int *ref_int = new int[ref_len];
        int *ref_comple_int = new int[ref_len];
        for (int j = 0; j < ref_len; j++){
            ref_int[j] = (int)ref_seq[j];
            ref_comple_int[j] = comple[ref_int[j]];
        }
        index_file.write((char *)(&ref_len), sizeof(unsigned int));
        for (int j = 0; j < ref_len-k+1; j++){
            support_coder = 0;
            for (int i = 0; i < 3; i++){
                kmer_index = 0;
                comple_kmer_index = 0;
                bool all_valid = true;
                for (int z = 0; z<k; z++){
                    m = coder[c*choose_coder[z*3+i]+ref_int[j+z]];
                    if (m == 5){
                        all_valid = false;
                        break;
                    }
                    kmer_index += m*base[z]; 
                    comple_kmer_index += coder[c*choose_coder[(k-1-z)*3+i]+ref_comple_int[j+z]]*base[(k-1-z)];  
                }
                if (kmer_index > comple_kmer_index){ //use a smaller index
                    real_index = comple_kmer_index;
                }   
                else{
                    real_index = kmer_index;
                }
                if (all_valid == false){
                    real_index = 0;
                }
                index_file.write((char *)(&real_index), sizeof(real_index));  
            }
        }        
        delete [] ref_int;
        delete [] ref_comple_int;
        time_t t1 = time(0);
        cout << "Last chr:\t" << chr_name<<"\t" << ref_index << "\t" <<ref_len <<" bp\t"<<slide_ref_len<<" bp\t" <<t1-t0<< endl;
    }

    cout << "Index is done."<< "\t" << extract_ref_len << endl;
    fa_file.close();
    index_file.close();
}

void read_index(bool* coder, int* base, int k, char* comple, string index_name, string interval_name,
                short *choose_coder, float hit_ratio, float perfect_hit_ratio, Peaks & MyPeak){
    bool *record_ref_hit; 
    unsigned int *record_ref_index;
    long total_peak_num = 0;

    ifstream index_file;
    ofstream interval_file;
    index_file.open(index_name, ios::in | ios::binary); 
    interval_file.open(interval_name, ios::out | ios::trunc);
    unsigned int real_index;
    unsigned int kmer_index;
    long i = 0;
    int ref_len, ref_index;
    long ref_end = 0;
    ref_index = 1;
    long ref_start = 0;
    long extract_ref_len = 0;
    long slide_ref_len = 0;
    // unsigned char *record_ref_hit; 
    time_t t0 = time(0);
 
    long long start_point = 400;
    int buffer_size;
    cout <<"Start slide ref..."<<endl;

    index_file.seekg (0, index_file.end);
    long long length = index_file.tellg();

    index_file.seekg(400, ios::beg);
    while (!index_file.eof()){
        index_file.read(reinterpret_cast<char*>(&real_index), sizeof(unsigned int));
        ref_len = real_index;
        buffer_size = (ref_len-k+1)*coder_num*4;
        char *each_ref_buffer = new char[buffer_size];
        index_file.read(each_ref_buffer, buffer_size);

        int n = 0;
        int j = 0;
        record_ref_hit = new bool[ref_len*coder_num];
        record_ref_index = new unsigned int [ref_len*coder_num];
        while (n < buffer_size){
            memcpy(&kmer_index, each_ref_buffer, sizeof(unsigned int));
            record_ref_index[j] = kmer_index;
            if ((int)kmer_count_table[kmer_index] == least_depth & real_index != 0){
                record_ref_hit[j] = true;
                // cout <<(int)kmer_count_table[real_index] << "\t"<<i<<endl;
            }
            else{
                record_ref_hit[j] = false;
            }
            j += 1;
            each_ref_buffer = each_ref_buffer + sizeof(unsigned int);
            n = n + sizeof(unsigned int);
        }
        extract_ref_len = slide_window(record_ref_hit, ref_len, ref_index, extract_ref_len, interval_file, hit_ratio, 
                                        perfect_hit_ratio,total_peak_num,MyPeak,record_ref_index);
        slide_ref_len += ref_len;
        start_point += 4;
        start_point += buffer_size;
        each_ref_buffer = each_ref_buffer - buffer_size; //move the pointer to the start
        delete [] each_ref_buffer;
        delete [] record_ref_hit;
        delete [] record_ref_index;
        ref_index += 1;
        if (start_point >= length){
            cout << "read index ending." <<endl;
            break;
        }
        if (ref_index % 1000 == 0){
            time_t t1 = time(0);
            cout << ref_index << "\t" <<slide_ref_len << " bp\t" << extract_ref_len <<" bp\t" <<t1-t0<< endl;
        } 
        
        i += 1;
        real_index = 0;
    }       
    index_file.close();  
    interval_file.close();
    time_t t1 = time(0);
    cout << slide_ref_len << " bp\t" << extract_ref_len <<" bp\t" <<total_peak_num<<"\t" <<t1-t0<< endl;
}

// vector<char> 
void read_fastq(string fastq_file, int k, bool* coder, int* base, char* comple, 
    short *choose_coder, int down_sam_ratio, long start, long end)
{
    time_t t0 = time(0);
    ifstream fq_file; 
    fq_file.open(fastq_file);

    long pos = 0;
    for (long i = start; i>0; i--){
        fq_file.seekg(i, ios::beg);
        char j;
        fq_file.get(j);
        if (j == '@'){ //only read name has this symbol.
            pos = i;
            break;
        }       
    }
    fq_file.seekg(pos, ios::beg);
    long add_size = start;


    string reads_seq;
    int reads_int [150];
    int reads_comple_int [150];

    // unsigned int i = 0;
    unsigned int lines = 0;
    int converted_reads [450];
    int complemented_reads [450];
    int m;
    int n;
    unsigned int kmer_index, comple_kmer_index, real_index, b;   
    int r ;
    short read_len = 0;
    // char abnormal_base[150];
    unsigned seed;
    seed = time(0);
    srand(seed);

    while (fq_file >> reads_seq)
    {
        if (add_size>=end){
            break;
        }
        add_size += reads_seq.length();

        if (lines % 4 == 1){
            time_t t1 = time(0);
            if (lines % 10000000 == 10000000-1){
                cout <<lines<<"reads\t" << t1-t0 <<endl;
            }
            if (lines == 1){
                read_len = reads_seq.length();//cal read length
            }
            // srand((unsigned)time(NULL));
            r = rand() % 100 ;
            // cout <<r << "r"<<endl;
            if (r < down_sam_ratio){
                for (int j = 0; j < read_len; j++){
                    reads_int[j] = (int)reads_seq[j];
                    reads_comple_int[j] = comple[reads_int[j]];
                }
                for (int j = 0; j < read_len-k+1; j++){
                    for (int i = 0; i < 3; i++){
                        kmer_index = 0;
                        comple_kmer_index = 0;
                        bool all_valid = true;
                        // cout <<j<<"\t"<<i<<"\t";
                        for (int z = 0; z<k; z++){
                            m = coder[c*choose_coder[z*3+i]+reads_int[j+z]]; // choose_coder[z*3+i] indicate which coder
                            n = coder[c*choose_coder[(k-1-z)*3+i]+reads_comple_int[j+z]];
                            // cout <<m<<" ";
                            if (m == 5){
                                all_valid = false;
                                break;
                            }
                            kmer_index += m*base[z]; 
                            comple_kmer_index += n*base[(k-1-z)];
                              
                        }
                        // cout<<endl;
                        // cout <<kmer_index<<"\t"<< j<<"\t"<<comple_kmer_index<<"\t"<<reads_seq[j] <<endl;
                        if (kmer_index > comple_kmer_index){ //use a smaller index
                            real_index = comple_kmer_index;
                        }   
                        else{
                            real_index = kmer_index;
                        }
                        if ((int)kmer_count_table[real_index] < least_depth & all_valid == true ){
                            kmer_count_table[real_index] += 1;
                            // cout << (int)kmer_count_table[real_index] << "\t" << endl;
                        }  
                    }
                }
            }         
        }
        lines++;
    }
    fq_file.close();
    // return kmer_count_table;

}
// */

bool * generate_coder(bool coder_num)
{
    // A:65 97 T:116 84 C:99 67 G: 103 71
    static bool coder [1000];
    for (int j = 0; j < 1000; j++){
        coder[j] = 5;
    }
    coder[65] = 1;
    coder[97] = 1;

    coder[84] = 1;
    coder[116] = 1;

    coder[67] = 0;
    coder[99] = 0;

    coder[71] = 0;
    coder[103] = 0;

    coder[65+c] = 1;
    coder[97+c] = 1;

    coder[84+c] = 0;
    coder[116+c] = 0;

    coder[67+c] = 1;
    coder[99+c] = 1;

    coder[71+c] = 0;
    coder[103+c] = 0;

    coder[65+2*c] = 1;
    coder[97+2*c] = 1;

    coder[84+2*c] = 0;
    coder[116+2*c] = 0;

    coder[67+2*c] = 0;
    coder[99+2*c] = 0;

    coder[71+2*c] = 1;
    coder[103+2*c] = 1;

    return coder;
}

int * generate_base(int k)
{
    static int base [32];
    for (int i = 0; i<k; i++){       
        base[i] = pow(2, k-i-1);
    }
    return base;
}

char * generate_complement(void)
{
    static char comple [256];
    for (int j = 0; j < 256; j++){
        comple[j] = 0;
    }
    comple[65] = 84;
    comple[97] = 84;
    comple[116] = 65;
    comple[84] = 65;
    comple[99] = 71;
    comple[67] = 71;
    comple[103] = 67;
    comple[71] = 67;
    return comple;
}

short * random_coder(int k)
{
    static short permu[18] = {0,1,2,0,2,1,1,2,0,1,0,2,2,0,1,2,1,0};
    static short choose_coder[100];
    int r;
    unsigned seed;
    seed = time(0);
    srand(seed);
    for (int j = 0; j < k; j++){
        r = rand() % 6;
        // cout << r << endl;
        for (int i = 0; i < coder_num; i++){
            choose_coder[j*3+i] = permu[r*3+i];
            // cout << choose_coder[j*3+i]<<"#"<<permu[r*3+i] << endl;
        }
    }
    // for (int j = 0; j < 100; j++){
    //         cout << j <<"\t" << choose_coder[j] << endl;
    // }
    return choose_coder;
}

short * saved_random_coder(string index_name){
    ifstream index_file;
    static short read_choose_coder[100];
    long i = 0;
    unsigned int real_index;
    index_file.open(index_name, ios::in | ios::binary);   
    while (!index_file.eof()){
        index_file.read(reinterpret_cast<char*>(&real_index), sizeof(unsigned int));
        if (i < 100){
            read_choose_coder[i] = real_index;
        }
        else{
            break;
        }
        i += 1;  
    }
    index_file.close(); 
    return read_choose_coder;
}

int cal_sam_ratio(string fq1, long down_sampling_size){
    int down_sam_ratio = 0;
    int read_len = 0;
    long i = 0;
    long sample_size = 0;
    string reads_seq;

    ifstream fq_file; 
    fq_file.open(fq1);
    while (fq_file >> reads_seq){
        if (i % 4 == 1 & read_len == 0){
            read_len = reads_seq.length();
            cout <<"read length is "<<read_len<<endl;
        }
        i += 1;
    }
    fq_file.close();
    sample_size = (i/2)*read_len;
    down_sam_ratio = 100*down_sampling_size/sample_size;
    cout<<i<<"\t"<<sample_size<<" Down-sampling ratio is "<<down_sam_ratio<< "%." <<endl;
    return down_sam_ratio;
}

long file_size(string filename)  
{  
    struct stat statbuf;  
    stat(filename.c_str(),&statbuf);  
    long size=statbuf.st_size;   
    return size;  
}  

int main( int argc, char *argv[])
{
    bool *coder;
    int *base;
    char *comple;
    short *choose_coder;
    coder = generate_coder(3);
    base = generate_base(k);
    comple = generate_complement();
    time_t now1 = time(0);

    string fasta_file = argv[3];
    string fq1 = argv[1];
    string fq2 = argv[2];   
    string interval_name = argv[4];
    string accept_hit_ratio = argv[5];
    string accept_perfect_hit_ratio = argv[6];
    float hit_ratio = stod(accept_hit_ratio);
    float perfect_hit_ratio = stod(accept_perfect_hit_ratio);
    long down_sampling_size = 2000000000; //2G bases

    int thread_num = 5;
    long start = 0;
    long end = 0;

    // int down_sam_ratio = cal_sam_ratio(fq1, down_sampling_size); //percent of downsampling ratio (1-100).
    int down_sam_ratio = 100;
    //index
    string index_name = fasta_file + ".index.dat";
    ifstream findex(index_name);
    if (! findex.good()){
        cout << "Reference index not detected, start index..." << endl;
        choose_coder = random_coder(k); 
        read_ref(fasta_file, coder, base, k, comple, index_name, choose_coder);
    }
    else{
        cout << "Reference index is detected." << endl;
    }

    // start
    cout << "Start extract..."<<endl;
    memset(kmer_count_table, 0, sizeof(char)*array_size);
    choose_coder = saved_random_coder(index_name);

    //multithread

    long size = file_size(fq1);
    long each_size = size/thread_num;
    cout <<size<<endl;
    
    std::vector<std::thread>threads;

    for (int i=0; i<thread_num; i++){
        start = i*each_size;
        end = (i+1)*each_size;
        if (i == thread_num-1){
            end = size;
        }
        threads.push_back(thread(read_fastq, fq1, k, coder, base, comple, choose_coder, down_sam_ratio, start, end));
    }
	for (auto&th : threads)
		th.join();
    threads.clear();

    
    for (int i=0; i<thread_num; i++){
        start = i*each_size;
        end = (i+1)*each_size;
        if (i == thread_num-1){
            end = size;
        }
        threads.push_back(thread(read_fastq, fq2, k, coder, base, comple, choose_coder, down_sam_ratio, start, end));
    }
	for (auto&th : threads)
		th.join();
    threads.clear();
    time_t now2 = time(0);
    cout << "reads finish.\t" << now2 - now1 << endl;
    Peaks MyPeak;
    read_index(coder, base, k, comple, index_name, interval_name, choose_coder, 
                hit_ratio, perfect_hit_ratio,MyPeak);

    // down_sam_ratio = 13;

    for (int i=0; i<thread_num; i++){
        start = i*each_size;
        end = (i+1)*each_size;
        if (i == thread_num-1){
            end = size;
        }
        threads.push_back(thread(&Peaks::slide_reads, MyPeak, fq1, coder, base, comple, choose_coder, down_sam_ratio, start, end));
    }
	for (auto&th : threads)
		th.join();
    threads.clear();

    for (int i=0; i<thread_num; i++){
        start = i*each_size;
        end = (i+1)*each_size;
        if (i == thread_num-1){
            end = size;
        }
        // cout <<start<<"\t"<<end<<endl;
        threads.push_back(thread(&Peaks::slide_reads, MyPeak, fq2, coder, base, comple, choose_coder, down_sam_ratio, start, end));
    }
	for (auto&th : threads)
		th.join();
    threads.clear();



    // MyPeak.slide_reads(fq1, coder, base, comple, choose_coder, down_sam_ratio, 0, 2000000000);
    // MyPeak.slide_reads(fq2, coder, base, comple, choose_coder, down_sam_ratio, 0, 2000000000);
    MyPeak.count_filtered_peak(interval_name);
    cout<< MyPeak.filter_peak_num<<"---------peak number:"<<MyPeak.my_peak_index<<endl;
    MyPeak.delete_array();
    delete [] kmer_count_table;
    time_t now3 = time(0);
    cout << "Finish with time:\t" << now3-now1<<endl;
    return 0;
}


