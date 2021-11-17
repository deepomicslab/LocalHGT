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


using namespace std;

const char coder_num = 1;
const char least_depth = 127;
const int k = 16;
long long array_size = pow(4, k);
char *kmer_count_table = new char[array_size];
// vector <char>kmer_count_table(array_size);
// char* kmer_count_table = (char*)malloc(array_size);


long slide_window(char* record_ref_hit, int ref_len, long ref_index, long extract_ref_len, ofstream & interval_file, float hit_ratio, float perfect_hit_ratio){ //find density hits regions
    int frag_index = 0;
    int window = 50;
    int three_coder_bases = 0; // the number of bases supported by three coder in a window
    int three_coder_min = window * perfect_hit_ratio;    
    bool good_window;
    bool conti_flag = false;
    int start = 0 ;
    int end = 0;
    int* save_good_intervals = new int[2*ref_len/window];
    short *trio_hit_num = new short[ref_len];
    bool peak_flag = false;
    int peak_win = 10;
    int diff = 0;

    for (int j = 0; j < ref_len; j++){
        if ((int)record_ref_hit[j] >= 3){
            trio_hit_num[j] = 1;
        }
        else{
            trio_hit_num[j] = 0;
        }

        peak_flag = false;
        if (j > 2* peak_win){
            diff = 0;
            for (int i = 0; i < peak_win; i++){
                diff = diff + ((int)record_ref_hit[j-i] - (int)record_ref_hit[j-peak_win-i]);
                // cout << (int)record_ref_hit[j-i] << "\t" << (int)record_ref_hit[j-peak_win-i] << endl;
            }
            if (abs(diff)/peak_win > 10){
                peak_flag = true;
                // cout << diff << "\t" << abs(diff)/peak_win<< "\t"<<(int)record_ref_hit[j]<<endl;
            }
        }
        // cout << (int)record_ref_hit[j]<< "\t" << peak_flag << endl;
        
        if (j < window){
            if ((int)record_ref_hit[j] >= 3){
                three_coder_bases += 1;
            }
        }
        else{
            three_coder_bases = three_coder_bases - trio_hit_num[j-window] + trio_hit_num[j];
        }

        if (three_coder_bases >= three_coder_min & peak_flag){
            good_window = true;
        }
        else{
            good_window = false;
        }  
        // cout << three_coder_bases << "\t" <<peak_flag<< "\t" << good_window<<endl;      
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
        extract_ref_len += (save_good_intervals[2*i+1] - save_good_intervals[2*i]);
        interval_file << ref_index << "\t" << save_good_intervals[2*i] << "\t"<< save_good_intervals[2*i+1] << endl;
    }
    // cout << frag_index << endl;
    delete [] save_good_intervals;
    delete [] trio_hit_num;
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

void find_ref(string fasta_file, short* coder, int* base, int k, char* comple, string interval_name, float hit_ratio, float perfect_hit_ratio)
{
    ofstream interval_file;
    interval_file.open(interval_name, ios::out | ios::trunc);
    ifstream fa_file;
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
    char *record_ref_hit; 

    // save random coder
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
                record_ref_hit = new char[ref_len*coder_num];
                kmer_index = 0; // start kmer
                comple_kmer_index = 0;
                int *ref_int = new int[ref_len];
                int *ref_comple_int = new int[ref_len];
                for (int j = 0; j < ref_len; j++){
                    ref_int[j] = (int)ref_seq[j];
                    ref_comple_int[j] = comple[ref_int[j]];
                }
                kmer_index = 0;
                comple_kmer_index = 0;
                for (int j = 0; j < ref_len-k+1; j++){
                    bool all_valid = true;
                    if (coder[ref_int[j]] == 5){
                        all_valid = false;
                    }
                    if (j == 0){
                        for (int z = 0; z<k; z++){
                            kmer_index += coder[ref_int[j+z]]*base[z];  
                            comple_kmer_index += coder[ref_comple_int[j+z]]*base[(k-1-z)];
                        }
                    }
                    else{
                        kmer_index = (kmer_index - coder[ref_int[j-1]]*base[0])*4 + coder[ref_int[j+k-1]];
                        comple_kmer_index = (comple_kmer_index - coder[ref_comple_int[j-1]])/4 + coder[ref_comple_int[j+k-1]]*base[0];
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
                    // if ((int)kmer_count_table[real_index] == least_depth & real_index != 0){
                    //     record_ref_hit[j] = true;
                        
                    // }
                    // else{
                    //     record_ref_hit[j] = false;
                    //     // cout << i <<endl;
                    // }
                    if (real_index != 0){
                        record_ref_hit[j] = kmer_count_table[real_index];
                    }
                    else{
                        record_ref_hit[j] = 0;
                    }
                }
                extract_ref_len = slide_window(record_ref_hit, ref_len, ref_index, extract_ref_len, interval_file, hit_ratio, perfect_hit_ratio);
                
                // cout << chr_name<<"\t"<<endl;
                delete [] ref_int;
                delete [] ref_comple_int;
                delete [] record_ref_hit;
            }
            if (ref_index % 1000 == 0){
                time_t t1 = time(0);
                cout << chr_name<<"\t" << ref_index << "\t" <<ref_len <<" bp\t"<<slide_ref_len<< " bp\t"<<extract_ref_len<<" bp\t" <<t1-t0<< endl;
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
        record_ref_hit = new char[ref_len];  
        int *ref_int = new int[ref_len];
        int *ref_comple_int = new int[ref_len];
        for (int j = 0; j < ref_len; j++){
            ref_int[j] = (int)ref_seq[j];
            ref_comple_int[j] = comple[ref_int[j]];
        }
        kmer_index = 0;
        comple_kmer_index = 0;
        for (int j = 0; j < ref_len-k+1; j++){
            bool all_valid = true;
            if (coder[ref_int[j]] == 5){
                all_valid = false;
            }
            if (j == 0){
                for (int z = 0; z<k; z++){
                    kmer_index += coder[ref_int[j+z]]*base[z];  
                    comple_kmer_index += coder[ref_comple_int[j+z]]*base[(k-1-z)];
                }
            }
            else{
                kmer_index = (kmer_index - coder[ref_int[j-1]]*base[0])*4 + coder[ref_int[j+k-1]];
                comple_kmer_index = (comple_kmer_index - coder[ref_comple_int[j-1]])/4 + coder[ref_comple_int[j+k-1]]*base[0];
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

            if ((int)kmer_count_table[real_index] == least_depth & real_index != 0){
                record_ref_hit[j] = true;
                
            }
            else{
                record_ref_hit[j] = false;
            }

        }     
        extract_ref_len = slide_window(record_ref_hit, ref_len, ref_index, extract_ref_len, interval_file, hit_ratio, perfect_hit_ratio);   
        delete [] ref_int;
        delete [] ref_comple_int;
        delete [] record_ref_hit;
        time_t t1 = time(0);
        cout << "Last chr:\t" << chr_name<<"\t" << ref_index << "\t" <<ref_len <<" bp\t"<<slide_ref_len<< " bp\t"<<extract_ref_len<<" bp\t" <<t1-t0<< endl;
    }

    cout << "Extracted ref length is"<< "\t" << extract_ref_len << endl;
    fa_file.close();
}

// vector<char> 
void read_fastq(string fastq_file, int k, short* coder, int* base, char* comple, int down_sam_ratio, long start, long end)
{
    time_t t0 = time(0);
    ifstream fq_file; 
    fq_file.open(fastq_file);

    long pos;
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

    unsigned int i = 0;
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

        if (i % 4 == 1){
            time_t t1 = time(0);
            if (i % 10000000 == 1){
                cout <<i<<"reads\t" << t1-t0 <<endl;
            }
            if (i == 1){
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
                kmer_index = 0;
                comple_kmer_index = 0;
                for (int j = 0; j < read_len-k+1; j++){

                    bool all_valid = true;

                    if (coder[reads_int[j]] == 5){
                        all_valid = false;
                    }

                    if (j == 0){
                        for (int z = 0; z<k; z++){
                            kmer_index += coder[reads_int[j+z]]*base[z];  
                            comple_kmer_index += coder[reads_comple_int[j+z]]*base[(k-1-z)];
                        }
                    }
                    else{
                        kmer_index = (kmer_index - coder[reads_int[j-1]]*base[0])*4 + coder[reads_int[j+k-1]];
                        comple_kmer_index = (comple_kmer_index - coder[reads_comple_int[j-1]])/4 + coder[reads_comple_int[j+k-1]]*base[0];
                    }


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
        i++;
    }
    fq_file.close();
    // return kmer_count_table;

}
// */

short * generate_coder()
{
    // A:65 97 T:116 84 C:99 67 G: 103 71
    static short coder [256];
    for (int j = 0; j < 256; j++){
        coder[j] = 5;
    }
    coder[65] = 0;
    coder[97] = 0;

    coder[84] = 1;
    coder[116] = 1;

    coder[67] = 2;
    coder[99] = 2;

    coder[71] = 3;
    coder[103] = 3;

    return coder;
}

int * generate_base(int k)
{
    static int base [16];
    for (int i = 0; i<k; i++){       
        base[i] = pow(4, k-i-1);
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
    // string ref_seq = read_ref();
    short *coder;
    int *base;
    char *comple;
    coder = generate_coder();
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

    int down_sam_ratio = cal_sam_ratio(fq1, down_sampling_size); //percent of downsampling ratio (1-100).
    // int down_sam_ratio = 100;
    // start
    cout << "Start all-matching extract..."<<endl;
    memset(kmer_count_table, 0, sizeof(char)*array_size);


    //multithread

    long size = file_size(fq1);
    long each_size = size/thread_num;
    // cout <<size<<endl;
    
    std::vector<std::thread>threads;

    for (int i=0; i<thread_num; i++){
        start = i*each_size;
        end = (i+1)*each_size;
        if (i == thread_num-1){
            end = size;
        }
        cout <<start<<"\t"<<end<<endl;
        threads.push_back(thread(read_fastq, fq1, k, coder, base, comple, down_sam_ratio, start, end));
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
        cout <<start<<"\t"<<end<<endl;
        threads.push_back(thread(read_fastq, fq2, k, coder, base, comple, down_sam_ratio, start, end));
    }
	for (auto&th : threads)
		th.join();
    threads.clear();
    time_t now2 = time(0);
    cout << "reads finish.\t" << now2 - now1 << endl;


    // read_index(coder, base, k, comple, index_name, interval_name, hit_ratio, perfect_hit_ratio);
    find_ref(fasta_file, coder, base, k, comple, interval_name, hit_ratio, perfect_hit_ratio);
    time_t now3 = time(0);
    cout << "Finish with time:\t" << now3-now1<<endl;
    
    delete [] kmer_count_table;
    return 0;
}


