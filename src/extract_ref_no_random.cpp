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

using namespace std;

const short c = 300;
const char coder_num = 3;
const char least_depth = 3;
const int k = 32;
long long array_size = pow(2, k);
char *kmer_count_table = new char[array_size];
// vector <char>kmer_count_table(array_size);


long slide_window(bool* record_ref_hit, int ref_len, string chr_name, long extract_ref_len, ofstream & interval_file){ //find density hits regions
    int frag_index = 0;
    int window = 500;
    int one_coder_bases = 0; // the number of bases supported by at least one coder in a window
    int three_coder_bases = 0; // the number of bases supported by three coder in a window
    int one_coder_min = window * 0.6;
    int three_coder_min = window * 0;    
    bool good_window;
    bool conti_flag = false;
    int start = 0 ;
    int end = 0;
    int* save_good_intervals = new int[2*ref_len/window];
    short *single_hit_num = new short[ref_len];
    short *trio_hit_num = new short[ref_len];
    // interval_file << "fhsakg"<<endl;
    for (int j = 0; j < ref_len; j++){
        short hit_coder_num = 0;
        for (int p = 0; p < 3; p++){
            if (record_ref_hit[ref_len*p+j] == true){
                // cout <<record_ref_hit[ref_len*p+j]<<endl;
                hit_coder_num += 1;
            }
        }
        // cout << hit_coder_num << endl;
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
        // cout << one_coder_bases << "\t" << three_coder_bases <<endl;
        if (one_coder_bases >= one_coder_min & three_coder_bases >= three_coder_min){
            good_window = true;
        }
        else{
            good_window = false;
        }
        // cout <<three_coder_bases<<"\t"<<hit_coder_num << endl;
        if (conti_flag == false & good_window == true){
            start = j - 2 * window;
            if (start < 1){
                start = 1;
            }
            conti_flag = true;
        }
        if (conti_flag == true & good_window == false){
            end = j + window;
            if (end > ref_len){
                end = ref_len;
            }
            if (frag_index > 0 & start - save_good_intervals[2*frag_index-1] < 200 ){
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
        if (frag_index > 0 & start - save_good_intervals[2*frag_index-1] < 200 ){
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
        interval_file << chr_name << "\t" << save_good_intervals[2*i] << "\t"<< save_good_intervals[2*i+1] << endl;
    }
    // cout << frag_index << endl;
    delete [] save_good_intervals;
    delete [] single_hit_num;
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

void read_ref(string fasta_file, bool* coder, int* base, int k, char* comple, string interval_name)
{
    ifstream fa_file;
    ofstream interval_file;
    fa_file.open(fasta_file, ios::in);
    interval_file.open(interval_name, ios::out | ios::trunc);
    string ref_seq, line_seq;
    ref_seq = "\0";
    int ref_len;
    char n;
    int m;
    int e;
    // int kmer_index, comple_kmer_index, real_index;
    unsigned int kmer_index[3], comple_kmer_index[3];
    unsigned int real_index;
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
    cout <<"Start slide ref..."<<endl;

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
            if (ref_len > 500){   
                for (int i = 0; i < 3; i++){
                    kmer_index[i] = 0; // start kmer
                    comple_kmer_index[i] = 0;
                }
                bool *record_ref_hit = new bool[ref_len*3];
                for (int j = 0; j < ref_len; j++){  
                    support_coder = 0; 
                    n = ref_seq[j];
                    m = n;
                    e = comple[m];  
                    for (int i = 0; i < 3; i++){
                        // covert the base to number
                        // cout << j <<"\t"<< i << "\t" << support_coder<<endl;
                        if (coder[c*i+m] == 5){ //abnormal bases like N, M
                            covert_num = 0;
                        }
                        else{
                            covert_num = coder[c*i+m]; 
                        }
                        if (coder[c*i+e] == 5){
                            comple_num = 0;
                        }
                        else{
                            comple_num = coder[c*i+e];
                        }  
                        // fast version
                        if (j == k-1){
                            convert_ref[50*i+j%k] = covert_num;
                            complemented_ref[50*i+j%k] = comple_num;
                            for (int z = 0; z<k; z++){
                                kmer_index[i] += convert_ref[50*i+z]*base[z];  
                                comple_kmer_index[i] += complemented_ref[50*i+z]*base[k-1-z];  
                            }
                        }
                        else{
                            if (j > k-1){
                                // cout << kmer_index[i]<< "\t" << convert_ref[50*i+j%k] << "\t no"<<base[0] <<endl;
                                kmer_index[i] = (kmer_index[i] - convert_ref[50*i+j%k]*base[0])*2 + covert_num;
                                comple_kmer_index[i] = (comple_kmer_index[i] - complemented_ref[50*i+j%k])/2 + comple_num*base[0];
                            }
                            convert_ref[50*i+j%k] = covert_num;
                            complemented_ref[50*i+j%k] = comple_num;
                        }

                        if (j >= k-1){
                            if (kmer_index[i] > comple_kmer_index[i]){
                                real_index = comple_kmer_index[i];
                            }   
                            else{
                                real_index = kmer_index[i];
                            }
                            // // find hit here, most time-consumption;
                            if (kmer_count_table[real_index] == least_depth){
                                support_coder += 1;
                                record_ref_hit[i*ref_len+j] = true;
                                // cout << "hit" << endl;
                            }
                            else{
                                record_ref_hit[i*ref_len+j] = false;
                                // cout << "no hit!" << endl;
                            }
                        }
                    }
                    // at locus j
                }
                // cout << chr_name<<"\t"<<endl;
                extract_ref_len = slide_window(record_ref_hit, ref_len, chr_name, extract_ref_len, interval_file);
                delete [] record_ref_hit;
            }
            if (ref_index % 1000 == 0){
                time_t t1 = time(0);
                cout << chr_name<<"\t" << ref_index << "\t" <<ref_len << "\t" << extract_ref_len <<"bp\t"<<slide_ref_len<<"bp\t" <<t1-t0<< endl;
            }  
            ref_index += 1;
            ref_seq = "\0";
        }
        else{
            ref_seq += line_seq;
            
        }
             
    }
    cout << "Extracted ref len"<< "\t" << extract_ref_len << endl;
    fa_file.close();
    interval_file.close();
}

// vector<char> 
void read_fastq(string fastq_file, int k, bool* coder, int* base, char* comple)
{
    ifstream fq_file; 
    fq_file.open(fastq_file);
    //fq_file.open("test.fq");
    string reads_seq;
    unsigned int i = 0;
    int converted_reads [450];
    int complemented_reads [450];
    char n;
    int m;
    int e;
    unsigned int kmer_index, comple_kmer_index, real_index, b;   
    int r ;
    short read_len = 150;
    // char abnormal_base[150];
    unsigned seed;
    seed = time(0);
    srand(seed);

    while (fq_file >> reads_seq)
    {
        if (i % 4 == 1){
            if (i % 1000000 == 1){
                cout <<"reads\t" << i <<endl;
            }
            // srand((unsigned)time(NULL));
            r = rand() % 10 ;
            // cout <<r << "r"<<endl;
            if (r < 10){
                for (int j = 0; j < 150; j++){
                    n = reads_seq[j];
                    m = n;
                    e = comple[m];
                    for (int p = 0; p < 3; p++){
                        if (coder[c*p+m] == 5){ //abnormal bases like N, M
                            converted_reads[read_len*p+j] = 0;
                        }
                        else{
                            converted_reads[read_len*p+j] = coder[c*p+m];  
                        }
                        if (coder[c*p+e] == 5){
                            complemented_reads[read_len*p+j] = 0;
                        }
                        else{
                            complemented_reads[read_len*p+j] = coder[c*p+e];
                        }
                    }                 
                }
                for (int p = 0; p < 3; p++){
                    kmer_index = 0;
                    comple_kmer_index = 0;
                    for (int j = 0; j < 150 - k + 1; j++){
                        if (j == 0){
                            for (int z = 0; z<k; z++){
                                kmer_index += converted_reads[read_len*p+j+z]*base[z];  
                                comple_kmer_index += complemented_reads[read_len*p+j+z]*base[k-1-z];
                            }
                        }
                        else{
                            kmer_index = (kmer_index - converted_reads[read_len*p+j-1]*base[0])*2 + converted_reads[read_len*p+j+k-1];
                            comple_kmer_index = (comple_kmer_index - complemented_reads[read_len*p+j-1])/2 + complemented_reads[read_len*p+j+k-1]*base[0];
                        }
                        
                        if (kmer_index > comple_kmer_index){ //use a smaller index
                            real_index = comple_kmer_index;
                        }   
                        else{
                            real_index = kmer_index;
                        }
                        
                        if ((int)kmer_count_table[real_index] < least_depth ){
                            kmer_count_table[real_index] += 1;
                            // cout << (int)kmer_count_table[real_index] << "\t" << endl;
                        }      
                        // else{
                        //     cout << (int)kmer_count_table[real_index] << "\t" << endl;
                        // }
                    }
                }
            }           
        }
        i++;
    }
    fq_file.close();
    // return kmer_count_table;

}

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
    // for (int j = 0; j < k; j++){
    //     for (int i = 0; i < coder_num; i++){
    //         cout << j <<"\t"<< i <<"\t" << choose_coder[j*3+i] << endl;
    //     }
    // }
    return choose_coder;
}

int main( int argc, char *argv[])
{
    // string ref_seq = read_ref();
    bool *coder;
    int *base;
    char *comple;
    short *choose_coder;
    coder = generate_coder(3);
    base = generate_base(k);
    comple = generate_complement();
    choose_coder = random_coder(k); 
    time_t now1 = time(0);

    string fq1 = argv[1];
    string fq2 = argv[2];
    string fasta_file = argv[3];
    string interval_name = argv[4];
    memset(kmer_count_table, 0, sizeof(char)*array_size);
    read_fastq(fq1, k, coder, base, comple);
    read_fastq(fq2, k, coder, base, comple);
    // float empty_rate = cal_tab_empty_rate();
    time_t now2 = time(0);
    cout << "reads finish.\t"<<now2 - now1<<endl;
    read_ref(fasta_file, coder, base, k, comple, interval_name);
        
    time_t now3 = time(0);
    cout << "Finish with time:\t" << now3-now1<<endl;
    delete [] kmer_count_table;

    return 0;
}


