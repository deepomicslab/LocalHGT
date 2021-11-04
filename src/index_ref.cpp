#include <iostream>
#include <fstream>
#include <cmath> 
#include <cstdlib> 
#include <string>
#include <vector>
#include <time.h>
#include <cstdlib>
#include <ctime>

using namespace std;

short c = 300;
char coder_num = 3;
char least_depth = 3;


void read_ref(string fasta_file, bool* coder, int* base, int k, char* comple, string index_file_name)
{
    ifstream fa_file;
    ofstream index_file;
    fa_file.open(fasta_file, ios::in);
    index_file.open(index_file_name, ios::out | ios::binary);
    string ref_seq, line_seq;
    ref_seq = "\0";
    int ref_len;
    int m, e, kmer_index, comple_kmer_index, real_index;
    int ref_index = 0;
    long extract_ref_len = 0;
    string chr_name ;
    string pre_name = "start";
    time_t t0 = time(0);
    char covert_num, comple_num, n;
    char convert_ref[50];
    char complemented_ref[50];

    while (fa_file >> line_seq){ //remember the last chr
        if (line_seq[0] == '>'){
            chr_name = pre_name;
            pre_name = line_seq.substr(1);
            if (ref_index > 1005){
                break;
            }
            ref_len= ref_seq.length();
            // cout << chr_name <<"\t"<< ref_index <<"\t" << ref_len << endl;
            if (ref_len > 500){   
                bool *record_ref_hit = new bool[ref_len*3];
                extract_ref_len += ref_len;
                for (int p = 0; p < 3; p++){
                    kmer_index = 0; // start kmer
                    comple_kmer_index = 0;
                    for (int j = 0; j < ref_len - k + 1; j++){   
                        // covert the base to number
                        n = ref_seq[j];
                        m = n;
                        e = comple[m];   
                        if (coder[c*p+m] == 5){ //abnormal bases like N, M
                            covert_num = 0;
                        }
                        else{
                            covert_num = coder[c*p+m]; 
                        }
                        if (coder[c*p+e] == 5){
                            comple_num = 0;
                        }
                        else{
                            comple_num = coder[c*p+e];
                        }  
                        if (j == k-1){
                            convert_ref[j%k] = covert_num;
                            complemented_ref[j%k] = comple_num;
                            for (int z = 0; z<k; z++){
                                kmer_index += convert_ref[z]*base[z];  
                                comple_kmer_index += complemented_ref[z]*base[k-1-z];  
                            }
                        }
                        else{
                            if (j > k-1){
                                kmer_index = (kmer_index - convert_ref[j%k]*base[0])*2 + covert_num;
                                comple_kmer_index = (comple_kmer_index - complemented_ref[j%k])/2 + comple_num*base[0];
                            }
                            convert_ref[j%k] = covert_num;
                            complemented_ref[j%k] = comple_num;
                        }


                        if (kmer_index > comple_kmer_index){
                            real_index = comple_kmer_index;
                        }   
                        else{
                            real_index = kmer_index;
                        }
                        // save index here
                        // index_file.write(real_index, 4);
                        // index_file.write((char *)(&real_index), sizeof(real_index));
                        // cout << real_index<<"\t"<<endl;
                    }
                }
                // cout << chr_name<<"\t"<<endl;
                if (ref_index % 200 == 0){
                    time_t t1 = time(0);
                    cout << chr_name<<"\t" << ref_index << "\t" <<ref_len << "\t" << extract_ref_len <<"bp\t"<<t1-t0<< endl;
                }   
                delete [] record_ref_hit;        
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
    index_file.close();
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
    static int base [30];
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
}

short * random_coder(int k)
{
    static short permu[18] = {0,1,2,0,2,1,1,2,0,1,0,2,2,0,1,2,1,0};
    static short choose_coder[90];
    int r;
    unsigned seed;
    seed = time(0);
    srand(seed);
    for (int j = 0; j < k; j++){
        r = rand() % 6;
        cout << r << endl;
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
    int k = 30;
    coder = generate_coder(3);
    base = generate_base(k);
    comple = generate_complement();
    // choose_coder = random_coder(k); 
    time_t now1 = time(0);
    // string fasta_file = "/mnt/d/breakpoints/HGT/UHGG/UHGG_reference.formate.fna";
    // string fq1 = "/mnt/d/breakpoints/HGT/small/small.1.fq";
    // string fq2 = "/mnt/d/breakpoints/HGT/small/small.2.fq";
    // string fq1 = "/mnt/d/breakpoints/HGT/cami_fq/species10_snp0.08_depth20_reads150_sample_0_high_HGT.1.fq";
    // string fq2 = "/mnt/d/breakpoints/HGT/cami_fq/species10_snp0.08_depth20_reads150_sample_0_high_HGT.2.fq";
    // string interval_name = "test.interval.txt";
    string fasta_file = argv[1];
    string index_file_name = "/mnt/d/breakpoints/HGT/test1/ref_index.dat";
    
    read_ref(fasta_file, coder, base, k, comple, index_file_name);
        
    time_t now3 = time(0);
    cout << "Finish with time:\t" << now3-now1<<endl;
    return 0;
}