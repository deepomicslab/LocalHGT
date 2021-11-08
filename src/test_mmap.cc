#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <time.h>
#include <cmath> 
#include <iostream>
#include <pthread.h>
#include <thread>
#include <vector>
#include <sys/stat.h> 
#include <sstream>

using namespace std;
// ifstream readin(fname.c_str());
 
long file_size2(string filename)  
{  
    // char filename [] = "/mnt/d/breakpoints/HGT/small/small.1.fq"; 
    struct stat statbuf;  
    stat(filename.c_str(),&statbuf);  
    long size=statbuf.st_size;  
  
    return size;  
}  

void read_fastq(string fastq_file, long start, long end){
    ofstream mid_file;
    string name = to_string(start);
    mid_file.open("/mnt/d/breakpoints/HGT/small/small."+name, ios::out | ios::trunc);

    ifstream fq_file; 
    fq_file.open(fastq_file);
    string reads_seq;
    long pos;

    for (long i = start; i>0; i++){
        fq_file.seekg(i, ios::beg);
        char j;
        fq_file.get(j);
        if (j == '@'){ //only read name has this symbol.
            pos = i;
            // cout <<"hi"<<j<<"\t"<<pos<<endl;
            break;
        }       
    }

    fq_file.seekg(pos, ios::beg);
    long add_size = start;
    while (fq_file >> reads_seq){
        
        if (add_size>=end){
            break;
        }
        mid_file <<start<<"\t"<< reads_seq<< "\t"<<reads_seq.length() <<endl;
        add_size += reads_seq.length();
    }
    fq_file.close();
    mid_file.close();
}


int main(){

    string index_name = "/mnt/d/breakpoints/HGT/test/ref.index.dat";
    string fai_name = "/mnt/d/breakpoints/HGT/UHGG/UHGG_reference.formate.fna.fai";
    ifstream index_file;
    ifstream fai_file;
    fai_file.open(fai_name, ios::in); 
    string line;
    string aa;
    int ref_len;


    int thread_num = 4;
    long index_size = file_size2(index_name);
    long each_index_size = index_size/thread_num + 1;
    cout <<index_size<<endl;
    cout <<each_index_size<<endl;

    unsigned int real_index = 0;

    index_file.open(index_name, ios::in | ios::binary); 
    long pos = 100 * 4;
    index_file.seekg(pos, ios::beg);
    long start, end;
    start = pos;
    int i = 0;
    int k = 32;
    int coder_num =3;
    long add;

    long a = 0;
    cout <<start <<endl;
    while(!fai_file.eof()){
        getline(fai_file,line);
        std::istringstream iss(line);
        if (!(iss >> aa >> ref_len)) { break; }
        a += ref_len;
        // cout <<ref_len<< "\t"<< a<<endl;
        index_file.read(reinterpret_cast<char*>(&real_index), sizeof(unsigned int));
        add = 4*((ref_len-k+1)*coder_num+1); //the size of the genome.
        if (pos -start > each_index_size){
            end = pos + add;
            cout << start <<"hh\t"<<end<<"\t"<<ref_len <<endl;
            // break;
            start = end;
        
        }
        pos += add;

        i += 1;
    }

    cout << start <<"hh\t"<<index_size<<"\t"<<ref_len <<endl;
    cout <<i<<"\tend\t"<<a<<endl;

    fai_file.close();
    index_file.close();
    return 0;
}

// int main(){

//     string index_name = "/mnt/d/breakpoints/HGT/test/ref.index.dat";
//     int thread_num = 10;
//     long index_size = file_size2(index_name);
//     long each_index_size = index_size/thread_num + 1;
//     cout <<index_size<<endl;
//     cout <<each_index_size<<endl;
//     ifstream index_file;
//     unsigned int real_index = 0;

//     index_file.open(index_name, ios::in | ios::binary); 
//     long pos = 100 * 4;
//     index_file.seekg(pos, ios::beg);
//     long start, end;
//     start = pos;
//     int i = 0;
//     int k = 32;
//     int coder_num =3;
//     long add;

//     // while (!index_file.eof())
//     while (pos < index_size ){
//         index_file.read(reinterpret_cast<char*>(&real_index), sizeof(unsigned int));
//         add = 4*((real_index-k+1)*coder_num+1); //the size of the genome.

//         if (pos -start > each_index_size){
//             end = pos + add;
//             cout << start <<"\t"<<end<<"\t"<<real_index <<endl;
//             // break;
//             start = end;

//         }
//         // cout<<pos<<"\tref_len\t"<<real_index<<endl;
//         // (ref_len-k+1)*coder_num;
//         pos += add;
//         // cout <<pos<<endl;
//         index_file.seekg(pos, ios::beg);
//         real_index = 0;
//         if (i % 1000000 == 0){
//             // break;
//             cout<<i<<endl;
//         }
        
//         i += 1;
        
//     }
//     index_file.close();

//     return 0;
// }




// int main(){
//     int thread_num = 2;
//     long start = 0;
//     long end = 0;

//     string fastq_file = "/mnt/d/breakpoints/HGT/small/small.1.fq";
//     long size = file_size2(fastq_file);
//     cout <<size<<endl;

//     long each_size = size/thread_num;

//     std::vector<std::thread>threads;
//     for (int i=0; i<thread_num; i++){
//         start = i*each_size;
//         end = (i+1)*each_size;
//         if (i == thread_num-1){
//             end = size;
//         }
//         cout <<start<<"\t"<<end<<endl;

//         threads.push_back(thread(read_fastq, fastq_file, start, end));
//     }

// 	for (auto&th : threads)
// 		th.join();

// }




/*

int cal_thread_line(long line_num, int thread_num){
    return line_num/thread_num;
}

void read_fastq(string fastq_file, long start, long end){
    ofstream mid_file;
    string name = to_string(start);
    mid_file.open("/mnt/d/breakpoints/HGT/small/small."+name, ios::out | ios::trunc);

    ifstream fq_file; 
    fq_file.open(fastq_file);
    string reads_seq;
    int pos = start * 33;
    fq_file.seekg(pos, ios::beg);
    long line_index = start;
    while (fq_file >> reads_seq){
        
        if (line_index>=end){
            break;
        }
        mid_file <<start<<"\t"<< reads_seq<< "\t"<<reads_seq.length() <<endl;
        

        line_index += 1;
    }
    fq_file.close();
    mid_file.close();
}


int main(){
    int thread_num = 2;
    long line_num = 20;
    long each_thread_line_num = cal_thread_line(line_num, thread_num);
    long start = 0;
    long end = 0;

    string fastq_file = "/mnt/d/breakpoints/HGT/small/small.1.fq";
    long size = file_size2(fastq_file);
    cout <<size<<endl;
    std::vector<std::thread>threads;
    for (int i=0; i<thread_num; i++){
        start = i*each_thread_line_num;
        end = (i+1)*each_thread_line_num;
        if (i == thread_num-1){
            end = line_num;
        }
        cout <<start<<"\t"<<end<<endl;

        threads.push_back(thread(read_fastq, fastq_file, start, end));
    }

	for (auto&th : threads)
		th.join();

}






int* times(int h){
    for (int i=0; i<1000000000; i++){
        // int index = random(); 
        int  j = i;
        // if (index % 3 == h){
            // kmer_count_table[index] = kmer_count_table[index] +1 -1;
        // }
        
    }    
}

int main (){

    // int nums[1000000];

    // for (int i=0; i<1000000; i++){
    //     nums[i] = i;
    // }
    time_t time0 = time(0);

    times(0);
    times(1);
    times(2);
    time_t time1 = time(0);
    cout <<time1-time0<<endl;   

    thread t1 (times, 0);
    thread t2 (times, 1);
    thread t3 (times, 2);
    t1.join();
    t2.join();
    t3.join();

    time_t time2 = time(0);
    cout <<time2-time1<<endl;

    // for (int i=0; i<100000000; i++){
    //     cout << nums[i] <<endl ;
    // }

    
}
*/