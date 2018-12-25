//
//  main.cpp
//  BWT
//
//  Created by 聂晨晞 on 2018/12/10.
//  Copyright © 2018 聂晨晞. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <time.h>
#define COMPARE_LENGTH 1000

using namespace std;

void read_genome(string path);
void swap(int &p, int &q);
int Partition(int arrayInput[], int nLow, int nHigh);
void Quick_sort(int arrayInput[], int nLow, int nHigh);
bool compare(int i, int j);
void writeToFile(string path, int BWT[]);

string genome; //为了不把这个基因组传来传去，直接放到全局变量里面
int genome_length;

int main(int argc, const char * argv[]) {
    clock_t start_time, end_time;
    start_time = clock();
    read_genome("/Users/billy/Junior/bioinformatics/Labs/Lab1_CSA/NC_008253.fna");
    // read_genome("/Users/billy/eclipse-workspace/Lab2/test.fa");
    //cout << genome << endl;
    //string a = genome;
    genome_length = int(genome.length());
    cout << genome_length;
    int a = genome_length;
    int *BWT = new int[genome_length];
    for(int i = 0; i < genome_length;i++) {
        BWT[i] = i;
    }
    Quick_sort(BWT, 0, genome_length - 1);
//    cout << endl;
//    for(int i = 0 ;i < genome_length;i++) {
//        //cout << genome[BWT[i] - 1];
//        if(BWT[i] == 0) {
//            cout << "$";
//        }
//        else {
//            cout << genome[BWT[i] - 1];
//        }
//
//    }
    cout << endl;
    writeToFile("/Users/billy/Junior/bioinformatics/Labs/BWT_2/answer.bwt", BWT);
    // writeToFile("/Users/billy/eclipse-workspace/Lab2/test.tally", BWT);
    delete(BWT);
    end_time = clock();
    cout << "运行时间" << (end_time - start_time) / CLOCKS_PER_SEC << endl;
    return 0;
}

void read_genome(string path) {
    string line;
    ifstream infile;
    infile.open(path);
    if(infile.is_open()) {
        getline(infile, line);
        while(getline(infile, line)) {
            //getline(infile, line);
            if(line.back() == '\n') {
                line.pop_back();
            }
            genome = genome + line;
        }
    }
    else{
        cout << "can not find the file to be opened\n";
    }
    genome = genome + "$";
    infile.close();
}

//利用C++的引用交换数组中的两个元素
void swap(int &p, int &q) {
    int temp = p;
    p = q;
    q = temp;
}

//照着这个网址和算法导论里面写的partition函数
//网址：https://www.cnblogs.com/pugang/archive/2012/06/27/2565093.html
int Partition(int arrayInput[], int nLow, int nHigh) {
    int nTemp = arrayInput[nHigh];
    int i = nLow - 1, j = nLow;
    for(;j < nHigh;j++) {
        if(compare(arrayInput[j], nTemp)) {
            i++;
            if(i != j) {
                swap(arrayInput[i], arrayInput[j]);
            }
        }
    }
    swap(arrayInput[i + 1], arrayInput[nHigh]);
    return (i + 1);
}

//Quick sort
void Quick_sort(int arrayInput[], int nLow, int nHigh) {
    if(nLow < nHigh) {
        int nIndex = Partition(arrayInput, nLow, nHigh);
        Quick_sort(arrayInput, nLow, nIndex - 1);
        Quick_sort(arrayInput, nIndex + 1, nHigh);
    }
}

//比较基因组中i和j为开头的两个小片段的大小
//为了计算方便仅仅只比较前1000个字符
//如果以i开头的小片段小于以j开头的小片段，那么返回true，否则返回false
bool compare(int i, int j) {
    for(int k = 0;k < COMPARE_LENGTH;k ++) {
        if(genome[i] == genome[j]) {
            i = (i + 1) % genome_length;
            j = (j + 1) % genome_length;
        } //如果两个字符串在这里是一样的，那么接着比
        else if((genome[i]) == '$' && (genome[j]) != '$') {
            return true;
        }
        else if ((genome[i] != '$') && (genome[j]) == '$') {
            return false;
        }
        else if(genome[i] < genome[j]) {
            return true;
        }
        else if(genome[i] > genome[j]) {
            return false;
        }
    }
    return false;
}

//将结果写入文件中
void writeToFile(string path, int BWT[]) {
    int A = 0, C = 0, G = 0, T = 0;
    ofstream out(path);
    if(out.is_open()) {
        for(int i = 0;i < genome_length;i++) {
            if(BWT[i] == 0) {
                out << "$\t"<< BWT[i] << "\t" << A << "\t" << C << "\t" << G << "\t" << T <<"\n";
            }
            else {
                switch(genome[BWT[i] - 1]){
                    case 'A':
                    case 'a':
                        A++;
                        break;
                    case 'C':
                    case 'c':
                        C++;
                        break;
                    case 'G':
                    case 'g':
                        G++;
                        break;
                    case 'T':
                    case 't':
                        T++;
                        break;
                }

                out << genome[BWT[i] - 1]<< "\t" << BWT[i] << "\t" << A << "\t" << C << "\t" << G << "\t" << T << endl;
            }
        }
    } else {
        cout << "fail to open file to be written";
    }
    out.close();
}
