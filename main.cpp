#include <iostream>
#include <fstream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <time.h>
#include <string>
#include <algorithm>
#include <cassert>

//****************** BAM FLAG*********************************************
/*! @abstract the read is paired in sequencing, no matter whether it is mapped in a pair */
#define BAM_FPAIRED        1
/*! @abstract the read is mapped in a proper pair */
#define BAM_FPROPER_PAIR   2
/*! @abstract the read itself is unmapped; conflictive with BAM_FPROPER_PAIR */
#define BAM_FUNMAP         4
/*! @abstract the mate is unmapped */
#define BAM_FMUNMAP        8
/*! @abstract the read is mapped to the reverse strand */
#define BAM_FREVERSE      16
/*! @abstract the mate is mapped to the reverse strand */
#define BAM_FMREVERSE     32
/*! @abstract this is read1 */
#define BAM_FREAD1        64
/*! @abstract this is read2 */
#define BAM_FREAD2       128
/*! @abstract not primary alignment */
#define BAM_FSECONDARY   256



using namespace std;

void alignment(string seq, string quality, vector<int>& result, int where_read_is_from);
int * findline(char lastbase1, char base, int baseF, int baseL);
int findone(char base, int baseF);
void alignOne(string seq,vector<int>& result);
void split(string target, vector<string>& strs);
string read_reverse(string read);
int flag(int read);

int A1;
int A2;
int C1;
int C2;
int G1;
int G2;
int T1;
int T2; //为了让alignment的参数列表不要太长，直接改成全局变量
int A, C, G, T, ref_length;

int baseF_vector[20];
int baseL_vector[20]; //由于AGCT对应的ASCII码值减去65之后分别为0，2，6，19，开20个数组方便后面的转换

vector<string> tally;//没有必要每次循环都读取tally，直接放全局变量就好

//************* FOR BAM FLAG ********************
bool paired = true;
bool paired2 = true;
bool proper_pair = true;
bool proper_pair2 = true;
bool unmap = true;
bool unmap2 = true;
bool mun_map = true;
bool mun_map2 = true;
bool reverse1 = false;
bool reverse2 = false;
bool mreverse = false;
bool mreverse2 = false;

//****************** FOR SAM CIGAR后面那个*******************
bool exact_match = true;


int main() {
    clock_t start_time, end_time;
    start_time = clock();
    string tally_bdx_path = "/Users/billy/Junior/bioinformatics/Labs/Lab2_Bowtie/test_files/NC_008253.tally.bdx";
    // string tally_bdx_path = "/Users/billy/Junior/bioinformatics/Labs/Lab2_Bowtie/test_files/test.tally.bdx";
    ifstream bdx(tally_bdx_path);

    string ref_name;
    if(bdx.is_open()) {
        bdx >> ref_name;
        string temp;
        bdx >> temp; bdx>> temp;
        A = stoi(temp);
        //cout << A;
        bdx >> temp; bdx>> temp;
        C = stoi(temp);
        //cout << C;
        bdx >> temp; bdx>> temp;
        G = stoi(temp);
        //cout << G;
        bdx >> temp; bdx>> temp;
        T = stoi(temp);
        //cout << T;
        ref_length = A + C + G + T;
    } else {
        cout << "could not open tally index file";
    }

    bdx.close();

    //定义AGCT的边界
    A1 = 1;
    A2 = A;
    C1 = A + 1;
    C2 = A + C;
    G1 = C2 + 1;
    G2 = G1 + G - 1;
    T1 = G2 + 1;
    T2 = T1 + T - 1;

    ofstream output_SAM("/Users/billy/Junior/bioinformatics/Labs/Lab2_Bowtie/test_files/NC_008253.sam");
    if(output_SAM.is_open()) {
        output_SAM << "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:" << ref_name << "\tLN:" << ref_length <<"\t\n@PG\tID:myBWT_Aligner\tPN:myBWT_Aligner\tVN:1.0\n";
    } else {
        cout << "could not open sam file to write";
    }


    ifstream read_file("/Users/billy/Junior/bioinformatics/Labs/Lab2_Bowtie/test_files/Ecoli_4x.fq1");
    ifstream read_file2("/Users/billy/Junior/bioinformatics/Labs/Lab2_Bowtie/test_files/Ecoli_4x.fq2");
    int i = 0;

    string tally_path = "/Users/billy/Junior/bioinformatics/Labs/Lab2_Bowtie/test_files/NC_008253.tally";
    ifstream tally_file(tally_path);
    while(!tally_file.eof()) {
        string temp;
        getline(tally_file, temp);
        tally.push_back(temp);
    }
    tally_file.close();//把tally里面的东西读到内存里面的一个可变数组里面

    string read_head;
    string read;
    string quality;
    string plus;
    i = 0;

    string read_head2;
    string read2;
    string quality2;
    string plus2;
    // int u = 1;
    if(read_file.is_open()) {
        while(!read_file.eof() && !read_file2.eof()) {
            cout << "正在处理第" << i << "对序列" <<endl;
            i++;
            getline(read_file, read_head);
            getline(read_file, read);
            getline(read_file, plus);
            getline(read_file, quality);

            getline(read_file2, read_head2);
            getline(read_file2, read2);
            getline(read_file2, plus2);
            getline(read_file2, quality2);

            vector<int> posi;
            vector<int> posi2;
            alignment(read, quality, posi, 1);
            alignment(read2, quality2, posi2, 2);

            //检查有没有比对上，如果这个正链没有比对上，那么试试反链
            bool posi_flag_empty = true;
            for(int i = 0;i < posi.size(); i++) {
                if(posi[i] > 0 && posi[i] <= ref_length) {
                    posi_flag_empty = false;
                }
            }

            if(posi_flag_empty) {
                read = read_reverse(read);
                posi.clear();
                alignment(read, quality, posi, 1);
            }

            for(int i = 0;i < posi.size();i++) {
                if(posi[i] > 0 && posi[i] <= ref_length) {
                    unmap = false;
                    mun_map2 = false;
                }
            }

            if(posi_flag_empty && !unmap) {
                reverse1 = true;
                mreverse2 = true;
            }

            bool posi_flag_empty2 = true;
            for(int i = 0;i < posi2.size();i++) {
                if(posi2[i] > 0 && posi2[i] <= ref_length) {
                    posi_flag_empty2 = false;
                }
            }

            if(posi_flag_empty2) {
                read2 = read_reverse(read2);
                posi2.clear();
                alignment(read2, quality2, posi2, 2);
            }


            for(int i = 0;i < posi.size();i++) {
                if(posi[i] > 0 && posi[i] <= ref_length) {
                    unmap2 = false;
                    mun_map = false;
                }
            }


            if (posi_flag_empty2) {
                reverse2 = true;
                mreverse = true;
            }


            for(int j = 0;j < posi.size(); j++) {
                int len = read.length();
                if(read_head.length() > 1) {
                    read_head = read_head.substr(1, read_head.length() - 1);
                }
                if(posi[j] > 0 && posi[j] < ref_length) {
                    if(exact_match){
                        output_SAM << read_head << "\t"<< flag(1) <<"\t" << ref_name << "\t" << posi[j] << "\t42\t" << len << "M\t=\t0\t0\t"
                                   << read << "\t" << quality << "\n";

                    } else {
                        output_SAM << read_head << "\t"<< flag(1) <<"\t" << ref_name << "\t" << posi[j] << "\t42\t" << len << "M\t*\t0\t0\t"
                                   << read << "\t" << quality << "\n";

                    }

                }

                if(posi2[j] > 0 && posi2[j] < ref_length) {
                    len = read2.length();
                    if(read_head2.length() > 1) {
                        read_head2 = read_head2.substr(1, read_head2.length() - 1);
                    }
                    if(exact_match) {
                        output_SAM << read_head2 << "\t" << flag(2) << "\t" << ref_name << "\t" << posi2[j] << "\t42\t" << len << "M\t=\t0\t0\t"
                                   << read2 << "\t" << quality2 << "\n";
                    } else {
                        output_SAM << read_head2 << "\t" << flag(2) << "\t" << ref_name << "\t" << posi2[j] << "\t42\t" << len << "M\t*\t0\t0\t"
                                   << read2 << "\t" << quality2 << "\n";
                    }

                }
            }
        }
    }
    read_file.close();

    output_SAM.close();
    end_time = clock();
    cout << "运行时间：" << (double)(end_time - start_time) / CLOCKS_PER_SEC << "秒" << endl;
    return 0;
}

//将read align到基因组对应的位置上去
void alignment(string seq, string quality, vector<int>& result, int where_read_is_from) {
    int count = 1;
    int seqlen = seq.length();
    char base = seq[seqlen - count];
    int check1, check2; //两个范围
    int base1, base2;

    if(base == 'A') {
        base1 = A1;
        base2 = A2;
    } else if (base == 'C') {
        base1 = C1;
        base2 = C2;
    } else if (base == 'G') {
        base1 = G1;
        base2 = G2;
    } else if (base == 'T') {
        base1 = T1;
        base2 = T2;
    } else {
        cout << "Base in alignment on line 150: characters other than A,G,C,T are not supported!" << endl;
        exit(-1);
    }
    check1 = base1;
    check2 = base2;

    vector<string> strs;
    // boost::split(strs, tally[base1], boost::is_any_of("\t"));
    split(tally[base1], strs);
    char lastbase1;

    int posi;
    lastbase1 = strs[0][0]; posi = stoi(strs[1]); baseF_vector[0] = stoi(strs[2]);
    baseF_vector[2] = stoi(strs[3]); baseF_vector[6] = stoi(strs[4]); baseF_vector[19] = stoi(strs[5]);

    strs.clear();
    // boost::split(strs, tally[base2], boost::is_any_of("\t"));
    split(tally[base2], strs);


    posi = stoi(strs[1]); baseL_vector[0] = stoi(strs[2]);
    baseL_vector[2] = stoi(strs[3]); baseL_vector[6] = stoi(strs[4]); baseL_vector[19] = stoi(strs[5]);

    base = seq[seqlen - count - 1];
    int baseF, baseL;
    int z = (int)base - 65;
    baseF = baseF_vector[z];
    baseL = baseL_vector[z];

    int fline, lline; int *findline_result;
    findline_result = findline(lastbase1, base, baseF, baseL);
    fline = *findline_result;
    lline = *(findline_result + 1); //C++函数返回数组
    count++;

    while(count < seqlen) {
        strs.clear();
        base = seq[seqlen - count - 1];
        // boost::split(strs, tally[fline], boost::is_any_of("\t"));
        split(tally[fline], strs);
        lastbase1 = strs[0][0];
        posi = stoi(strs[1]);

        baseF_vector[0] = stoi(strs[2]); baseF_vector[2] = stoi(strs[3]);
        baseF_vector[6] = stoi(strs[4]); baseF_vector[19] = stoi(strs[5]);


        if (fline < lline) {
            strs.clear();
            // boost::split(strs, tally[lline], boost::is_any_of("\t"));
            split(tally[lline], strs);
            posi = stoi(strs[1]);

            baseL_vector[0] = stoi(strs[2]);
            baseL_vector[2] = stoi(strs[3]);
            baseL_vector[6] = stoi(strs[4]);
            baseL_vector[19] = stoi(strs[5]);

            z = (int)base - 65;
            baseF = baseF_vector[z];
            baseL = baseL_vector[z];

            findline_result = findline(lastbase1, base, baseF, baseL);
            fline = *findline_result;
            lline = *(findline_result + 1);

            count++;
        } else {
            if (base == lastbase1) {
                z = (int)base - 65;
                baseF = baseF_vector[z];
                fline = findone(base, baseF);
                lline = fline;
                count++;
            } else {
                fline = 0;
                lline = 0;
                break;
            }
        }
    }

    if(fline != 0) {
        int all = lline - fline + 1;
        while(fline <= lline) {
            strs.clear();
            // boost::split(strs, tally[fline], boost::is_any_of("\t"));
            split(tally[fline], strs);
            posi = stoi(strs[1]);
            posi++;
            result.push_back(posi);
            fline++;
        }
    } else {
        if(where_read_is_from == 1)
            proper_pair = false;
        else if(where_read_is_from == 2)
            proper_pair2 = false;
        int start_check_point = seqlen - count - 1;
        int check_point = start_check_point;
        if(check_point < 0) {
            result.push_back(0);
            return;
        }

        exact_match = false;
        string mut = seq.substr(check_point, 1);

        string seq1,seq2,seq3;
        if (mut == "A") {
            seq[check_point] = 'C';
            seq1 = seq;
            seq[check_point] = 'G';
            seq2 = seq;
            seq[check_point] = 'T';
            seq3 = seq;
        } else if (mut == "C") {
            seq[check_point] = 'A';
            seq1 = seq;
            seq[check_point] = 'G';
            seq2 = seq;
            seq[check_point] = 'T';
            seq3 = seq;
        } else if (mut == "G") {
            seq[check_point] = 'A';
            seq1 = seq;
            seq[check_point] = 'C';
            seq2 = seq;
            seq[check_point] = 'T';
            seq3 = seq;
        } else {
            seq[check_point] = 'A';
            seq1 = seq;
            seq[check_point] = 'C';
            seq2 = seq;
            seq[check_point] = 'G';
            seq3 = seq;
        }

        vector<int> result1;
        alignOne(seq1, result1);
        vector<int> result2;
        alignOne(seq2, result2);
        vector<int> result3;
        alignOne(seq3, result3);

        result.insert(result.end(), result1.begin(), result1.end());
        result.insert(result.end(), result2.begin(), result2.end());
        result.insert(result.end(), result3.begin(), result3.end());
    }
}


int * findline(char lastbase1, char base, int baseF, int baseL) {
    if(lastbase1 == base) {
        baseF--;
    }
    static int result[2]; //result数组的0表示fline，1表示lline
    if(baseF != baseL) {
        if(base == 'A') {
            result[0] = baseF + 1;
            result[1] = baseL;
        } else if(base == 'C') {
            result[0] = A2 + baseF + 1;
            result[1] = A2 + baseL;
        } else if(base == 'G') {
            result[0] = C2 + baseF + 1;
            result[1] = C2 + baseL;
        } else if(base == 'T') {
            result[0] = G2 + baseF + 1;
            result[1] = G2 + baseL;
        }
        return result;
    } else {
        result[0] = 0;
        result[1] = 0;
        return result;
    }
}

int findone(char base, int baseF) {
    int fline;
    if(base == 'A') {
        fline = baseF;
    } else if(base == 'C') {
        fline = A + baseF;
    } else if(base == 'G') {
        fline = A + C + baseF;
    } else if (base == 'T') {
        fline = A + C + G + baseF;
    }
    return fline;
}

void alignOne(string seq,vector<int>& result) {
    int count = 1;
    int seqlen = seq.length();
    char base = seq[seqlen - count];

    vector<string> strs;
    int base1, base2;
    if(base == 'A') {
        base1 = A1;
        base2 = A2;
    } else if (base == 'C') {
        base1 = C1;
        base2 = C2;
    } else if (base == 'G') {
        base1 = G1;
        base2 = G2;
    } else if (base == 'T') {
        base1 = T1;
        base2 = T2;
    } else {
        cout << "Base in alignment on line 331: characters other than A,G,C,T are not supported!" << endl;
        exit(-1);
    }

    strs.clear();
    // boost::split(strs, tally[base1], boost::is_any_of("\t"));
    split(tally[base1], strs);
    char lastbase1 = strs[0][0];
    int posi;
    posi = stoi(strs[1]); baseF_vector[0] = stoi(strs[2]); baseF_vector[2] = stoi(strs[3]);
    baseF_vector[6] = stoi(strs[4]);  baseF_vector[19] = stoi(strs[5]);

    strs.clear();
    // boost::split(strs, tally[base2], boost::is_any_of("\t"));
    split(tally[base2], strs);

    posi = stoi(strs[1]); baseL_vector[0] = stoi(strs[2]); baseL_vector[2] = stoi(strs[3]);
    baseL_vector[6] = stoi(strs[4]); baseL_vector[19] = stoi(strs[5]);

    base = seq[seqlen - count - 1];
    int baseF, baseL;

    int z = (int)base - 65;
    baseF = baseF_vector[z];
    baseL = baseL_vector[z];

    int *findline_result = findline(lastbase1, base, baseF, baseL);
    int fline, lline;
    fline = *findline_result;
    lline = *(findline_result + 1);
    count ++;

    while(count < seqlen) {
        base = seq[seqlen - count - 1];
            strs.clear();
            // boost::split(strs, tally[fline], boost::is_any_of("\t"));
            split(tally[fline], strs);
            lastbase1 = strs[0][0];

            posi = stoi(strs[1]); baseF_vector[0] = stoi(strs[2]); baseF_vector[2] = stoi(strs[3]);
            baseF_vector[6] = stoi(strs[4]); baseF_vector[19] = stoi(strs[5]);

            if(fline < lline) {
                strs.clear();
                // boost::split(strs, tally[lline], boost::is_any_of("\t"));
                split(tally[lline], strs);

                posi = stoi(strs[1]); baseL_vector[0] = stoi(strs[2]); baseL_vector[2] = stoi(strs[3]);
                baseL_vector[6] = stoi(strs[4]); baseL_vector[19] = stoi(strs[5]);

                int z = (int)base - 65;
                baseF = baseF_vector[z];
                baseL = baseL_vector[z];

                findline_result = findline(lastbase1, base, baseF, baseL);
                fline = *findline_result;
                lline = *(findline_result + 1);
                count++;
        } else {
            if(base == lastbase1) {
                baseF = baseF_vector[(int)base - 65];
                fline = findone(base, baseF);
                lline = fline;
                count ++;
            } else {
                fline = 0;
                lline = 0;
                break;
            }
        }
    }

    if(fline != 0) {
        int all = lline - fline + 1;
        while(fline <= lline) {
            strs.clear();
            // boost::split(strs, tally[fline], boost::is_any_of("\t"));
            split(tally[fline], strs);
            posi = stoi(strs[1]);
            posi++;
            result.push_back(posi);
            fline++;
        }
    } else {
        result.push_back(0);
    }
}

void split(string target, vector<string>& strs) {
    //string sentence = "C\t4938920\t0\t1\t0\t0";

    istringstream iss(target);

    copy(istream_iterator<string>(iss),
         istream_iterator<string>(),
         back_inserter(strs));
}

string read_reverse(string read) {
    reverse(read.begin(), read.end());
    for(int i = 0; i < read.length(); i++) {
        switch (read[i]) {
            case 'A':
                read[i] = 'T';
                break;
            case 'C':
                read[i] = 'G';
                break;
            case 'G':
                read[i] = 'C';
                break;
            case 'T':
                read[i] = 'A';
                break;
            default:
                cout << "只有AGCT被支持!" << endl;
        }
    }
    return read;
}

int flag(int read) {

    int flag = 0;
    if(read == 1) {
        if(paired) {
            flag += BAM_FPAIRED;
        }
        if(proper_pair) {
            flag += BAM_FPROPER_PAIR;
        }
        if(unmap) {
            flag += BAM_FUNMAP;
        }
        if(mun_map) {
            flag += BAM_FMUNMAP;
        }
        if(reverse1) {
            flag += BAM_FREVERSE;
        }
        if(mreverse) {
            flag += BAM_FMREVERSE;
        }
        flag += BAM_FREAD1;
    }
    else if (read == 2) {
        if(paired2) {
            flag += BAM_FPAIRED;
        }
        if(proper_pair2) {
            flag += BAM_FPROPER_PAIR;
        }
        if(unmap2) {
            flag += BAM_FUNMAP;
        }
        if(mun_map2) {
            flag += BAM_FMUNMAP;
        }
        if(reverse2) {
            flag += BAM_FREVERSE;
        }
        if(mreverse2) {
            flag += BAM_FMREVERSE;
        }
        flag += BAM_FREAD2;
    }
    return flag;
}