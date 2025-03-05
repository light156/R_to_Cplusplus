#include "trans_ancestry_cojo.h"


void taCOJO::read_bimfile(string bimFile, bool isFirst) 
{
    // Read bim file: recombination rate is defined between SNP i and SNP i-1
    int ibuf = 0;
    string A1_buf = "0", A2_buf = "0";
    double dbuf = 0.0;
    string snp_buf;
    ifstream Bim(bimFile.c_str());

    if (!Bim) LOGGER.e(0, "cannot open the file [" + bimFile + "] to read");
    LOGGER << "Reading PLINK BIM file from [" + bimFile + "]..." << endl;

    int snp_num = 0;
    string A1, A2;
    unordered_set<string> commonSNP_buf;
    unordered_map<string, ItemBim>::iterator item_iter;

    while (Bim) {
        Bim >> ibuf;
        if (Bim.eof()) break;
        Bim >> snp_buf;
        Bim >> dbuf;
        Bim >> ibuf;
        Bim >> A1_buf;
        StrFunc::to_upper(A1_buf);
        Bim >> A2_buf;
        StrFunc::to_upper(A2_buf);

        if (isFirst) {
            bimSNP1.push_back(snp_buf);

            commonSNP_buf.insert(snp_buf);
            if (snp_num == commonSNP_buf.size())
                LOGGER.e(0, "Duplicate SNP " + snp_buf + " found in the file [" + bimFile + "]");

            ItemBim newItem = ItemBim(A1_buf, A2_buf, false);
            bimData.insert(pair<string, ItemBim>(snp_buf, newItem));
            snp_num++;
        } else {
            bimSNP2.push_back(snp_buf);
            
            if((item_iter = bimData.find(snp_buf)) != bimData.end()) {
                A1 = item_iter->second.A1;
                A2 = item_iter->second.A2;

                if ((A1_buf == A1 && A2_buf == A2) || (A1_buf == A2 && A2_buf == A1)) {
                    commonSNP_buf.insert(snp_buf);
                    if (snp_num == commonSNP_buf.size())
                        LOGGER.e(0, "Duplicate SNP " + snp_buf + " found in the file [" + bimFile + "]");
                    if (A1_buf == A2)
                        item_iter->second.swap = true;

                    snp_num++;
                } 
            }
        }
    }
    Bim.close();

    commonSNP.swap(commonSNP_buf);
    unordered_set<string>().swap(commonSNP_buf);

    LOGGER << commonSNP.size() << " SNPs included" << endl;
}


void taCOJO::read_cojofile(string cojoFile, unordered_map<string, ItemCojo> &cojoData) 
{    
    LOGGER << "Reading GWAS summary-level statistics from [" + cojoFile + "] ..." << endl;
    ifstream Meta(cojoFile.c_str());
    if (!Meta) LOGGER.e(0, "cannot open the file [" + cojoFile + "] to read");

    // SNP A1 A2 freq b se p N 
    double freq_buf = 0.0, b_buf = 0.0, se_buf = 0.0, p_buf = 0.0, N_buf = 0.0;
    string A1_buf, A2_buf, snp_buf, str_buf0, str_buf;

    vector<string> vs_buf;
    getline(Meta, str_buf); // the header line
    if (StrFunc::split_string(str_buf, vs_buf) < 7) LOGGER.e(0, "format error in the input file [" + cojoFile + "]");

    string A1, A2;
    unordered_set<string> commonSNP_buf;
    unordered_map<string, ItemBim>::iterator item_iter;

    while (Meta) {
        getline(Meta, str_buf0);
        stringstream iss(str_buf0);
        iss >> snp_buf >> A1_buf >> A2_buf;
        
        StrFunc::to_upper(A1_buf);
        StrFunc::to_upper(A2_buf);
        iss >> str_buf;
        freq_buf = atof(str_buf.c_str());
        if (freq_buf<0.01 || freq_buf>0.99) continue;
        iss >> str_buf;
        if (str_buf == "NA" || str_buf == ".") continue;
        b_buf = atof(str_buf.c_str());
        iss >> str_buf;
        if (str_buf == "NA" || str_buf == "." || str_buf == "0") continue;
        se_buf = atof(str_buf.c_str());
        iss >> str_buf;
        if (str_buf == "NA" || str_buf == ".") continue;
        p_buf = atof(str_buf.c_str());
        iss >> str_buf;
        if (str_buf == "NA" || str_buf == ".") continue;
        N_buf = atof(str_buf.c_str());
        if (N_buf < 10) LOGGER.e(0, "invalid sample size in line:\n\"" + str_buf0 + "\"");
        if (Meta.eof()) break;

        if((item_iter = bimData.find(snp_buf)) != bimData.end()) {
            A1 = item_iter->second.A1;
            A2 = item_iter->second.A2;

            if (A1_buf == A1 && A2_buf == A2) {
                commonSNP_buf.insert(snp_buf);
                ItemCojo newItem = ItemCojo(freq_buf, b_buf, se_buf, p_buf, N_buf);
                cojoData.insert(pair<string, ItemCojo>(snp_buf, newItem));
            } else if (A1_buf == A2 && A2_buf == A1) {
                commonSNP_buf.insert(snp_buf);
                ItemCojo newItem = ItemCojo(1-freq_buf, -b_buf, se_buf, p_buf, N_buf);
                cojoData.insert(pair<string, ItemCojo>(snp_buf, newItem));
            } 
        }
    }
    Meta.close();
    
    commonSNP.swap(commonSNP_buf);
    unordered_set<string>().swap(commonSNP_buf);
    LOGGER << commonSNP.size() << " SNPs included" << endl;
}


void taCOJO::read_famfile(string famFile, int &indi_num) 
{
    ifstream Fam(famFile.c_str());
    if (!Fam) LOGGER.e(0, "cannot open the file [" + famFile + "] to read");
    LOGGER << "Reading PLINK FAM file from [" + famFile + "]" << endl;

    indi_num = 0;
    string str_buf1, str_buf2, indi_id_buf;
    unordered_set<string> indi_ids;

    while (Fam) {
        Fam >> str_buf1;
        if (Fam.eof()) break;
        Fam >> str_buf2;

        indi_ids.insert(str_buf1+':'+str_buf2);
        if (indi_num == indi_ids.size())
            LOGGER.e(0, "Duplicate individual ID found: \"" + str_buf1 + "\t" + str_buf2 + "\".");

        Fam >> str_buf1;
        Fam >> str_buf1;
        Fam >> str_buf1;
        Fam >> str_buf1;
        indi_num++;
    }
    Fam.clear();
    Fam.close();

    unordered_set<string>().swap(indi_ids);
    LOGGER << indi_num << " individuals included" << endl;
}


void taCOJO::read_bedfile(string bedFile, int indi_num, vector<string> &bimSNP, unordered_map<string, vector<char>> &bedData, 
    unordered_map<string, double> &bedSNPavg, unordered_map<string, double> &bedSNPstd) 
{       
    int snp_num = commonSNP.size();

    // some code are adopted from PLINK with modifications
    int i = 0, j = 0, k = 0;
    bool SNP1, SNP2, flag_all_NA;
    double SNP12, SNP_sum, SNP_square_sum, not_NA_indi_num, SNP_avg, SNP_std;
    string SNP_buf;

    // Read bed file
    char ch[1];
    bitset<8> b;

    fstream BIT(bedFile.c_str(), ios::in | ios::binary);
    if (!BIT) LOGGER.e(0, "cannot open the file [" + bedFile + "] to read.");
    LOGGER << "Reading PLINK BED file from [" + bedFile + "] in SNP-major format ..." << endl;
    for (i = 0; i < 3; i++) {BIT.read(ch, 1);} // skip the first three bytes
    
    for (j = 0; j < snp_num; j++) {
        SNP_buf = bimSNP[j];        
        if (commonSNP.find(SNP_buf) == commonSNP.end()) {
            for (i = 0; i < indi_num; i += 4) BIT.read(ch, 1);
            continue;
        }
        
        flag_all_NA = true;
        SNP_sum = 0;
        SNP_square_sum = 0;
        not_NA_indi_num = 0;
        vector<char> single_SNP_temp(indi_num, 0);

        // Read genotype in SNP-major mode, 00: homozygote AA; 11: homozygote BB; 01: hetezygote; 10: missing
        for (i = 0; i < indi_num;) {

            BIT.read(ch, 1);
            if (!BIT) LOGGER.e(0, "problem with the BED file ... has the FAM/BIM file been changed?");
            b = ch[0];
            k = 0;
            while (k < 7 && i < indi_num) { // change code: 11 for AA; 00 for BB;
                SNP2 = (!b[k++]);
                SNP1 = (!b[k++]);
                if (SNP1 && !SNP2)
                    single_SNP_temp[i] = 10;
                else {
                    flag_all_NA = false;
                    SNP12 = SNP1+SNP2;
                    if (bimData[SNP_buf].swap) 
                        SNP12 = 2-SNP12;

                    single_SNP_temp[i] = SNP12;
                    SNP_sum += SNP12;
                    SNP_square_sum += SNP12*SNP12;
                    not_NA_indi_num += 1;
                }
                i++;
            }
        }
        
        if (!flag_all_NA) {
            SNP_avg = SNP_sum/not_NA_indi_num;
            SNP_std = sqrt((SNP_square_sum-SNP_avg*SNP_avg*not_NA_indi_num)/(indi_num-1));
            bedData.insert(pair<string, vector<char>>(SNP_buf, single_SNP_temp));
            bedSNPavg.insert(pair<string, double>(SNP_buf, SNP_avg));
            bedSNPstd.insert(pair<string, double>(SNP_buf, SNP_std));
        } else {
            LOGGER.w(0, "all values are NA for SNP " + SNP_buf + " in bedfile " + bedFile);
            commonSNP.erase(SNP_buf);
        }
    }
    BIT.clear();
    BIT.close();

    vector<string>().swap(bimSNP);
    LOGGER << "Genotype data for " << indi_num << " individuals and " << snp_num << " SNPs included from [" + bedFile + "]." << endl;
}


void taCOJO::process_bedfile(ArrayXXd &X, int indi_num, unordered_map<string, vector<char>> &bedData, 
    unordered_map<string, double> &bedSNPavg, unordered_map<string, double> &bedSNPstd) 
{
    X.resize(commonSNP_num, indi_num);

    # pragma omp parallel for 
    for (int row = 0; row < commonSNP_num; row++) {
        string SNP = commonSNP_ordered[row];
        for (int col = 0; col < indi_num; col++) {
            if (bedData[SNP][col] == 10) 
                X(row, col) = 0;
            else
                X(row, col) = (double(bedData[SNP][col])-bedSNPavg[SNP])/bedSNPstd[SNP];
        }
    };

    unordered_map<string, vector<char>>().swap(bedData);
    unordered_map<string, double>().swap(bedSNPavg);
    unordered_map<string, double>().swap(bedSNPstd);
}


double taCOJO::median(const ArrayXd &eigen_vector) 
{
    if (eigen_vector.size()==1) 
        return eigen_vector(0);

    int size = eigen_vector.size();
    vector<double> b(eigen_vector.data(), eigen_vector.data() + size);
    double b_median; 

    stable_sort(b.begin(), b.end());
    if (size%2==1)
        b_median = b[(size-1)/2];
    else 
        b_median = (b[size/2]+b[size/2-1])/2;

    vector<double>().swap(b);
    return b_median;
}


double taCOJO::process_cojofile(ArrayXXd &sumstat, unordered_map<string, ItemCojo> &cojoData) 
{   
    sumstat.resize(commonSNP_num, 7);
    // col 0:b, 1:se2, 2:p, 3:freq, 4:N, 5:V, 6:D 

    # pragma omp parallel for 
    for (int i = 0; i < commonSNP_num; i++) {
        ItemCojo item = cojoData[commonSNP_ordered[i]];

        sumstat(i, 0) = item.b;
        sumstat(i, 1) = item.se*item.se;
        sumstat(i, 2) = item.p;
        sumstat(i, 3) = item.freq;
        sumstat(i, 4) = item.N;
    };
    
    sumstat.col(5) = sumstat.col(3) * (1-sumstat.col(3)) * 2;
    ArrayXd Vp_gcta_list = sumstat.col(5) * sumstat.col(4) * (sumstat.col(1) + square(sumstat.col(0))/(sumstat.col(4)-1));
    double Vp = median(Vp_gcta_list);

    sumstat.col(4) = (Vp - sumstat.col(5)*square(sumstat.col(0))) / (sumstat.col(5)*sumstat.col(1)) + 1;
    sumstat.col(6) = sumstat.col(4) * sumstat.col(5);

    Vp_gcta_list = sumstat.col(6) * (sumstat.col(1) + square(sumstat.col(0))/(sumstat.col(4)-1));
    Vp = median(Vp_gcta_list);
    
    unordered_map<string, ItemCojo>().swap(cojoData);
    return Vp;
}


void taCOJO::read_files(string cojoFile1, string cojoFile2, string PLINK1, string PLINK2) 
{
    // read cojo and PLINK files 
    read_bimfile(PLINK1+".bim", true);
    read_bimfile(PLINK2+".bim", false);
    cout << endl;

    read_cojofile(cojoFile1, cojoData1);
    read_cojofile(cojoFile2, cojoData2);
    cout << endl;
    
    read_famfile(PLINK1+".fam", indi_num1);
    read_bedfile(PLINK1+".bed", indi_num1, bimSNP1, bedData1, bedSNPavg1, bedSNPstd1);
    cout << endl;

    read_famfile(PLINK2+".fam", indi_num2);
    read_bedfile(PLINK2+".bed", indi_num2, bimSNP2, bedData2, bedSNPavg2, bedSNPstd2);
    cout << endl;

    commonSNP_num = commonSNP.size();
    if (commonSNP_num==0)
        LOGGER.e(0, "Input data has no common SNPs.");

    cout << commonSNP_num << " SNPs included for two cohorts" << endl;

    // assign SNP with fixed order (ID number)
    for (auto iter = commonSNP.begin(); iter != commonSNP.end(); iter++) {
        commonSNP_ordered.push_back(*iter);
    } 

    // clear bim file and two bim file SNP lists
    for (auto iter = bimData.begin(); iter != bimData.end(); iter++) {
        if (commonSNP.find(iter->first) == commonSNP.end())
            bimData.erase(iter);
    } 
    
    unordered_set<string>().swap(commonSNP);

    // initialize matrices and vectors for calculation
    process_bedfile(X1, indi_num1, bedData1, bedSNPavg1, bedSNPstd1);
    process_bedfile(X2, indi_num2, bedData2, bedSNPavg2, bedSNPstd2);

    Vp1 = process_cojofile(sumstat1, cojoData1);
    Vp2 = process_cojofile(sumstat2, cojoData2);
    
    cout << Vp1 << " " << Vp2 << endl;
}
