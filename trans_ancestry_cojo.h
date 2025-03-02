#ifndef _TACOJO_H
#define _TACOJO_H

#include "Logger.h"
#include "CommFunc.h"
#include "StrFunc.h"
#include <omp.h>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <map>
#include <bitset>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SpecialFunctions>

using namespace Eigen;


struct ItemCojo {
    ItemCojo() = default;
    ItemCojo(double freq, double b, double se, double p, double N)
         : freq(freq), b(b), se(se), p(p), N(N) {}

    double freq, b, se, p, N; 
};


struct ItemBim{
    ItemBim() = default;
    ItemBim(string A1, string A2, bool swap): A1(A1), A2(A2), swap(swap) {}

    string A1, A2;
    bool swap;
};


class taCOJO {

public:
    void read_files(string cojoFile1, string cojoFile2, string PLINK1, string PLINK2);
    void main_loop();

// Step 1: read files and prepare matrices and vectors for calculation
public:
    void read_cojofile(string cojoFile, unordered_map<string, ItemCojo> &cojoData);
    void read_bimfile(string bimFile, bool isFirst);
    void read_famfile(string famFile, int &indi_num);
    void read_bedfile(string bedFile, int indi_num, vector<string> &bimSNP, unordered_map<string, vector<char>> &bedData, 
        unordered_map<string, double> &bedSNPavg, unordered_map<string, double> &bedSNPstd);
    void process_bedfile(MatrixXd &X_cohort, int indi_num, unordered_map<string, vector<char>> &bedData, 
        unordered_map<string, double> &bedSNPavg, unordered_map<string, double> &bedSNPstd);    
    void process_cojofile(ArrayXXd &sumstat_cohort, double &Vp_cohort, unordered_map<string, ItemCojo> &cojoData);

private: // memory all freed after Step 1
    unordered_map<string, ItemCojo> cojoData1, cojoData2;    
    unordered_map<string, vector<char>> bedData1, bedData2;
    unordered_map<string, double> bedSNPavg1, bedSNPstd1, bedSNPavg2, bedSNPstd2;
    vector<string> bimSNP1, bimSNP2;
    unordered_set<string> commonSNP;


// Step 2: prepare matrices and vectors for calculation
public:
    int inverse_var_meta(); 
    void calc_conditional_effects();

    void remove_SNP(MatrixXd &matrix);
    void prepare_sumstat(const ArrayXXd &sumstat, ArrayXXd &sumstat_candidate, ArrayXXd &sumstat_screened);
    void record_original_indices();

    void calc_colinear_SNP(const MatrixXd &X, MatrixXd &X_candidate, MatrixXd &r);


public: // all necessary information and data for calculation
    int indi_num1, indi_num2, commonSNP_num;
    vector<string> commonSNP_ordered;
    unordered_map<string, ItemBim> bimData;

    set<int> candidate_SNP, backward_removed_SNP, all_excluded_SNP;
    vector<int> screened_SNP_original_indices;

    MatrixXd X1, X2;
    double Vp1, Vp2;

    // 0:freq, 1:b, 2:se2, 3:p, 4:N, 5:V, 6:D
    ArrayXXd sumstat1, sumstat2;
    
    // 0:freq, 1:b, 2:se2, 3:N, 4:Z, 5:p
    ArrayXXd sumstat_loop;
    
    ArrayXXd conditional_effect_meta;

    MatrixXd X1_candidate, X2_candidate, r1, r2, R1_inv_pre, R2_inv_pre;
    ArrayXXd sumstat1_candidate, sumstat2_candidate, sumstat1_screened, sumstat2_screened;
    



// hyperparameters for users to predefine and adjust
public: 
    double threshold = 5e-8;
    double colinearity_threshold = 0.8;
    double R2_incremental_threshold = 0.0;
    double R2_incremental_threshold_backwards = -0.5;
    int max_iter_num = 200;
};

#endif
