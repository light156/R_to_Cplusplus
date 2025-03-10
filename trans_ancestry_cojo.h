#ifndef _TACOJO_H
#define _TACOJO_H

#include "Logger.h"
#include "StrFunc.h"
#include <omp.h>
#include <unordered_map>
#include <unordered_set>
#include <fstream>
#include <bitset>
#include <Eigen/Dense>
#include <unsupported/Eigen/SpecialFunctions>
#include <iomanip>
#include <numeric>

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
    void save_results(string filepath);
    void main_loop();

// Step 1: read files and prepare matrices and vectors for calculation
public:
    void read_cojofile(string cojoFile, bool isFirst);
    void read_bimfile(string bimFile, bool isFirst);
    int read_famfile(string famFile);
    void read_bedfile(string bedFile, int indi_num, bool isFirst, 
        vector<string> &bimSNP, unordered_map<string, vector<double>> &bedData);
    double generate_sumstat_matrix(ArrayXXd &sumstat, unordered_map<string, ItemCojo> &cojoData);
    double median(const ArrayXd &vector);

private: // memory all freed after Step 1
    unordered_map<string, ItemCojo> cojoData1, cojoData2;    
    unordered_map<string, vector<double>> bedData1, bedData2;
    vector<string> bimSNP1, bimSNP2;
    unordered_set<string> commonSNP;

// Step 2: main loop
public:
    void initialize_hyperparameters();
    void initialize_matrices();
    
    void inverse_var_meta(const ArrayXd &b_cohort1, const ArrayXd &b_cohort2, 
        const ArrayXd &se2_cohort1, const ArrayXd &se2_cohort2, ArrayXXd &merge);
    void calc_conditional_effects(ArrayXd &conditional_beta1, ArrayXd &conditional_beta2);
    bool calc_joint_effects(const ArrayXXd &sumstat_candidate, const MatrixXd &R_inv_post, 
        double Vp, ArrayXd &beta, ArrayXd &beta_var, double &R2, bool flag);
        
    void append_row(ArrayXXd &matrix, const ArrayXXd &matrix_new);
    void append_row(MatrixXd &matrix, const MatrixXd &matrix_new);
    void append_column(MatrixXd &matrix, const MatrixXd &matrix_new); 
    void remove_row(ArrayXXd &matrix, int index);
    void remove_row(MatrixXd &matrix, int index);
    void remove_column(MatrixXd &matrix, int index); 
    void fast_inv(const MatrixXd &R_inv_pre, const VectorXd &new_column, MatrixXd &R_inv_post);
    
    void remove_new_colinear_SNP();
    void refuse_max_SNP_as_candidate();
    
    void initialize_backward_selection();
    void adjust_SNP_according_to_backward_selection();

public: // all necessary data and results during the loop
    int indi_num1, indi_num2, commonSNP_num, max_SNP_index;
    vector<string> commonSNP_ordered;
    unordered_map<string, ItemBim> bimData;
    vector<int> candidate_SNP, backward_removed_SNP, screened_SNP, excluded_SNP;
    
    // bed matrix
    MatrixXd X1, X2, X1_candidate, X2_candidate, X1_screened, X2_screened, X1_new_model, X2_new_model;
    // cols 0:b, 1:se2, 2:p, 3:freq, 4:N, 5:V, 6:D 
    ArrayXXd sumstat1, sumstat2, sumstat1_candidate, sumstat2_candidate, 
        sumstat1_screened, sumstat2_screened, sumstat1_new_model, sumstat2_new_model;
    
    // merge: 0:b, 1:se2, 2:Zabs, 3:p
    ArrayXXd sumstat_merge, new_model_joint;

    // outputted joint vectors
    ArrayXd output_b_cohort1, output_b_cohort2, output_se2_cohort1, output_se2_cohort2;
    
    MatrixXd r1, r2;
    MatrixXd R1_inv_pre, R2_inv_pre, R1_inv_post, R2_inv_post;
    double Vp1, Vp2;


// hyperparameters for users to predefine and adjust
public: 
    double threshold;
    double colinearity_threshold, colinearity_threshold_sqrt, iter_colinearity_threshold;

    double R2_incremental_threshold;
    double R2_incremental_threshold_backwards;

    int max_iter_num;
};

#endif
