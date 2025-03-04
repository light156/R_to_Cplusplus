#include "trans_ancestry_cojo.h"


void taCOJO::initialize()
{   
    threshold = 5e-8;
    colinearity_threshold = 0.8;
    colinearity_threshold_sqrt = sqrt(colinearity_threshold);
    iter_colinearity_threshold = 1 / (1 - colinearity_threshold);
    
    R2_incremental_threshold = 0.0;
    R2_incremental_threshold_backwards = -0.5;
    max_iter_num = 200;  

    R1_inv_pre = MatrixXd::Identity(1,1);
    R2_inv_pre = MatrixXd::Identity(1,1);
}
    

void taCOJO::inverse_var_meta(const ArrayXd &b_cohort1, const ArrayXd &b_cohort2, 
    const ArrayXd &se2_cohort1, const ArrayXd &se2_cohort2, ArrayXXd &merge) 
{
    // merge: 0:b, 1:se2, 2:Z, 3:p
    merge.resize(b_cohort1.size(), 4);

    merge.col(0) = (b_cohort1 * se2_cohort2 + b_cohort2 * se2_cohort1) / (se2_cohort1 + se2_cohort2);
    merge.col(1) = se2_cohort1 * se2_cohort2 / (se2_cohort1 + se2_cohort2);
    merge.col(2) = merge.col(0) / sqrt(merge.col(1));
    merge.col(3) = erfc(abs(merge.col(2)/sqrt(2)));

    abs(merge.col(2)).maxCoeff(&max_SNP_index);
}


void taCOJO::calc_colinear_SNP(const MatrixXd &X, MatrixXd &X_candidate, ArrayXXd &r) 
{
    int i = 0;
    X_candidate.resize(candidate_SNP.size(), X.cols());

    for (auto iter = candidate_SNP.begin(); iter != candidate_SNP.end(); iter++, i++)
        X_candidate.row(i) = X.row(*iter);

    r = X * X_candidate.transpose() / (X.cols()-1);
    ArrayXd temp = abs(r.array()).rowwise().maxCoeff();

    // remove colinear SNPs
    for (i = 0; i < commonSNP_num; i++) 
        if (temp(i) >= colinearity_threshold_sqrt) all_excluded_SNP.insert(i);
}


void taCOJO::remove_SNP(ArrayXXd &matrix)
{
    int numRows = matrix.rows(), numCols = matrix.cols(), rowToRemove;

    for (auto riter = all_excluded_SNP.rbegin(); riter != all_excluded_SNP.rend(); riter++) {
        numRows--;
        rowToRemove = *riter;
        if (rowToRemove < numRows)
            matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.bottomRows(numRows-rowToRemove);
    }

    matrix.conservativeResize(numRows, numCols);
}


void taCOJO::prepare_sumstat(const ArrayXXd &sumstat, ArrayXXd &sumstat_candidate, ArrayXXd &sumstat_screened)
{   
    int i = 0, numRows = sumstat.rows(), numCols = sumstat.cols(), rowToRemove;
    sumstat_candidate.resize(candidate_SNP.size(), numCols);
    sumstat_screened = sumstat;

    for (auto iter = candidate_SNP.begin(); iter != candidate_SNP.end(); iter++, i++)
        sumstat_candidate.row(i) = sumstat.row(*iter);

    for (auto riter = all_excluded_SNP.rbegin(); riter != all_excluded_SNP.rend(); riter++) {
        numRows--;
        rowToRemove = *riter;
        if (rowToRemove < numRows)
            sumstat_screened.block(rowToRemove,0,numRows-rowToRemove,numCols) = sumstat_screened.bottomRows(numRows-rowToRemove);
    }

    sumstat_screened.conservativeResize(numRows, numCols);
}


void taCOJO::record_original_indices() 
{
    vector<int>().swap(screened_SNP_original_indices);
    screened_SNP_original_indices.resize(commonSNP_num-all_excluded_SNP.size());
    int original_index = 0, i = 0;

    for (auto iter = all_excluded_SNP.begin(); iter != all_excluded_SNP.end(); original_index++) {
        if (original_index == *iter) {
            iter++;
        } else {
            screened_SNP_original_indices[i] = original_index;
            i++;
        }
    }

    for (; i<commonSNP_num-all_excluded_SNP.size(); i++, original_index++)
        screened_SNP_original_indices[i] = original_index;
}


void taCOJO::calc_conditional_effects(const ArrayXXd &r, const ArrayXXd &sumstat_candidate, 
    const ArrayXXd &sumstat_screened, const MatrixXd &R_inv_pre, ArrayXd &conditional_effect_beta) 
{   
    MatrixXd temp1, temp2;
    
    temp1 = r1.colwise() * sqrt(inverse(sumstat_screened.col(5)));
    temp2 = sumstat_candidate.col(0) * sqrt(sumstat_candidate.col(5));
    conditional_effect_beta = sumstat_screened.col(0) - (temp1 * R_inv_pre * temp2).array();
}


bool taCOJO::calc_joint_effects(ArrayXXd &sumstat_candidate, const ArrayXXd &sumstat_screened, const MatrixXd &X, 
    MatrixXd &X_candidate, MatrixXd &R_inv_post, double Vp, ArrayXd &beta, ArrayXd &beta_var, double &R2) 
{   
    int M = candidate_SNP.size()+1;

    sumstat_candidate.conservativeResize(M, NoChange);
    sumstat_candidate.row(M-1) = sumstat_screened.row(max_SNP_index);

    X_candidate.conservativeResize(M, NoChange);
    X_candidate.row(M-1) = X.row(screened_SNP_original_indices[max_SNP_index]);

    double Neff = median(sumstat_candidate.col(4));
    R_inv_post.noalias() = (X_candidate*X_candidate.transpose() / (X_candidate.cols()-1)).inverse();

    if ((abs(R_inv_post.minCoeff()) > iter_colinearity_threshold) || (abs(R_inv_post.maxCoeff()) > iter_colinearity_threshold)) 
        return false;

    ArrayXXd temp = (R_inv_post.array().colwise() * sqrt(inverse(sumstat_candidate.col(5)))).rowwise() 
        * sqrt(sumstat_candidate.col(5)).transpose();
    beta = temp.matrix() * sumstat_candidate.col(0).matrix();
    
    double sigma_J_squared = Vp - (sumstat_candidate.col(0) * sumstat_candidate.col(5) * beta).sum();
                
    temp = (R_inv_post.array().colwise() * sqrt(inverse(sumstat_candidate.col(6)))).rowwise() 
        * sqrt(inverse(sumstat_candidate.col(6))).transpose();
    beta_var = temp.matrix().diagonal() * sigma_J_squared;
    
    if (beta_var.minCoeff() <= 0)
        return false;

    R2 = 1 - sigma_J_squared * (Neff-1) / Vp / (Neff-M-1);
    return true;
}


void taCOJO::main_loop() 
{   
    int iter_num = 0;

    initialize();
    inverse_var_meta(sumstat1.col(0), sumstat2.col(0), sumstat1.col(1), sumstat2.col(1), sumstat_merge);
    
    if (sumstat_merge(max_SNP_index, 3) > threshold)
        LOGGER.e(0, "Input data has no significant SNPs.");

    candidate_SNP.insert(max_SNP_index);
    cout << max_SNP_index << " " << commonSNP_ordered[max_SNP_index] << endl;

    loop_break_indicator = false;
    conditional_effect_search_indicator = false;
    bool no_NA_flag = true;

    while (!loop_break_indicator && iter_num<max_iter_num) {
        
        // Step 1: Get submatrices according to candidate_SNP, colinear_SNP and backward_removed_SNP
        set<int>().swap(all_excluded_SNP);
        all_excluded_SNP.insert(backward_removed_SNP.begin(), backward_removed_SNP.end());

        calc_colinear_SNP(X1, X1_candidate, r1);
        calc_colinear_SNP(X2, X2_candidate, r2);
        record_original_indices();

        remove_SNP(r1);
        remove_SNP(r2);
        
        prepare_sumstat(sumstat1, sumstat1_candidate, sumstat1_screened);
        prepare_sumstat(sumstat2, sumstat2_candidate, sumstat2_screened);

        // Step 2: Calculate conditional effects
        calc_conditional_effects(r1, sumstat1_candidate, sumstat1_screened, R1_inv_pre, conditional_effect_beta1);
        calc_conditional_effects(r2, sumstat2_candidate, sumstat2_screened, R2_inv_pre, conditional_effect_beta2);
        
        inverse_var_meta(conditional_effect_beta1, conditional_effect_beta2, sumstat1_screened.col(1), sumstat2_screened.col(1), conditional_effect_meta);

        if (conditional_effect_meta(max_SNP_index, 3) > threshold)
            loop_break_indicator = true;
        else
            conditional_effect_search_indicator = false;
    
        if (loop_break_indicator) {
            LOGGER.w(0, "Calculation finished, NoNewSNP-aboveSigThreshold");
            break;
        }

        /*
        for (int rowToRemove = numRows-1; rowToRemove >= 0; rowToRemove--) {
            if (conditional_effect_meta(rowToRemove, 3) > threshold) {
                numRows--;
                screened_SNP_original_indices.erase(screened_SNP_original_indices.begin()+rowToRemove);

                if (rowToRemove < numRows)
                    conditional_effect_meta.block(rowToRemove,0,numRows-rowToRemove,4) = conditional_effect_meta.bottomRows(numRows-rowToRemove);
            }
        }
        conditional_effect_meta.conservativeResize(numRows, 4);
        */

        while (!conditional_effect_search_indicator) {
            // omit some code here
            
            double R2_cohort1, R2_cohort2;

            no_NA_flag = calc_joint_effects(sumstat1_candidate, sumstat1_screened, X1, X1_candidate, R1_inv_post, Vp1, beta1, beta_var1, R2_cohort1);
            no_NA_flag = calc_joint_effects(sumstat2_candidate, sumstat2_screened, X2, X2_candidate, R2_inv_post, Vp2, beta2, beta_var2, R2_cohort2);
                
            inverse_var_meta(beta1, beta2, beta_var1, beta_var2, joint_effect_meta);
            cout << R2_cohort1 << " " << R2_cohort2 << endl;
            break;  
        }

        break;
    }
}
