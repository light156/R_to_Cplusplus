#include "trans_ancestry_cojo.h"
#include <iomanip>


int taCOJO::inverse_var_meta() 
{
    // 0:freq, 1:b, 2:se2, 3:N, 4:Z, 5:p
    sumstat_loop.resize(6, commonSNP_num);
    
    // 0:freq, 1:b, 2:se2, 3:p, 4:N, 5:V, 6:D
    sumstat_loop.row(0) = (sumstat1.row(0)*sumstat2.row(2) + sumstat2.row(0)*sumstat1.row(2)) / (sumstat1.row(2)+sumstat2.row(2));
    sumstat_loop.row(1) = (sumstat1.row(1)*sumstat2.row(2)+sumstat2.row(1)*sumstat1.row(2)) / (sumstat1.row(2)+sumstat2.row(2));
    sumstat_loop.row(2) = sumstat1.row(2)*sumstat2.row(2) / (sumstat1.row(2)+sumstat2.row(2));
    sumstat_loop.row(3) = sumstat1.row(4) + sumstat2.row(4);
    sumstat_loop.row(4) = sumstat_loop.row(1) / sqrt(sumstat_loop.row(2));
    sumstat_loop.row(5) = erfc(abs((sumstat_loop.row(4)/sqrt(2))));

    int min_p_value, min_p_index;
    min_p_value = abs(sumstat_loop.row(4)).minCoeff(&min_p_index);
    if (min_p_value>threshold)
        LOGGER.e(0, "Input data has no significant SNPs.");

    return min_p_index;
}


void taCOJO::calc_colinear_SNP(const MatrixXd &X, MatrixXd &X_candidate, MatrixXd &r) 
{
    int i=0;
    X_candidate.resize(X.rows(), candidate_SNP.size());

    for (auto iter = candidate_SNP.begin(); iter != candidate_SNP.end(); iter++, i++) {
        X_candidate.col(i) = X.col(*iter);
    } 

    r = X_candidate.transpose()*X / (X.rows()-1);

    VectorXd temp = r.cwiseAbs().colwise().maxCoeff();
    double colinearity_threshold_sqrt = sqrt(colinearity_threshold);

    for (i = 0; i < commonSNP_num; i++) {
        // remove colinear SNPs
        if (temp(i) >= colinearity_threshold_sqrt) all_excluded_SNP.insert(i);
    }

    temp.resize(0);
}


void taCOJO::remove_SNP(MatrixXd &matrix)
{
    int numRows = matrix.rows(), numCols = matrix.cols(), colToRemove;

    for (auto riter = all_excluded_SNP.rbegin(); riter != all_excluded_SNP.rend(); riter++) {
        numCols--;
        colToRemove = *riter;
        if (colToRemove < numCols)
            matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.rightCols(numCols-colToRemove);
    }

    matrix.conservativeResize(numRows,numCols);
}


void taCOJO::prepare_sumstat(const ArrayXXd &sumstat, ArrayXXd &sumstat_candidate, ArrayXXd &sumstat_screened)
{   
    int i = 0, numRows = sumstat.rows(), numCols = sumstat.cols(), colToRemove;
    sumstat_candidate.resize(numRows, candidate_SNP.size());

    for (auto iter = candidate_SNP.begin(); iter != candidate_SNP.end(); iter++, i++) {
        sumstat_candidate.col(i) = sumstat.col(*iter);
    } 

    sumstat_screened = sumstat;

    for (auto riter = all_excluded_SNP.rbegin(); riter != all_excluded_SNP.rend(); riter++) {
        numCols--;
        colToRemove = *riter;
        if (colToRemove < numCols)
            sumstat_screened.block(0,colToRemove,numRows,numCols-colToRemove) = sumstat_screened.rightCols(numCols-colToRemove);
    }

    sumstat_screened.conservativeResize(numRows,numCols);
}


void taCOJO::record_original_indices() {

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


void taCOJO::calc_conditional_effects() 
{   
    VectorXd V_sqrt, V2_inv_sqrt, conditional_effect_beta1, conditional_effect_beta2;

    V_sqrt = sqrt(sumstat1_candidate.row(5));
    V2_inv_sqrt = sqrt(inverse(sumstat1_screened.row(5)));
    conditional_effect_beta1 = sumstat1_screened.row(1).matrix() - V2_inv_sqrt.asDiagonal()*r1*R1_inv_pre*V_sqrt.asDiagonal()*sumstat1_candidate.row(1).matrix();

    V_sqrt = sqrt(sumstat2_candidate.row(5));
    V2_inv_sqrt = sqrt(inverse(sumstat2_screened.row(5)));
    conditional_effect_beta2 = sumstat2_screened.row(1).matrix() - V2_inv_sqrt.asDiagonal()*r2*R2_inv_pre*V_sqrt.asDiagonal()*sumstat2_candidate.row(1).matrix();

    // 0:bC_ma, 1:se2C_ma, 2:zC_ma, 3:pC_ma
    conditional_effect_meta.resize(4, commonSNP_num);

    conditional_effect_meta.row(0) = (conditional_effect_beta1.array()*sumstat2.row(2) + sumstat2.row(0)*sumstat1.row(2)) / (sumstat1.row(2)+sumstat2.row(2));
    conditional_effect_meta.row(1) = sumstat1_candidate.row(2)*sumstat2_candidate.row(2) / (sumstat1_candidate.row(2)+sumstat2_candidate.row(2));
    conditional_effect_meta.row(2) = conditional_effect_meta.row(0) / sqrt(conditional_effect_meta.row(1));
    conditional_effect_meta.row(3) = erfc(abs((conditional_effect_meta.row(2)/sqrt(2))));
}


void taCOJO::main_loop() {
    
    cout << fixed << setprecision(10) << Vp1 << " " << Vp2 << endl;
    
    int min_p_index = inverse_var_meta();
    candidate_SNP.insert(min_p_index);

    int iter_num = 0;

    while(iter_num < max_iter_num) {
        
        set<int>().swap(all_excluded_SNP);

        calc_colinear_SNP(X1, X1_candidate, r1);
        calc_colinear_SNP(X2, X2_candidate, r2);

        all_excluded_SNP.insert(backward_removed_SNP.begin(), backward_removed_SNP.end());
        record_original_indices();
        
        remove_SNP(r1);
        remove_SNP(r2);
        prepare_sumstat(sumstat1, sumstat1_candidate, sumstat1_screened);
        prepare_sumstat(sumstat2, sumstat2_candidate, sumstat2_screened);
        
        // cout << r1.rows() << " " << r1.cols() << endl;
        // cout << sumstat1_candidate.rows() << " " << sumstat1_candidate.cols() << endl;
        // cout << sumstat1_screened.rows() << " " << sumstat1_screened.cols() << endl;


        break;
    }
}
