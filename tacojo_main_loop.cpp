#include "trans_ancestry_cojo.h"


void taCOJO::initialize_hyperparameters()
{   
    threshold = 5e-8;
    colinearity_threshold = 0.8;
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
    // merge: 0:b, 1:se2, 2:Zabs, 3:p
    merge.resize(b_cohort1.size(), 4);

    merge.col(0) = (b_cohort1 * se2_cohort2 + b_cohort2 * se2_cohort1) / (se2_cohort1 + se2_cohort2);
    merge.col(1) = se2_cohort1 * se2_cohort2 / (se2_cohort1 + se2_cohort2);
    merge.col(2) = abs(merge.col(0) / sqrt(merge.col(1)));
    merge.col(3) = erfc(merge.col(2)/sqrt(2));
}


void taCOJO::calc_colinear_SNP(const ArrayXXd &X, ArrayXXd &X_candidate, ArrayXXd &r) 
{
    int i = 0;
    X_candidate.resize(candidate_SNP.size(), X.cols());
    
    for (auto iter = candidate_SNP.begin(); iter != candidate_SNP.end(); iter++, i++)
        X_candidate.row(i) = X.row(*iter);

    r = X.matrix() * X_candidate.transpose().matrix() / (X.cols()-1);
    ArrayXd temp = square(r).rowwise().maxCoeff();

    // remove colinear SNPs
    for (i = 0; i < commonSNP_num; i++) {
        if (temp(i) >= colinearity_threshold) {
            all_excluded_SNP.insert(i);
        }
    }
}


void taCOJO::remove_SNP(ArrayXXd &matrix, int rowToRemove)
{
    int numRows = matrix.rows()-1, numCols = matrix.cols();

    if (rowToRemove < numRows)
        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.bottomRows(numRows-rowToRemove);

    matrix.conservativeResize(numRows, numCols);
}


void taCOJO::prepare_sumstat(const ArrayXXd &sumstat, ArrayXXd &sumstat_candidate, ArrayXXd &sumstat_screened)
{   
    int i = 0, numRows = sumstat.rows(), numCols = sumstat.cols(), rowToRemove;
    sumstat_candidate.resize(candidate_SNP.size(), numCols);
    sumstat_screened = sumstat;

    for (auto iter = candidate_SNP.begin(); iter != candidate_SNP.end(); iter++, i++)
        sumstat_candidate.row(i) = sumstat.row(*iter);

    for (auto riter = all_excluded_SNP.rbegin(); riter != all_excluded_SNP.rend(); riter++) 
        remove_SNP(sumstat_screened, *riter);
}


void taCOJO::initialize_candidate_and_screened_matrices() 
{
    set<int>().swap(all_excluded_SNP);
    calc_colinear_SNP(X1, X1_candidate, r1);
    calc_colinear_SNP(X2, X2_candidate, r2);
    all_excluded_SNP.insert(backward_removed_SNP.begin(), backward_removed_SNP.end());

    // remove the excluded rows in r and sumstat matrices
    for (auto riter = all_excluded_SNP.rbegin(); riter != all_excluded_SNP.rend(); riter++) {
        remove_SNP(r1, *riter);
        remove_SNP(r2, *riter);
    }
    
    prepare_sumstat(sumstat1, sumstat1_candidate, sumstat1_screened);
    prepare_sumstat(sumstat2, sumstat2_candidate, sumstat2_screened);

    // keep track of original indices of all screened indices
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

    temp1 = r.colwise() * inverse(sqrt(sumstat_screened.col(5)));
    temp2 = sumstat_candidate.col(0) * sqrt(sumstat_candidate.col(5));
    conditional_effect_beta = sumstat_screened.col(0) - (temp1 * R_inv_pre * temp2).array();
}


bool taCOJO::calc_joint_effects(ArrayXXd &sumstat_candidate, const ArrayXXd &sumstat_screened, const ArrayXXd &X, 
    ArrayXXd &X_candidate, MatrixXd &R_inv_post, double Vp, ArrayXd &beta, ArrayXd &beta_var, double &R2) 
{   
    int M = sumstat_candidate.rows()+1;

    sumstat_candidate.conservativeResize(M, NoChange);
    sumstat_candidate.row(M-1) = sumstat_screened.row(max_SNP_index);

    X_candidate.conservativeResize(M, NoChange);
    X_candidate.row(M-1) = X.row(screened_SNP_original_indices[max_SNP_index]);

    double Neff = median(sumstat_candidate.col(4));
    R_inv_post.noalias() = (X_candidate.matrix() * X_candidate.transpose().matrix() / (X_candidate.cols()-1)).inverse();

    if ((abs(R_inv_post.minCoeff()) > iter_colinearity_threshold) || (abs(R_inv_post.maxCoeff()) > iter_colinearity_threshold)) 
        return true;

    ArrayXXd temp = (R_inv_post.array().colwise() * inverse(sqrt(sumstat_candidate.col(5)))).rowwise() 
        * sqrt(sumstat_candidate.col(5)).transpose();
    beta = temp.matrix() * sumstat_candidate.col(0).matrix();
    
    double sigma_J_squared = Vp - (sumstat_candidate.col(0) * sumstat_candidate.col(5) * beta).sum();
                
    temp = (R_inv_post.array().colwise() * inverse(sqrt(sumstat_candidate.col(6)))).rowwise() 
        * inverse(sqrt(sumstat_candidate.col(6))).transpose();
    beta_var = temp.matrix().diagonal() * sigma_J_squared;
    
    if (beta_var.minCoeff() <= 0)
        return true;

    R2 = 1 - sigma_J_squared * (Neff-1) / Vp / (Neff-M-1);
    return false;
}


void taCOJO::remove_max_SNP_from_matrices(bool both_cohorts) 
{
    if (sumstat_merge.rows() <= 1) {
        max_SNP_index = -1;
        return;
    }

    int M = sumstat1_candidate.rows();

    screened_SNP_original_indices.erase(screened_SNP_original_indices.begin()+max_SNP_index);
    remove_SNP(sumstat_merge, max_SNP_index);

    remove_SNP(sumstat1_screened, max_SNP_index);
    remove_SNP(r1, max_SNP_index);
    remove_SNP(sumstat1_candidate, M);
    remove_SNP(X1_candidate, M);

    if (both_cohorts) {
        remove_SNP(sumstat2_screened, max_SNP_index);
        remove_SNP(r2, max_SNP_index);
        remove_SNP(sumstat2_candidate, M);
        remove_SNP(X2_candidate, M);
    }

    sumstat_merge.col(2).maxCoeff(&max_SNP_index);
}


void taCOJO::save_file(string bimfile) {
    ofstream Bim(bimfile.c_str());
    if (!Bim) LOGGER.e(0, "cannot open the file [" + bimfile + "] to write.");
    LOGGER << "Writing PLINK BIM file to [" + bimfile + "] ..." << endl;

    ItemBim temp_item;

    for (auto iter=commonSNP.begin(); iter!=commonSNP.end(); iter++) {
        temp_item = bimData[*iter]; 
        Bim << *iter << "\t" << temp_item.A1 << "\t" << temp_item.A2 << endl;
    }

    Bim.close();
}


void taCOJO::main_loop() 
{   
    initialize_hyperparameters();

    output_b_cohort1 = sumstat1.col(0);
    output_b_cohort2 = sumstat2.col(0);
    output_se2_cohort1 = sumstat1.col(1);
    output_se2_cohort2 = sumstat2.col(1);

    inverse_var_meta(output_b_cohort1, output_b_cohort2, output_se2_cohort1, output_se2_cohort2, sumstat_merge);

    sumstat_merge.col(2).maxCoeff(&max_SNP_index);
    if (sumstat_merge(max_SNP_index, 3) > threshold)
        LOGGER.e(0, "Input data has no significant SNPs.");

    candidate_SNP.push_back(max_SNP_index);
    cout << max_SNP_index << " " << commonSNP_ordered[max_SNP_index] << endl;

    ArrayXd conditional_effect_beta1, conditional_effect_beta2;
    ArrayXd beta1, beta2, beta_var1, beta_var2;
    double R2_cohort1, R2_cohort2, previous_R2_cohort1 = 0.0, previous_R2_cohort2 = 0.0;

    bool NA_flag = false, loop_break_indicator = false;
    int iter_num = 0;

    while (!loop_break_indicator && iter_num<20) {
        
        initialize_candidate_and_screened_matrices();

        calc_conditional_effects(r1, sumstat1_candidate, sumstat1_screened, R1_inv_pre, conditional_effect_beta1);
        calc_conditional_effects(r2, sumstat2_candidate, sumstat2_screened, R2_inv_pre, conditional_effect_beta2);

        inverse_var_meta(conditional_effect_beta1, conditional_effect_beta2, 
            sumstat1_screened.col(1), sumstat2_screened.col(1), sumstat_merge);
        
        sumstat_merge.col(2).maxCoeff(&max_SNP_index);
        if (sumstat_merge(max_SNP_index, 3) > threshold) {
            LOGGER.w(0, "Calculation finished, NoNewSNP-aboveSigThreshold");
            break;
        }
        
        while (true) {
            if (max_SNP_index == -1 || sumstat_merge(max_SNP_index, 3) > threshold) {
                LOGGER.w(0, "Calculation finished, NoNewSNP-notMeetingCriteria");
                loop_break_indicator = true;
                break;
            }
            
            cout << "Newly identified SNP: " << commonSNP_ordered[screened_SNP_original_indices[max_SNP_index]] << endl;
            NA_flag = calc_joint_effects(sumstat1_candidate, sumstat1_screened, X1, X1_candidate, R1_inv_post, Vp1, beta1, beta_var1, R2_cohort1);
            
            if (NA_flag) {
                LOGGER.w(0, "NA produced in new joint model, potentially due to colinearity, searching next SNP");
                remove_max_SNP_from_matrices(false);
                continue;
            }
            
            NA_flag = calc_joint_effects(sumstat2_candidate, sumstat2_screened, X2, X2_candidate, R2_inv_post, Vp2, beta2, beta_var2, R2_cohort2);
            
            if (NA_flag) {
                LOGGER.w(0, "NA produced in new joint model, potentially due to colinearity, searching next SNP");
                remove_max_SNP_from_matrices(true);
                continue;
            }
                
            inverse_var_meta(beta1, beta2, beta_var1, beta_var2, new_model_joint);
            
            if ((R2_cohort1 < (1+R2_incremental_threshold)*previous_R2_cohort1) || (R2_cohort2 < (1+R2_incremental_threshold)*previous_R2_cohort2)) {
                LOGGER.w(0, "R2 increment unsatisfactory, searching next SNP");
                remove_max_SNP_from_matrices(true);
                continue;
            } 
            
            if (new_model_joint.col(3).maxCoeff() < threshold) {
                LOGGER.i(0, "Candidate SNP has passed all checks in this iteration");

                previous_R2_cohort1 = R2_cohort1;
                previous_R2_cohort2 = R2_cohort2;
                R1_inv_pre = R1_inv_post;
                R2_inv_pre = R2_inv_post;

                // cout << "Added R vector cohort 1: " << r1.row(max_SNP_index) << endl;
                // cout << "Added R vector cohort 2: " << r2.row(max_SNP_index) << endl;

                cout << "Added diagnoal value cohort 1: " << R1_inv_post(R1_inv_post.rows()-1, R1_inv_post.cols()-1) << endl;
                cout << "Added diagnoal value cohort 2: " << R2_inv_post(R2_inv_post.rows()-1, R2_inv_post.cols()-1) << endl;

                cout << "Joint b: " << new_model_joint(new_model_joint.rows()-1, 0) << endl;
                cout << "Joint se: " << sqrt(new_model_joint(new_model_joint.rows()-1, 1)) << endl;
                cout << "Joint p-value: " << scientific << new_model_joint(new_model_joint.rows()-1, 3) << endl;

                cout << "Adjusted R2 for cohort 1: " << fixed << R2_cohort1 << endl;
                cout << "Adjusted R2 for cohort 2: " << R2_cohort2 << endl;
                
                output_b_cohort1 = beta1;
                output_b_cohort2 = beta2;
                output_se2_cohort1 = beta_var1;
                output_se2_cohort2 = beta_var2; 

                candidate_SNP.push_back(screened_SNP_original_indices[max_SNP_index]);
                remove_max_SNP_from_matrices(true);
                break; 
            }
        }

        iter_num++;
        cout << "iter " << iter_num << " finished" << endl;
        cout << "--------------------------------" << endl;
    }

    inverse_var_meta(output_b_cohort1, output_b_cohort2, output_se2_cohort1, output_se2_cohort2, sumstat_merge);
    cout << "bJ.ma: " << sumstat_merge.col(0).transpose() << endl;
    cout << "seJ.ma: " << sqrt(sumstat_merge.col(1)).transpose() << endl;
    cout << "zJ: " << (sumstat_merge.col(0)/sqrt(sumstat_merge.col(1))).transpose() << endl;
    cout << "pJ: " << scientific << sumstat_merge.col(3).transpose() << endl;
}
