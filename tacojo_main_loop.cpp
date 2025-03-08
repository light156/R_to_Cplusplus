#include "trans_ancestry_cojo.h"


void taCOJO::initialize_hyperparameters()
{   
    threshold = 5e-8;

    colinearity_threshold = 0.8;
    colinearity_threshold_sqrt = sqrt(colinearity_threshold);
    iter_colinearity_threshold = 1 / (1 - colinearity_threshold);
    
    R2_incremental_threshold = 0.0;
    R2_incremental_threshold_backwards = -0.5;
    max_iter_num = 200;  
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


void taCOJO::remove_row(ArrayXXd &matrix, int index=-1)
{   
    // -1 indicates the last row
    int numRows = matrix.rows()-1, numCols = matrix.cols();

    if (index != -1 && index < numRows)
        matrix.block(index, 0, numRows-index, numCols) = matrix.bottomRows(numRows-index);

    matrix.conservativeResize(numRows, numCols);
}


void taCOJO::remove_column(MatrixXd &matrix, int index=-1)
{   
    // -1 indicates the last column
    int numRows = matrix.rows(), numCols = matrix.cols()-1;

    if (index != -1 && index < numCols)
        matrix.block(0, index, numRows, numCols-index) = matrix.rightCols(numCols-index);

    matrix.conservativeResize(numRows, numCols);
}


void taCOJO::initialize_matrices() 
{   
    vector<int>().swap(screened_SNP_original_indices);
    screened_SNP_original_indices.resize(commonSNP_num);
    iota(screened_SNP_original_indices.begin(), screened_SNP_original_indices.end(), 0);

    int M = candidate_SNP.size(), temp_index;

    sumstat1_candidate.resize(M, sumstat1.cols());
    sumstat2_candidate.resize(M, sumstat2.cols());
    X1_candidate.resize(indi_num1, M);
    X2_candidate.resize(indi_num2, M);

    for (int i = 0; i < M; i++) {
        temp_index = candidate_SNP[i];
        sumstat1_candidate.row(i) = sumstat1.row(temp_index);
        sumstat2_candidate.row(i) = sumstat2.row(temp_index);
        X1_candidate.col(i) = X1.col(temp_index);
        X2_candidate.col(i) = X2.col(temp_index);
    }

    sumstat1_screened = sumstat1;
    sumstat2_screened = sumstat2;
    X1_screened = X1;
    X2_screened = X2;

    auto riter = backward_removed_SNP.rbegin();
    while (riter != backward_removed_SNP.rend()) {
        temp_index = *riter;
        remove_row(sumstat1_screened, temp_index);
        remove_row(sumstat2_screened, temp_index);
        remove_column(X1_screened, temp_index);
        remove_column(X2_screened, temp_index);
        screened_SNP_original_indices.erase(screened_SNP_original_indices.begin()+temp_index);
        riter++;
    }

    r1 = X1_candidate.transpose() * X1_screened / (indi_num1-1);
    r2 = X2_candidate.transpose() * X2_screened / (indi_num2-1);

    // calculate colinear SNP
    ArrayXd temp1 = r1.cwiseAbs().colwise().maxCoeff(), temp2 = r2.cwiseAbs().colwise().maxCoeff();

    for (int i = r1.cols()-1; i >= 0; i--) {
        if (temp1(i) >= colinearity_threshold_sqrt || temp2(i) >= colinearity_threshold_sqrt) {
            remove_row(sumstat1_screened, i);
            remove_row(sumstat2_screened, i);
            remove_column(X1_screened, i);
            remove_column(X2_screened, i);            
            remove_column(r1, i);
            remove_column(r2, i);
            screened_SNP_original_indices.erase(screened_SNP_original_indices.begin()+i);
        }
    }
}


void taCOJO::calc_conditional_effects(ArrayXd &conditional_beta1, ArrayXd &conditional_beta2) 
{   
    MatrixXd temp1, temp2;
    ArrayXd temp3;

    temp1 = sumstat1_candidate.col(0) * sqrt(sumstat1_candidate.col(5));
    temp2.noalias() = R1_inv_pre * temp1;
    temp3 = r1.transpose() * temp2;
    conditional_beta1 = sumstat1_screened.col(0) - temp3 / sqrt(sumstat1_screened.col(5));
    
    temp1 = sumstat2_candidate.col(0) * sqrt(sumstat2_candidate.col(5));
    temp2.noalias() = R2_inv_pre * temp1;
    temp3 = r2.transpose() * temp2;
    conditional_beta2 = sumstat2_screened.col(0) - temp3 / sqrt(sumstat2_screened.col(5));
}


void taCOJO::fast_inv(const MatrixXd &R_inv_pre, const VectorXd &new_column, MatrixXd &R_inv_post) {
    double temp_element = 1 / (1 - new_column.transpose() * R_inv_pre * new_column);
    VectorXd temp_vector = R_inv_pre * new_column;

    int dim = R_inv_pre.rows();
    R_inv_post.resize(dim+1, dim+1);

    R_inv_post.block(0, 0, dim, dim) = R_inv_pre + temp_element * temp_vector * temp_vector.transpose();
    R_inv_post.block(0, dim, dim, 1) = -temp_element * temp_vector;
    R_inv_post.block(dim, 0, 1, dim) = -temp_element * temp_vector.transpose();
    R_inv_post(dim, dim) = temp_element;
}


bool taCOJO::calc_joint_effects(const ArrayXXd &sumstat_candidate, const MatrixXd &R_inv_post, 
    double Vp, ArrayXd &beta, ArrayXd &beta_var, double &R2, bool flag) 
{       
    if (flag && ((abs(R_inv_post.minCoeff()) > iter_colinearity_threshold) || (abs(R_inv_post.maxCoeff()) > iter_colinearity_threshold)))
        return true;

    VectorXd temp1 = sqrt(sumstat_candidate.col(5)) * sumstat_candidate.col(0);
    ArrayXd temp2 = R_inv_post * temp1;
    beta = temp2 / sqrt(sumstat_candidate.col(5));
    
    double sigma_J_squared = Vp - (sumstat_candidate.col(0) * sumstat_candidate.col(5) * beta).sum();
                
    beta_var = R_inv_post.diagonal().array() / sumstat_candidate.col(6) * sigma_J_squared;

    if (flag && beta_var.minCoeff() <= 0)
        return true;

    double Neff = median(sumstat_candidate.col(4));
    int M = R_inv_post.rows();
    R2 = 1 - sigma_J_squared * (Neff-1) / Vp / (Neff-M-1);
    return false;
}


void taCOJO::refuse_max_SNP_as_candidate() 
{
    // cohort 1
    remove_row(sumstat1_candidate);

    // cohort 2
    if (sumstat1_candidate.rows() != sumstat2_candidate.rows())
        remove_row(sumstat2_candidate);

    // erfc(0) = 1
    sumstat_merge(max_SNP_index, 2) = 0;
}


void taCOJO::accept_max_SNP_as_candidate() 
{       
    int M = candidate_SNP.size();
    X1_candidate.conservativeResize(NoChange, M);
    X2_candidate.conservativeResize(NoChange, M);
    X1_candidate.col(M-1) = X1_screened.col(max_SNP_index);
    X2_candidate.col(M-1) = X2_screened.col(max_SNP_index);

    MatrixXd r1_temp_vec = X1_screened.col(max_SNP_index).transpose() * X1_screened / (indi_num1 - 1); 
    MatrixXd r2_temp_vec = X2_screened.col(max_SNP_index).transpose() * X2_screened / (indi_num2 - 1); 
    
    r1.conservativeResize(M, NoChange);
    r2.conservativeResize(M, NoChange);
    r1.row(M-1) = r1_temp_vec;
    r2.row(M-1) = r2_temp_vec;
    
    // calculate colinear SNP
    for (int i = r1.cols()-1; i >= 0; i--) {
        if (abs(r1_temp_vec(i)) >= colinearity_threshold_sqrt || abs(r2_temp_vec(i)) >= colinearity_threshold_sqrt) {
            remove_row(sumstat1_screened, i);
            remove_row(sumstat2_screened, i);
            remove_column(r1, i);
            remove_column(r2, i);
            remove_column(X1_screened, i);
            remove_column(X2_screened, i);
            screened_SNP_original_indices.erase(screened_SNP_original_indices.begin()+i);
        }
    }
}


void taCOJO::initialize_backward_selection() 
{   
    int M = candidate_SNP.size()+1;

    X1_new_model.resize(indi_num1, M);
    X2_new_model.resize(indi_num2, M);
    X1_new_model << X1_candidate, X1_screened.col(max_SNP_index);
    X2_new_model << X2_candidate, X2_screened.col(max_SNP_index);

    sumstat1_new_model = sumstat1_candidate;
    sumstat2_new_model = sumstat2_candidate;

    for (int i = new_model_joint.rows()-1; i >= 0; i--) {
        if (new_model_joint(i, 3) > threshold) {
            remove_row(sumstat1_new_model, i);
            remove_row(sumstat2_new_model, i);
            remove_column(X1_new_model, i);
            remove_column(X2_new_model, i);
        }
    }

    R1_inv_post.noalias() = (X1_new_model.transpose() * X1_new_model).inverse() * (indi_num1-1);
    R2_inv_post.noalias() = (X2_new_model.transpose() * X2_new_model).inverse() * (indi_num2-1);
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
    output_b_cohort1 = sumstat1.col(0);
    output_b_cohort2 = sumstat2.col(0);
    output_se2_cohort1 = sumstat1.col(1);
    output_se2_cohort2 = sumstat2.col(1);
    R1_inv_pre = MatrixXd::Identity(1,1);
    R2_inv_pre = MatrixXd::Identity(1,1);

    inverse_var_meta(output_b_cohort1, output_b_cohort2, output_se2_cohort1, output_se2_cohort2, sumstat_merge);

    sumstat_merge.col(2).maxCoeff(&max_SNP_index);
    if (sumstat_merge(max_SNP_index, 3) > threshold)
        LOGGER.e(0, "Input data has no significant SNPs.");

    candidate_SNP.push_back(max_SNP_index);
    cout << "First SNP: " << commonSNP_ordered[max_SNP_index] << endl;

    ArrayXd conditional_beta1, conditional_beta2;
    ArrayXd beta1, beta2, beta_var1, beta_var2;
    double R2_cohort1, R2_cohort2, previous_R2_cohort1 = 0.0, previous_R2_cohort2 = 0.0;

    bool NA_flag = false, loop_break_indicator = false;
    int iter_num = 0;

    initialize_matrices();

    while (!loop_break_indicator && iter_num<max_iter_num) {
        
        calc_conditional_effects(conditional_beta1, conditional_beta2);
        inverse_var_meta(conditional_beta1, conditional_beta2, sumstat1_screened.col(1), sumstat2_screened.col(1), sumstat_merge);

        while (true) {
            // select maximal SNP
            sumstat_merge.col(2).maxCoeff(&max_SNP_index);
            if (sumstat_merge(max_SNP_index, 3) > threshold) {
                LOGGER.w(0, "Calculation finished, NoNewSNP-aboveSigThreshold");
                loop_break_indicator = true;
                break;
            }

            string temp_SNP_name = commonSNP_ordered[screened_SNP_original_indices[max_SNP_index]]; 
            
            int M = candidate_SNP.size()+1;
            sumstat1_candidate.conservativeResize(M, NoChange);
            sumstat2_candidate.conservativeResize(M, NoChange);
            sumstat1_candidate.row(M-1) = sumstat1_screened.row(max_SNP_index);
            sumstat2_candidate.row(M-1) = sumstat2_screened.row(max_SNP_index);
        
            fast_inv(R1_inv_pre, X1_screened.col(max_SNP_index).transpose() * X1_candidate / (indi_num1-1), R1_inv_post);
            fast_inv(R2_inv_pre, X2_screened.col(max_SNP_index).transpose() * X2_candidate / (indi_num2-1), R2_inv_post);
            
            NA_flag = calc_joint_effects(sumstat1_candidate, R1_inv_post, Vp1, beta1, beta_var1, R2_cohort1, true);            
            if (NA_flag) {
                // LOGGER.w(0, "NA produced, potentially due to colinearity", temp_SNP_name);
                refuse_max_SNP_as_candidate();
                continue;
            }
            
            NA_flag = calc_joint_effects(sumstat2_candidate, R2_inv_post, Vp2, beta2, beta_var2, R2_cohort2, true);
            if (NA_flag) {
                // LOGGER.w(0, "NA produced, potentially due to colinearity", temp_SNP_name);
                refuse_max_SNP_as_candidate();
                continue;
            }

            inverse_var_meta(beta1, beta2, beta_var1, beta_var2, new_model_joint);
            if ((R2_cohort1 < (1+R2_incremental_threshold)*previous_R2_cohort1) || 
                (R2_cohort2 < (1+R2_incremental_threshold)*previous_R2_cohort2)) {
                // LOGGER.w(0, "R2 increment unsatisfactory", temp_SNP_name);
                refuse_max_SNP_as_candidate();
                continue;
            }

            if (new_model_joint.col(3).maxCoeff() <= threshold) {
                LOGGER.i(0, "All checks passed in this iteration", temp_SNP_name);
                candidate_SNP.push_back(screened_SNP_original_indices[max_SNP_index]);
                accept_max_SNP_as_candidate();

                // cout << "Added R vector cohort 1: " << r1.row(max_SNP_index) << endl;
                // cout << "Added R vector cohort 2: " << r2.row(max_SNP_index) << endl;
                cout << "Added diagnoal value cohort 1: " << R1_inv_post(M-1, M-1) << endl;
                cout << "Added diagnoal value cohort 2: " << R2_inv_post(M-1, M-1) << endl;
                cout << "Joint b: " << new_model_joint(M-1, 0) << endl;
                cout << "Joint se: " << sqrt(new_model_joint(M-1, 1)) << endl;
                cout << "Joint p-value: " << scientific << new_model_joint(M-1, 3) << endl;
                cout << "Adjusted R2 for cohort 1: " << fixed << R2_cohort1 << endl;
                cout << "Adjusted R2 for cohort 2: " << R2_cohort2 << endl;
                break; 
            }

            LOGGER.w(0, "Backward selection", temp_SNP_name);
            clock_t tStart = clock();
            initialize_backward_selection();
            calc_joint_effects(sumstat1_new_model, R1_inv_post, Vp1, beta1, beta_var1, R2_cohort1, false);
            calc_joint_effects(sumstat2_new_model, R2_inv_post, Vp2, beta2, beta_var2, R2_cohort2, false);
            printf("Step 1: Time taken: %.5fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

            cout << "Adjusted R2 for cohort 1: " << fixed << R2_cohort1 << endl;
            cout << "Adjusted R2 for cohort 2: " << R2_cohort2 << endl;

            if ((R2_cohort1 < (1+R2_incremental_threshold_backwards)*previous_R2_cohort1) || 
                (R2_cohort2 < (1+R2_incremental_threshold_backwards)*previous_R2_cohort2)) {
                LOGGER.d(0, "Adjusted R2 lower than backward selection threshold", temp_SNP_name);
                refuse_max_SNP_as_candidate();
                continue;
            }  
            
            LOGGER.i(0, "Backward selection complete", temp_SNP_name);
            candidate_SNP.push_back(screened_SNP_original_indices[max_SNP_index]);

            // move bad candidate SNPs into backward removed SNP list
            for (int i = candidate_SNP.size()-1; i >= 0; i--) {
                if (new_model_joint(i, 3) > threshold) {
                    cout << new_model_joint.rows() << " " << i << " " << commonSNP_ordered[candidate_SNP[i]] << endl;
                    backward_removed_SNP.insert(candidate_SNP[i]);
                    candidate_SNP.erase(candidate_SNP.begin()+i);
                }
            }

            tStart = clock();
            vector<int>().swap(screened_SNP_original_indices);
            screened_SNP_original_indices.resize(commonSNP_num);
            iota(screened_SNP_original_indices.begin(), screened_SNP_original_indices.end(), 0);

            M = candidate_SNP.size();
            int temp_index;

            sumstat1_candidate.resize(M, sumstat1.cols());
            sumstat2_candidate.resize(M, sumstat2.cols());
            X1_candidate.resize(indi_num1, M);
            X2_candidate.resize(indi_num2, M);

            for (int i = 0; i < M; i++) {
                temp_index = candidate_SNP[i];
                sumstat1_candidate.row(i) = sumstat1.row(temp_index);
                sumstat2_candidate.row(i) = sumstat2.row(temp_index);
                X1_candidate.col(i) = X1.col(temp_index);
                X2_candidate.col(i) = X2.col(temp_index);
            }

            sumstat1_screened = sumstat1;
            sumstat2_screened = sumstat2;
            X1_screened = X1;
            X2_screened = X2;

            auto riter = backward_removed_SNP.rbegin();
            while (riter != backward_removed_SNP.rend()) {
                temp_index = *riter;
                remove_row(sumstat1_screened, temp_index);
                remove_row(sumstat2_screened, temp_index);
                remove_column(X1_screened, temp_index);
                remove_column(X2_screened, temp_index);
                screened_SNP_original_indices.erase(screened_SNP_original_indices.begin()+temp_index);
                riter++;
            }

            r1 = X1_candidate.transpose() * X1_screened / (indi_num1-1);
            r2 = X2_candidate.transpose() * X2_screened / (indi_num2-1);
            printf("Step 2: Time taken: %.5fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

            tStart = clock();
            // calculate colinear SNP
            ArrayXd temp1 = r1.cwiseAbs().colwise().maxCoeff(), temp2 = r2.cwiseAbs().colwise().maxCoeff();

            for (int i = r1.cols()-1; i >= 0; i--) {
                if (temp1(i) >= colinearity_threshold_sqrt || temp2(i) >= colinearity_threshold_sqrt) {
                    remove_row(sumstat1_screened, i);
                    remove_row(sumstat2_screened, i);
                    remove_column(X1_screened, i);
                    remove_column(X2_screened, i);            
                    remove_column(r1, i);
                    remove_column(r2, i);
                    screened_SNP_original_indices.erase(screened_SNP_original_indices.begin()+i);
                }
            }
            printf("Step 3: Time taken: %.5fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
            exit(-1);
            break;
        }

        if (!loop_break_indicator) {
            previous_R2_cohort1 = R2_cohort1;
            previous_R2_cohort2 = R2_cohort2;
            R1_inv_pre = R1_inv_post;
            R2_inv_pre = R2_inv_post;

            output_b_cohort1 = beta1;
            output_b_cohort2 = beta2;
            output_se2_cohort1 = beta_var1;
            output_se2_cohort2 = beta_var2; 
        }

        cout << "iter " << ++iter_num << " finished" << endl;
        cout << "--------------------------------" << endl;
    }

    inverse_var_meta(output_b_cohort1, output_b_cohort2, output_se2_cohort1, output_se2_cohort2, sumstat_merge);
    cout << "bJ.ma: " << sumstat_merge.col(0).transpose() << endl;
    cout << "seJ.ma: " << sqrt(sumstat_merge.col(1)).transpose() << endl;
    cout << "zJ: " << (sumstat_merge.col(0)/sqrt(sumstat_merge.col(1))).transpose() << endl;
    cout << "pJ: " << scientific << sumstat_merge.col(3).transpose() << endl;
}
