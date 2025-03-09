#include "trans_ancestry_cojo.h"


int main() {

    cout << fixed << setprecision(15) << endl;

    clock_t tStart = clock();
    
    string cojoFile1, cojoFile2, PLINK1, PLINK2;
    /*
    cojoFile1 = "test_data_old/Chr22_LD_complexity_Simple_Casual_SNPs_5_Scenario_9_AFR.Pheno.PHENO1.glm.linear_gcta_format";
    cojoFile2 = "test_data_old/Chr22_LD_complexity_Simple_Casual_SNPs_5_Scenario_9_EUR.Pheno.PHENO1.glm.linear_gcta_format";
    PLINK1 = "test_data_old/test";
    PLINK2 = "test_data_old/test";
    */
    cojoFile1 = "Height_AFR_GIANT_QCed_withoutUKB_Chr22.sumstat";
    cojoFile2 = "Height_EUR_GIANT_QCed_withoutUKB_Chr22.sumstat";
    PLINK1 = "1000G_Chr22_AFR_QCed";
    PLINK2 = "1000G_Chr22_EUR_QCed";

    taCOJO g = taCOJO();
    g.initialize_hyperparameters();
    g.read_files(cojoFile1, cojoFile2, PLINK1, PLINK2);
    g.main_loop();
    
    printf("Time taken: %.5fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    return 0;
}