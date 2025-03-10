#include "trans_ancestry_cojo.h"


int main() {

    cout << fixed << setprecision(15) << endl;

    clock_t tStart = clock();
    
    string cojoFile1, cojoFile2, PLINK1, PLINK2;

    // cojoFile1 = "Height_AFR_GIANT_QCed_withoutUKB_Chr22.sumstat";
    // cojoFile2 = "Height_EUR_GIANT_QCed_withoutUKB_Chr22.sumstat";

    int num = 40;
    cojoFile1 = "Simstat/Chr22_LD_Block_1_Casual_SNPs_20_Scenario_"+to_string(num)+"_AFR.PHENO1.glm.linear_gcta_format";
    cojoFile2 = "Simstat/Chr22_LD_Block_1_Casual_SNPs_20_Scenario_"+to_string(num)+"_EUR1.PHENO1.glm.linear_gcta_format";
    
    PLINK1 = "1000G_Chr22_AFR_QCed";
    PLINK2 = "1000G_Chr22_EUR_QCed";

    taCOJO g = taCOJO();
    g.initialize_hyperparameters();
    g.read_files(cojoFile1, cojoFile2, PLINK1, PLINK2);
    g.main_loop();      
    g.save_results("results.jma.cojo");
    
    printf("Time taken: %.5fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    return 0;
}