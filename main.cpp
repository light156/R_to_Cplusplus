#include "trans_ancestry_cojo.h"
#include <time.h>
#include <numeric>


int main() {

    cout << fixed << setprecision(15) << endl;

    string cojoFile1, cojoFile2, PLINK1, PLINK2;
    /*
    cojoFile1 = "Chr22_LD_complexity_Simple_Casual_SNPs_5_Scenario_9_AFR.Pheno.PHENO1.glm.linear_gcta_format";
    cojoFile2 = "Chr22_LD_complexity_Simple_Casual_SNPs_5_Scenario_9_EUR.Pheno.PHENO1.glm.linear_gcta_format";
    PLINK1 = "test";
    PLINK2 = "test";
    */
    cojoFile1 = "Height_AFR_GIANT_QCed_withoutUKB_Chr22.sumstat";
    cojoFile2 = "Height_EUR_GIANT_QCed_withoutUKB_Chr22.sumstat";
    PLINK1 = "1000G_Chr22_AFR_QCed";
    PLINK2 = "1000G_Chr22_EUR_QCed";

    taCOJO g = taCOJO();
    g.read_files(cojoFile1, cojoFile2, PLINK1, PLINK2);
    g.main_loop();

    return 0;
}