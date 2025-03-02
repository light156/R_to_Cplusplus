#include "trans_ancestry_cojo.h"
#include <time.h>
#include <numeric>


int main() {
    /*
    clock_t tStart = clock();

    vector<int> v(100000) ; // vector with 100 ints.
    iota (begin(v), end(v), 0);

    for (int i=0;i<5000;i++) {
        v.erase(v.begin()+i*8);
    }

    printf("Time taken: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
    return 0;
    */

    string cojoFile1, cojoFile2, PLINK1, PLINK2;
    cojoFile1 = "Chr22_LD_complexity_Simple_Casual_SNPs_5_Scenario_9_AFR.Pheno.PHENO1.glm.linear_gcta_format";
    cojoFile2 = "Chr22_LD_complexity_Simple_Casual_SNPs_5_Scenario_9_EUR.Pheno.PHENO1.glm.linear_gcta_format";
    PLINK1 = "test";
    PLINK2 = "test";

    taCOJO g = taCOJO();
    g.read_files(cojoFile1, cojoFile2, PLINK1, PLINK2);
    g.main_loop();

    return 0;
}