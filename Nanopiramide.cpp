#include <vector>
#include <fstream>
#include "Nanopiramide.h"

void Nanopiramide::lettura_piramide(ifstream ifile){
    int n;
    ifile>>n;
    double** siti=new double*[n_siti];
    for(int i=0;i<n_siti;i++){
        siti[i]=new double[3];
        for(int k=0;k<3;k++){
            file_siti>>siti[i][k];
        }
    }
}