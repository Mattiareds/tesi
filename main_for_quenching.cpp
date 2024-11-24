#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <random>
#include "Reticolo.h"

using namespace std;

int main(){

    Reticolo piramide;
    Reticolo reticolo;

    ifstream ifile("dati_piramide.txt");
    reticolo.lettura_dati_piramide(ifile);
    reticolo.dati_guscio();
    ifstream file_siti("siti_jmol.xyz");
    reticolo.lettura_file(file_siti);
    ifstream base("piramide.xyz");
    piramide.lettura_file(base);

    ifile.close();
    file_siti.close();
    base.close();

    //verifica piani
    reticolo.costruttore_piani();
    reticolo.associa_piano();


    vector<int> indici;

    for(int i=0;i<reticolo.N;i++){
        cout<<reticolo.siti[i][7]<<endl;
        if(reticolo.siti[i][7]==1){
            cout<<"pino"<<endl;
            if(reticolo.siti[i][3]==5){
                indici.push_back(i);
                
            }
        }
    }

    ofstream ofile;
    ofile.open("input.xyz");
    piramide.scrittura_file(piramide.N +indici.size(),ofile);
    for(int j=0;j<indici.size();j++){
        reticolo.seleziona_sito(indici[j]);
        reticolo.aggiunte_atomo(ofile);
    }

}

