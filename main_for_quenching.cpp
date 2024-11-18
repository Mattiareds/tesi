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

    //verifica piani
    reticolo.costruttore_piani();
    reticolo.associa_piano();

    //posiziono casualmente, la posizione è accettata solo se non è su uno spigolo
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> distribuzione_posizione(0,reticolo.N-1);
    bool start=false;
    while(start==false){
        int initial=distribuzione_posizione(gen);
        if(reticolo.siti[initial][7]==1){
            start= true;
            reticolo.seleziona_sito(initial);
        }else{
            start= false;
        }
    }

    ofstream ofile;
    ofile.open("input.xyz");
    piramide.scrittura_file(piramide.N +1,ofile);
    reticolo.aggiunte_atomo(ofile);
}

