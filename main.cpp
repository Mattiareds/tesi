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
    ifstream dati_moto("dati_moto.txt");
    reticolo.lettura_dati_moto(dati_moto);
    
    //verifica piani
    reticolo.costruttore_piani();
    reticolo.associa_piano();


    //moto, la posizione è accettata solo se non è su uno spigolo
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> distribuzione_posizione(0,reticolo.N-1);
    bool start=false;
    int posizione;
    while(start==false){
        int initial=distribuzione_posizione(gen);
        if(reticolo.siti[initial][7]==1){
            start= true;
            reticolo.seleziona_sito(initial);
        }else{
            start= false;
        }
        posizione=initial;
    }

    
    cout<<reticolo.mc_step;
    for(int i=0;i<reticolo.mc_step;i++){
        vector<int> pv=reticolo.pv(posizione);

        //seleziono il primo vicino
        uniform_int_distribution<> distribuzione_pv(0,pv.size()-1);
        int mossa= distribuzione_pv(gen);

        uniform_real_distribution<> r(0.,1.);
        double g=r(gen);
        posizione=reticolo.metropolis(posizione,pv[mossa],g);//nuova posizione, in funzione di quello che esce fuori dal metropolis
        //cout<<"marco "<<posizione<<endl;
        reticolo.seleziona_sito(posizione);
        reticolo.output_moto(i,piramide);

    }
}

