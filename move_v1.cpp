#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <vector>
#include <random>

using namespace std;

/* DATI PLATINO
*passo reticolare: 3,9242 ang
distanza tra primi vicini: 2.7748 ang

QUESTA VERSIONE FA MUOVERE UN ATOMO SUL RETICOLO FUORI.
VERSIONE FUNZIONANTE CHE NONN METTE ALCUNA BARRIERA NEL MOTO, IL MONTECARLO È FATTO
ACCETTANDO SEMPRE LE MOSSE

ATTENZIONE I SITI SONO RUOTATI RISPETTO ALLA COSTRUZIONE INIZIALE
*/

vector <int> pv(int n_siti,int posizione,double** siti,double d,double tolerance){
    vector <int> p_v;
    for(int i=0;i<n_siti;i++){
        if(i!=posizione){
            if(sqrt(pow((siti[posizione][0]-siti[i][0]),2)+pow((siti[posizione][1]-siti[i][1]),2)+pow((siti[posizione][2]-siti[i][2]),2))-d<tolerance){
                p_v.push_back(i);
            }
        }
    }

    return p_v;
}

int main(){
    ifstream ifile("dati_piramide.txt");
    double d_pv;
    int n_l;
    int n_cut[3]; 
    ifile>>d_pv>>n_l;
    for(int i=0;i<3;i++){
        ifile>>n_cut[i];
    }

    //registro in un vettore tutti i siti possibili
    ifstream file_siti("siti.xyz");
    int n_siti;
    file_siti>>n_siti;
    double** siti=new double*[n_siti];
    for(int i=0;i<n_siti;i++){
        siti[i]=new double[3];
        for(int k=0;k<3;k++){
            file_siti>>siti[i][k];
        }
    }
    

    ifstream base("piramide.xyz");
    int n_atoms;
    base>>n_atoms;
    string line;
    vector <string> lines;
    while (std::getline(base, line)) {
        lines.push_back(line);
    }

    //random
    random_device rd;
    mt19937 gen(rd());

    double* adatom=new double[3];
    //estraggo la posizione iniziale
    uniform_int_distribution<> distribuzione_posizione(0,n_siti-1);
    int initial= distribuzione_posizione(gen);
    cout<<initial<<endl;
    for(int l=0;l<3;l++){
        adatom[l]=siti[initial][l];
        cout<<adatom[l]<<endl;
    }


    double tolerance=0.0005;
    int mc_step=1000;
    int posizione=initial;
    for(int i=0;i<mc_step;i++){
        vector<int> pv1=pv(n_siti,posizione,siti,d_pv,tolerance);

        //for(int q=0;q<pv1.size();q++) cout<<pv1[q]<<endl;

        uniform_int_distribution<> distribuzione_pv(0,pv1.size()-1);
        //cout<<"test"<<endl;
        int mossa= distribuzione_pv(gen);
        for(int l=0;l<3;l++) adatom[l]=siti[pv1[mossa]][l];
        posizione=pv1[mossa];

        ostringstream filename;
        string num= to_string(i+1);
        unsigned int number_of_zeros = 7 - num.length();
        num.insert(0,number_of_zeros,'0');
        filename<<"singolo"<<num<<".xyz";
        ofstream ofile;
        ofile.open(filename.str());
        ofile<<n_atoms+1;

        for(int j=0;j<lines.size();j++){
            ofile<<lines[j]<<endl;
        }
        ofile<<"Au "<<adatom[0]<<" "<<adatom[1]<<" "<<adatom[2]<<endl;
    }
}