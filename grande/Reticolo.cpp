#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
#include "Reticolo.h"
#include <random>
#include <string>
#include <sstream>

using namespace std;

//lettura del classico file .xyz di jmol
void Reticolo::lettura_file(ifstream& ifile){
    ifile>>N;
    ifile>>specie_chimica>>specie_chimica>>a>>b;
    siti.resize(N, std::vector<double>(8,0.0));
    for(int i=0;i<N;i++){
        ifile>>specie_chimica>>siti[i][0]>>siti[i][1]>>siti[i][2];
        //assegno -1 agli altri indici così non mi incasino coi piani
        for(int j=3;j<8;j++){
            siti[i][j]=-1.;
        }
    }
    
}

void Reticolo::lettura_dati_piramide(ifstream& ifile){
    ifile>>d_pv>>spigolo;
    for(int i=0;i<3;i++){
        ifile>>n_cut[i];
    }
    ifile>>mc_step;
}

//prende in ingresso il numero di atomi che dipende da caso a caso
void Reticolo::scrittura_file(int number_of_atoms,ofstream& ofile){
    ofile<<number_of_atoms<<endl;
    ofile<<specie_chimica<<" "<<specie_chimica<<" "<<a<<" "<<b<<endl;
    for(int i=0;i<siti.size();i++){
        ofile<<specie_chimica<<" "<<siti[i][0]<<" "<<siti[i][1]<<" "<<siti[i][2]<<endl;
    }
}

//questo mi scrive con certi atomi nelle prime posizioni e gli altri sotto, perché serve per tenerli fermi durante il drag method
void Reticolo::scrittura_preferenziale_file(int number_of_atoms,vector <int> fermi, std::ofstream& ofile){
    ofile<<number_of_atoms<<endl;
    ofile<<specie_chimica<<" "<<specie_chimica<<" "<<a<<" "<<b<<endl;
    for(int i=0;i<fermi.size();i++){
        ofile<<specie_chimica<<" "<<siti[fermi[i]][0]<<" "<<siti[fermi[i]][1]<<" "<<siti[fermi[i]][2]<<endl; 
    }
    for(int i=0;i<siti.size();i++){
        bool test=false;
        for(int j=0;j<fermi.size();j++){
            if(i==fermi[j]) test=true;
        }
        if(test==false) ofile<<specie_chimica<<" "<<siti[i][0]<<" "<<siti[i][1]<<" "<<siti[i][2]<<endl;
    }
}

void Reticolo::aggiunte_file(ofstream& ofile){
    for(int i=0;i<siti.size();i++){
        ofile<<specie_chimica<<" "<<siti[i][0]<<" "<<siti[i][1]<<" "<<siti[i][2]<<endl;
    }
}

void Reticolo::seleziona_sito(int i){
    sito=siti[i];
}
void Reticolo::aggiunte_atomo(ofstream& ofile){

    ofile<<specie_chimica<<" "<<sito[0]<<" "<<sito[1]<<" "<<sito[2]<<endl;
}

void Reticolo::costruttore_piani(){
    //100 tipo 0
    piani.push_back({1,0,0,d_pv/(sqrt(2))*(spigolo-1-n_cut[2])});
    piani.push_back({-1,0,0,d_pv/(sqrt(2))*(spigolo-1-n_cut[2])});
    piani.push_back({0,1,0,d_pv/(sqrt(2))*(spigolo-1-n_cut[2])});
    piani.push_back({0,-1,0,d_pv/(sqrt(2))*(spigolo-1-n_cut[2])});
    piani.push_back({0,0,1,d_pv/(sqrt(2))*(spigolo-1-n_cut[0])});
    piani.push_back({0,0,-1,d_pv/(sqrt(2))*(spigolo-1-n_cut[1])});
    //111 tipo 1
    piani.push_back({1,1,1,d_pv/(sqrt(2))*(spigolo-1)});
    piani.push_back({-1,1,1,d_pv/(sqrt(2))*(spigolo-1)});
    piani.push_back({-1,-1,1,d_pv/(sqrt(2))*(spigolo-1)});
    piani.push_back({1,-1,1,d_pv/(sqrt(2))*(spigolo-1)});
    piani.push_back({-1,-1,-1,d_pv/(sqrt(2))*(spigolo-1)});
    piani.push_back({-1,1,-1,d_pv/(sqrt(2))*(spigolo-1)});
    piani.push_back({1,-1,-1,d_pv/(sqrt(2))*(spigolo-1)});
    piani.push_back({1,1,-1,d_pv/(sqrt(2))*(spigolo-1)});
}

bool Reticolo::test_piano(int i,int j){
    if(abs(siti[i][0]*piani[j][0]+siti[i][1]*piani[j][1]+siti[i][2]*piani[j][2]-piani[j][3])<0.05){
        return true;
    }
    return false;
}

void Reticolo::associa_piano(){
    int count[siti.size()] = { };
    for(int i=0;i<siti.size();i++){
        for(int j=0;j<piani.size();j++){
            if(test_piano(i,j)==true) {
                if(count[i]==0){
                    siti[i][3]=(double) j;
                }
                else if(count[i]==1){
                    siti[i][4]=(double) j;
                }
                else if(count[i]==2){
                    siti[i][5]=(double) j;
                }
                else if(count[i]==3){
                    siti[i][6]=(double) j;
                }
                count[i]++;
            }
        }
        siti[i][7]=count[i];//so quanti piani
    }
}
int plane_type(int i){
    if(i<6) return 0;
    return 1;
}

vector <int> Reticolo::pv(int pos){
    vector <int> p_v;
    for(int i=0;i<siti.size();i++){
        if(i!=pos){
            if(sqrt(pow((siti[pos][0]-siti[i][0]),2)+pow((siti[pos][1]-siti[i][1]),2)+pow((siti[pos][2]-siti[i][2]),2))-d_pv<0.05){
                p_v.push_back(i);
            }
        }
    }
    return p_v;
}

void Reticolo::dati_guscio(){
    spigolo=spigolo+2;
    for(int i=0;i<3;i++) n_cut[i]=n_cut[i]+1;
}

void Reticolo::dati_secondo_guscio(){
    spigolo=spigolo+4;
    for(int i=0;i<3;i++) n_cut[i]=n_cut[i]+2;
}

ofstream Reticolo::output_moto(int i, Reticolo& piramide){
    ostringstream filename;
    string num=to_string(i+1);
    unsigned int number_of_zeros= 7 - num.length();
    num.insert(0,number_of_zeros,'0');
    filename<<"singolo"<<num<<".xyz";
    ofstream ofile;
    ofile.open(filename.str());
    ofile<<piramide.N+1<<endl;
    ofile<<piramide.specie_chimica<<" "<<piramide.specie_chimica<<" "<<1<<" "<<1<<endl;

    piramide.aggiunte_file(ofile);
    this->aggiunte_atomo(ofile);

    return ofile;
}

void Reticolo::lettura_dati_moto(ifstream& ifile){
    for(int i=0;i<5;i++){
        ifile>>barriere[i];
    }
    ifile>>metodo;
}

//funzione easy che ritorna se fare o no la mossa
int confronto(int pos, int mossa, double random, double b){
    if(b>0 && random > exp(-b)){
        return pos;
    }
    return mossa;
}

int Reticolo::metropolis(int pos,int mossa,double random){
    //caso 1: mossa sempre accettata:
    if(metodo==1) return mossa;

    //caso 2: 5 barriere a seconda del moto
    //NON POSSO FINIRE SUGLI SPIGOLI, NON LO CONSIDERO COME PASSO NEL MC STEP
    //unica opzione è appartenere ad un piano, quindi cambia la casistica
    if(metodo==2){
        //se arrivo su un piano solo vuol dire che mi sto muovendo sul piano
        if(siti[mossa][7]==1){
            current_plane=siti[pos][3];
            return confronto(pos,mossa,random,barriere[1 - plane_type((int) siti[pos][3])]);
        }
        //se arrivo su spigolo
        else if(siti[mossa][7]!=1){
            vector<int> ppv=pv(mossa);
            double d_min=N*d_pv;
            for(int l=0;l<ppv.size();l++){
                if(siti[ppv[l]][7]==1){ //solo se non è sullo spigolo
                    if(siti[ppv[l]][3]!=current_plane){ //solo se non è sullo stesso piano
                        double d= sqrt(pow((siti[ppv[l]][0]-siti[pos][0]),2)+pow((siti[ppv[l]][1]-siti[pos][1]),2)+pow((siti[ppv[l]][2]-siti[pos][2]),2));
                        if(d<d_min) {
                            d_min=d;
                            mossa=ppv[l];
                        }
                    }
                }
            }
            //ora che ho determinato il piu vicino decido
            int barriera;
            int jump= plane_type(current_plane)-plane_type((int) siti[mossa][3]);
            if(jump==0) barriera=2; 
            else if(jump>0) barriera=3;
            else if(jump<1) barriera=4;
            int result= confronto(pos,mossa,random,barriere[barriera]);
            current_plane=siti[result][3];
            return result;
        }
        
    }
    return pos;
}

/*
int Reticolo::metropolis(int pos,int mossa,double random){
    //caso 1: mossa sempre accettata:
    if(metodo==1) return mossa;
    //caso 2: 5 barriere a seconda del moto
    if(metodo==2){
        //se appartengono ad un piano solo vado sicuramente nello stesso
        if(siti[pos][7]==1){
            int t=plane_type((int) siti[pos][3]);
            if(t==0) {
                current_plane=siti[pos][3];
                return confronto(pos, mossa,random, barriere[1]);
            }
            else if(t==1) {
                current_plane=siti[pos][3];
              return confronto(pos,mossa,random,barriere[0]);
            }

        }
        //se finisco sullo spigolo/vertice mi muovo sicuramente altrove
        else{
            //se un solo piano vuol dire che ho cambiato piano e quindi vedo se applicare la barriera
            if(siti[mossa][7]==1){
                int jump=plane_type(current_plane)-plane_type((int) siti[mossa][3]);
                if(current_plane==siti[mossa][3]) return confronto(pos,mossa,random,barriere[1-plane_type((int) siti[mossa][3])]);
                else if(jump==0){
                    int result=confronto(pos,mossa,random,barriere[2]);
                    if(result==mossa) current_plane=siti[mossa][3]; //cambia piano se ci spostiamo
                    return result;
                }
                else if(jump>0){
                    int result=confronto(pos,mossa,random,barriere[3]);
                    if(result==mossa) current_plane=siti[mossa][3];
                    return result;
                }
                else if(jump<0){
                    int result=confronto(pos,mossa,random,barriere[4]);
                    if(result==mossa) current_plane=siti[mossa][3];
                    return result;
                }
            }
            //essendo una posizione instabile non ha senso che viaggi lungo lo spigolo
            else{
                return pos;
            }
            
        }
    }
    return pos;
}

*/