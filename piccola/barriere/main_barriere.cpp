#include "Reticolo.h"
#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>

using namespace std;

int main(){
    Reticolo piramide;
    Reticolo reticolo;

    double start[3];
    double finish[3];
    vector<int> fermi;

    ifstream posizioni("start.xyz");
    piramide.lettura_file(posizioni);

    ifstream file_siti("siti_jmol.xyz");
    reticolo.lettura_file(file_siti);

    //ATTENZIONE: NEL FILE DI INPUT LE POSIZIONI FORNITE SONO QUELLE CON GLI INDICI DI JMOL CHE PARTE DA 1
    //MENTRE IL VETTORE DEI SITI PARTE DA 0 QUINDI DEVO RIDURRRE TUTTO DI 1 
    ifstream ifile("pos_barriere.txt");
    ifile>>start[0]>>start[1]>>start[2];
    cout<<start[0]<<" "<<start[1]<<" "<<start[2]<<endl;
    ifile>>finish[0]>>finish[1]>>finish[2];
    cout<<finish[0]<<" "<<finish[1]<<" "<<finish[2]<<endl;
    int i;
    while(ifile>>i) fermi.push_back(i-1);

    ofstream testfile;
    testfile.open("test.xyz");
    testfile<<piramide.N<<endl;
    testfile<<" Pt Pt 1 1"<<endl;

    //traslo nella posizione dell'atomo di partenza 
    for(int i=0;i<piramide.N;i++){
        for(int j=0;j<3;j++){
            piramide.siti[i][j]=piramide.siti[i][j] - start[j];
        }
    }
    for(int i=0;i<reticolo.N;i++){
        for(int j=0;j<3;j++){
            reticolo.siti[i][j] = reticolo.siti[i][j] -start[j];
        }
    }
    for(int i=0;i<3;i++){
        finish[i]=finish[i]-start[i];
        start[i]=0;
    }
    //calcolo modulo distanza
    double mod;
    for(int i=0;i<3;i++){
        mod += finish[i]*finish[i];
    }
    mod=sqrt(mod);
    //divido e ricavo angoli
    vector <double> v;
    for (int i=0;i<3;i++) v.push_back(finish[i]/mod);
    cout<<"congiungente "<<v[0]<<" "<<v[1]<<" "<<v[2]<<endl;

    //SE LA X DEL VETTORE CONGIUNGENTE è NEGATIVA CAMBIO IL SEGNO A X E Y CHE è COME RUOTARE DI 180
    if(v[0]<0){
        cout<<"pino"<<endl;
        for(int j=0;j<2;j++){
            for(int i=0;i<piramide.N;i++) piramide.siti[i][j]= (-1.)*piramide.siti[i][j];
            for(int i=0;i<reticolo.N;i++) reticolo.siti[i][j] = (-1.)*reticolo.siti[i][j];
            v[j]=(-1.)*v[j];
        }
    }

    double theta,phi,alpha;
    theta = atan(v[1]/v[0]);
    phi = acos(v[2]);
    cout<<"theta "<<theta<<"  phi "<<phi<<endl;
    cout<<"prova "<<cos(theta)<<" "<<sin(theta)<<endl;

    //ruoto prima attorno a z in senso antiorario di theta
    //NOTA BENE: IO VOGLIO RUOTARE GLI ASSI DI THETA QUINDI TUTTO IL SISTEMA SARà RUOTATO IN SENSO OPPOSTO, OVVERO IN SENSO ORARIO
    for(int i=0;i<piramide.N;i++){
        double provvisorio[3];
        provvisorio[0]= piramide.siti[i][0]*cos(theta) + piramide.siti[i][1]*sin(theta);
        provvisorio[1]= -piramide.siti[i][0]*sin(theta) + piramide.siti[i][1]*cos(theta);
        provvisorio[2]= piramide.siti[i][2];  
        for(int j=0;j<3;j++)    piramide.siti[i][j]=provvisorio[j]; 
    }
    for(int i=0;i<reticolo.N;i++){
        double provvisorio[3];
        provvisorio[0]= reticolo.siti[i][0]*cos(theta) + reticolo.siti[i][1]*sin(theta);
        provvisorio[1]= -reticolo.siti[i][0]*sin(theta) + reticolo.siti[i][1]*cos(theta);
        provvisorio[2]= reticolo.siti[i][2];  
        for(int j=0;j<3;j++) reticolo.siti[i][j]=provvisorio[j];
    }

    double provvisorio[3];
    provvisorio[0] = v[0]*cos(theta) + v[1]*sin(theta);
    provvisorio[1] = -v[0]*sin(theta) + v[1]*cos(theta);
    provvisorio[2] = v[2];
    for(int j=0;j<3;j++) v[j]=provvisorio[j];

    //ruoto gli assi attorno a y in senso orario di 90-phi, a meno che non siamo già lungo il piano xy
    //ovvero ruoto il sistema in senso antiorario di 90-phi
    double b=1.;
    double a= (M_PI/2.)-phi;
    cout<<"a "<<a<<" "<<sin(a)<<endl;
    if(v[2]<0) b=-1.;
    
    for(int i=0;i<piramide.N;i++){
        double provvisorio[3];
        provvisorio[0]= piramide.siti[i][0]*cos(a) - b*piramide.siti[i][2]*sin(a);
        provvisorio[1]= piramide.siti[i][1];
        provvisorio[2]= + b* piramide.siti[i][0]*sin(a) + piramide.siti[i][2]*cos(a); 
        for(int j=0;j<3;j++) piramide.siti[i][j]=provvisorio[j];
        testfile<<"Pt "<<piramide.siti[i][0]<<" "<<piramide.siti[i][1]<<" "<<piramide.siti[i][2]<<endl; 
    }
    for(int i=0;i<reticolo.N;i++){
        double provvisorio[3];
        provvisorio[0]= reticolo.siti[i][0]*cos(a) - b*reticolo.siti[i][2]*sin(a);
        provvisorio[1]= reticolo.siti[i][1];
        provvisorio[2]= + b*reticolo.siti[i][0]*sin(a) + reticolo.siti[i][2]*cos(a); 
        for(int j=0;j<3;j++) reticolo.siti[i][j]=provvisorio[j];
    }
    double provv[3];
    provv[0] = v[0]*cos(a) - b* v[2]*sin(a);
    provv[1] = v[1];
    provv[2] = + b* v[0]*sin(a) + v[2]*cos(a);
    for(int j=0;j<3;j++) v[j]=provv[j];
    

    for(int i=0;i<3;i++){
        finish[i]=v[i]*mod;
        cout<<finish[i]<<" ";
    }
    cout<<endl;

    ofstream ofile;
    ofile.open("intermedie.txt");

    double p=1/100.;
    //if(m<0){
        for(int y=1;y<101;y++){
            ofile<< finish[0]*p*y<<" "<<0<<" "<<0<<endl;
        }
    //}
    /*
    //calcolo i vettori direttori, moltiplicando le componenti per le parti infinitesime ho vettori uguali ma più piccoli che saranno le mie posizioni
    else if(m>0){
        vector<double> intermedio=reticolo.siti[m-1];
        vector<double> r1 ;
        for(int i=0;i<3;i++) r1.push_back(intermedio[i]-start[i]);
        for(int y=1;y<51;y++){
            ofile<< r1[0]*p*(y*2.)<<" "<<r1[1]*p*(y*2.)<<" "<<r1[2]*p*(y*2.)<<endl;
        }
        vector<double> r2 ;
        for(int i=0;i<3;i++) r2.push_back(finish[i]-intermedio[i]);
        for(int y=1;y<51;y++){
            ofile<< r2[0]*p*(y*2.) +intermedio[0]<<" "<<r2[1]*p*(y*2.) +intermedio[1]<<" "<<r2[2]*p*(y*2.) +intermedio[2]<<endl;
        }
    }
    */
    ofstream ofile1;
    ofile1.open("input.xyz");
    piramide.scrittura_preferenziale_file(piramide.N,fermi,ofile1);

}