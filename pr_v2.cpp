#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;

/* DATI PLATINO
*passo reticolare: 3,9242 ang
distanza tra primi vicini: 2.7748 ang

QUESTA VERSIONE FA UN OTTAEDRO TRONCATO E BASTA, SERVE COME BACKUP
*/

vector<vector<double>> square_maker(int l,double* c, double d_pv){
    vector<vector<double>> square;

    //creo nell'origine e poi lo traslo
    for(int i=0;i<l;i++){
        for(int j=0;j<l;j++){
            vector<double> coord={0,d_pv*j,d_pv*i};
            square.push_back(coord);
        }
    }

    //traslo nell'origine
    vector <double> cm={0,0,0};

    for(int t=0;t<(l*l);t++){
        cm[0]+=(square[t][0]/double(l*l));
        cm[1]+=(square[t][1]/double(l*l));
        cm[2]+=(square[t][2]/double(l*l));
    }
    for(int i=0;i<(l*l);i++){
        square[i][0]=square[i][0]-cm[0];
        square[i][1]=square[i][1]-cm[1];
        square[i][2]=square[i][2]-cm[2];        
    }

    for(int i=0;i<(l*l);i++){
        square[i][0]=square[i][0]+c[0];
        square[i][1]=square[i][1]+c[1];
        square[i][2]=square[i][2]+c[2];        
    }
    return square;
}

vector<vector<double>> cutted_square_maker(int l,double* c, double d_pv,int n_cut){
    vector<vector<double>> square=square_maker(l,c,d_pv);
    vector<vector<double>> cutted;
    for(int i=0;i<square.size();i++){
        if(((fabs(square[i][1])+fabs(square[i][2]))<(double(l-n_cut)*(d_pv-(0.05*d_pv))))){ //traslazione linea di una frazione del passo per sicurezza, senno a causa del troncamento dava problemi
            cutted.push_back(square[i]);
        }
    }
    return cutted;
}

int main(){

    //DATIIIII
    ifstream ifile("dati_piramide.txt");
    double d_pv;
    int n_l;
    int n_cut[3]; 
    ifile>>d_pv>>n_l;
    for(int i=0;i<3;i++){
        ifile>>n_cut[i];
    }
    
    //setto il vettore centro del quadrato che verrà posizionato
    double* c=new double[3];
    for (int i=0;i<3;i++) c[i]=0;

    vector<vector<double>> atoms;

    //creo il primo centrato nell'origine
    vector<vector<double>> q=cutted_square_maker(n_l,c,d_pv,n_cut[2]);

    for(int i=0;i<q.size();i++){
        atoms.push_back(q[i]);
    }

    //creo gli altri nel semiasse positivo
    double h= d_pv*sin(((54.74)*M_PI)/180.);
    int j=1;
    while(j<(n_l-n_cut[0])){
        int l=n_l-j;
        c[0]=j*h;

        //metto quadrati tagliati finché ha senso
        if(l>(n_l-n_cut[2])){
            vector<vector<double>> bmb=cutted_square_maker(l,c,d_pv,n_cut[2]-j);
            for(int i=0;i<bmb.size();i++) atoms.push_back(bmb[i]);
            j++;
        }else{
            vector<vector<double>> bmb=square_maker(l,c,d_pv);
            for(int i=0;i<bmb.size();i++) atoms.push_back(bmb[i]);
            j++;
        }
    }

    //creo asse negativo
    int y=1;
    while(y<(n_l-n_cut[1])){
        int q=n_l-y;
        c[0]=(-1)*y*h;

        //metto quadrati tagliati 
        if(q>(n_l-n_cut[2])){
            vector<vector<double>> bmb=cutted_square_maker(q,c,d_pv,n_cut[2]-y);
            for(int i=0;i<bmb.size();i++) atoms.push_back(bmb[i]);
            y++;
        }else{
            vector<vector<double>> bmb=square_maker(q,c,d_pv);
            for(int i=0;i<bmb.size();i++) atoms.push_back(bmb[i]);
            y++;
        }
    }
    cout<<atoms.size()<<endl;

    ofstream ofile;
    ofile.open("piramide.xyz");
    ofile<<atoms.size()<<endl;
    ofile<<"Pt "<<"Pt "<<"1 "<<"1"<<endl;
    
    for(int i=0;i<atoms.size();i++){
        //cout<<i<<endl;
        ofile<<"Pt "<<atoms[i][0]<<" "<<atoms[i][1]<<" "<<atoms[i][2]<<endl;
    }




}
