#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

/* DATI PLATINO
*passo reticolare: 3,9242 ang
distanza tra primi vicini: 2.7748 ang

QUESTA VERSIONE FA UN OTTAEDRO TRONCATO E RESTITUISCE UN SECONDO FILE CON I SITI SUI QUALI SI MUOVE L'ADATOM
RISPETTO ALLA 3 SI CREA CON ASSE Z NON X
*/

vector<vector<double>> square_maker(int l,double* c, double d_pv){
    vector<vector<double>> square;

    //creo nell'origine e poi lo traslo
    for(int i=0;i<l;i++){
        for(int j=0;j<l;j++){
            vector<double> coord={d_pv*j,d_pv*i,0};
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
        if(((fabs(square[i][1])+fabs(square[i][0]))<(double(l-n_cut)*(d_pv-(0.05*d_pv))))){ //traslazione linea di una frazione del passo per sicurezza, senno a causa del troncamento dava problemi
            cutted.push_back(square[i]);
        }
    }
    return cutted;
}

bool equal_rows(vector<double> a,vector<double> b,double tolerance){
    for(int i=0;i<a.size();i++){
        if(abs(a[i]-b[i])>tolerance){
            return false;
        }
    }
    return true;
}

vector<vector<double>> rotate_xy_45deg(vector<vector<double>> v){
    vector<vector<double>> ruotato(v.size(),vector<double>(3,0.0));
    double a=sqrt(2)/2.;

    for(int i=0;i<v.size();i++){
        ruotato[i][0]= a*v[i][0]-a*v[i][1];
        ruotato[i][1]=a*v[i][0]+a*v[i][1];
        ruotato[i][2]=v[i][2];
    }
    return ruotato;
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

    //per creare il bordo faccio girare il programma due volte, alla seconda aggiungo due 
    //alla fine per il bordo confronto i due risultati e scarto tutti gli atomi in comune
    
    //vettori dove vado a salvare tutto
    vector<vector<double>> piccolo;
    vector<vector<double>> grande;

    for(int ciclo=0;ciclo<2;ciclo++){

        if(ciclo==1) {
            n_l=n_l+2;
            for(int p=0;p<3;p++) n_cut[p]=n_cut[p]+1;
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

        //creo gli altri nel semiasse positivo STAVOLTA DELLE Z NON DELLE X
        double h= d_pv*sin(M_PI/4.);
        int j=1;
        while(j<(n_l-n_cut[0])){
            int l=n_l-j;
            c[2]=j*h;

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
            c[2]=(-1)*y*h;

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

        if(ciclo==0) piccolo=atoms;
        if(ciclo==1) grande=atoms;
    }
    

    //mi salvo gli indici da eliminare
    vector<int> index;
    for(int j=0;j<grande.size();j++){
        for(int i=0;i<piccolo.size();i++){
            if(equal_rows(grande[j],piccolo[i],0.001)==true){
                index.push_back(j);
                cout<<j<<endl;
                cout<<i<<endl;
                cout<<grande[j][0]<<" "<<grande[j][1]<<" "<<grande[j][2]<<endl;
                cout<<piccolo[i][0]<<" "<<piccolo[i][1]<<" "<<piccolo[i][2]<<endl;
            }else{
                //cout<<grande[j][0]<<" "<<grande[j][1]<<" "<<grande[j][2]<<endl;
                //cout<<piccolo[i][0]<<" "<<piccolo[i][1]<<" "<<piccolo[i][2]<<endl;
            }
        }
    }

    //ordino gli indici in modo descrescente perché conviene
    //riumuoverli dal più alto sennò mi cambia l'indice che mi interessa durante l'operazione
    sort(index.rbegin(),index.rend());

    vector<vector<double>>::iterator it = grande.begin();
    for(int i=0;i<index.size();i++){
        grande.erase(it+index[i]);
    }


    //PRIMA DI STAMPARE RUOTO DI 45 GRADI NEL PIANO XY COSI VIENE CON GLI ANGOLI GIUSTI
    vector<vector<double>> r_grande=rotate_xy_45deg(grande);
    vector<vector<double>> r_piccolo=rotate_xy_45deg(piccolo);    

    ofstream ofile;
    ofile.open("piramide.xyz");
    ofile<<r_piccolo.size()<<endl;
    ofile<<"Pt "<<"Pt "<<"1 "<<"1"<<endl;
    for(int i=0;i<r_piccolo.size();i++){
        //cout<<i<<endl;
        ofile<<"Pt "<<r_piccolo[i][0]<<" "<<r_piccolo[i][1]<<" "<<r_piccolo[i][2]<<endl;
    }
    ofstream ofile1;
    ofile1.open("siti_jmol.xyz");
    ofile1<<r_grande.size()<<endl;
    ofile1<<"Au "<<"Au "<<"1 "<<"1"<<endl;
    for(int i=0;i<r_grande.size();i++){
        ofile1<<"Au "<<r_grande[i][0]<<" "<<r_grande[i][1]<<" "<<r_grande[i][2]<<endl;
    }

    ofstream ofile2;
    ofile2.open("siti.xyz");
    ofile2<<r_grande.size()<<endl;
    for(int i=0;i<r_grande.size();i++){
        ofile2<<r_grande[i][0]<<" "<<r_grande[i][1]<<" "<<r_grande[i][2]<<endl;
    }

//piccolo[i][0]==grande[j][0] && piccolo[i][1]==grande[j][1] && piccolo[i][2]==grande[j][2]

}
