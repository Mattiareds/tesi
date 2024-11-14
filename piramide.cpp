#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>

using namespace std;

/* DATI PLATINO
*passo reticolare: 3,9242 ang
distanza tra primi vicini: 2.7748 ang
*/

//DISCLAIMER: CI SONO UN MILIONE DI CASINI PERCHÉ HO USATO GLI ARRAY, USANDO I VECTOR CI SONO MODI PIÙ FURBI

//CREO I QUADRATI LUNGO ASSE X
//crea un quadrato di lato l particelle con centro nel punto desiderato 
//lo creo nel piano x=0 e poi lo traslo nell'origine e poi nel punto desiderato
double** square_maker(int l, double* c,double d_pv){
    double** square=new double* [l*l];
    for(int i=0;i<(l*l);i++){
        square[i]= new double[3];
    }
    //lo creo nell'origine e poi lo traslo
    int count=0;
    for(int i=0;i<l;i++){
        for(int j=0;j<l;j++){
            square[count][0]=0;
            square[count][1]=d_pv*j;
            square[count][2]=d_pv*i;
            count++;
        }
    }

    //lo traslo nell'origine
    double* cm=new double[3];
    cm[0]=0;
    cm[1]=0;
    cm[2]=0;

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

//creazione quadrati cutted, assumo che il cut sia uguale per tutti i vertici
double** cutted_square_maker(int l, double* c,double d_pv,int n_cut){
    double** quadrato=square_maker(l,c,d_pv); //questo è trasalto nell'origine    
    //calcolo il numero di particelle rimosse 
    int rm=0;
    for(int p=0;p<n_cut;p++) rm+=4*(p+1);

    double** cutted=new double* [(l*l)-rm];
    for(int i=0;i<((l*l)-rm);i++) cutted[i]= new double[3];    

    int added=0;
    for(int i=0;i<(l*l);i++){
        if(((fabs(quadrato[i][1])+fabs(quadrato[i][2]))<(double(l-n_cut)*(d_pv-(0.05*d_pv))))){ //traslazione linea di una frazione del passo per sicurezza, senno a causa del troncamento dava problemi
            for(int y=0;y<3;y++) {
                cutted[added][y]=quadrato[i][y];
            }
            added++;
        }
    }

return cutted;
}

int n_atoms_square(int l,int n_cut){
    int tolti=0;
    for(int i=0;i<n_cut;i++) tolti+= 4*(i+1);
    int size=l*l-tolti;
    return size;
}

int n_holes_square(int l,int n_cut){
    int uno=n_atoms_square(l,n_cut);
    int due=n_atoms_square(l-1,n_cut);
    return uno-due;
}

//creazione del reticolo di bordo: ne creo uno più grande e aggiungo solo se 
double** edge_maker(int l, double* c,double d_pv, int n_cut){
    double** grande=cutted_square_maker(l+1,c,d_pv,n_cut);
    
    int size1=n_atoms_square(l+1,n_cut);

    int size2=n_atoms_square(l,n_cut);

    double** vec=new double* [size1-size2];
    for(int i=0;i<(size1-size2);i++) vec[i]= new double[3];    

    //tutti quelli con x e y massimi, poi quelli la cui vale la relazione sopra
    int count=0;
    double z_max=0;
    double z_min=0;
    double y_max=0;
    double y_min=0;
    for(int i=0;i<(size1-size2);i++){
        if(grande[i][2]>z_max) z_max=grande[i][2];
        if(grande[i][2]<z_min) z_min=grande[i][2];
        if(grande[i][1]>y_max) y_max=grande[i][2];
        if(grande[i][1]<y_min) y_min=grande[i][2];        
    }
    
    bool test=false;
    for(int i=0;i<(size1-size2);i++){
        if(grande[i][2]==z_max) {
            vec[count]=grande[i];
            count++;
            test=true;
        }
        if(grande[i][2]==z_min) {
            vec[count]=grande[i];
            count++;
            test=true;
        }if(grande[i][1]==y_max) {
            vec[count]=grande[i];
            count++;
            test=true;
        }if(grande[i][1]==y_max) {
            vec[count]=grande[i];
            count++;
            test=true;
        }
        if(test==false){
            if(((fabs(grande[i][1])+fabs(grande[i][2]))>(double(l-n_cut)*(d_pv-(0.05*d_pv))))){
                vec[count]=grande[i];
                count++;
            }
        }
    }

    return vec;
}


//numero di atomi in  otteadro troncato
int number_of_atoms(int n_l, int* n_cut){
    int n_atoms=(1./3.)*(2*n_l*n_l*n_l+n_l);
    int n_cuts=0;
    for(int j=0;j<3;j++){
        if(j<2) for(int i=0;i<n_cut[j];i++) n_cuts+= (i+1)*(i+1);
        if(j==2) for(int i=0;i<n_cut[2];i++) n_cuts+= 4*(i+1)*(i+1);
    }
    return n_atoms-n_cuts;
}

//CREO IL QUADRATO CENTRALE, CON CM NELL'ORIGINE, POI CREO UN QUADRATO PIÙ CORTO DI 1 E LO SPOSTO SULL'ASSE X MA A X= DISTANZA TRA I DUE QUADRATI 
//angolo tra quadrato e triangolo 54,75 ovvero 0.9553

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
    
    //CREAZIONE DELLA PIRAMIDE:POSIZIONAMENTO DEI QUADRATI
    //setto il vettore centro del quadrato che verrà posizionato
    double* c=new double[3];
    for (int i=0;i<3;i++) c[i]=0;

    //creo il vettore che mi terrà conto di tutti gli atomi
    int n_atoms=number_of_atoms(n_l,n_cut);

    double** atoms=new double*[n_atoms];
    for(int i=0;i<(n_atoms);i++){
        atoms[i]= new double[3];
    }

    //vettore che mi dice le coordinate del reticolo oltre al primo
    int n_holes=number_of_atoms(n_l+1,n_cut)-number_of_atoms(n_l,n_cut);
    double** holes=new double*[n_holes];
    for(int i=0;i<(n_holes);i++){
        holes[i]= new double[3];
    }
    
    //indice che mi dice quanti atomi ho gia creato
    int atom_counter=0;
    int holes_counter=0;

    int n_1=n_atoms_square(n_l,n_cut[2]);
    int h_1=n_1-n_atoms_square(n_l+1,n_cut[2]);

    double** q=cutted_square_maker(n_l,c,d_pv,n_cut[2]);
    double** qh=edge_maker(n_l,c,d_pv,n_cut[2]);
    //creo il primo quadrato centrato nell'origine
    for(int i=0;i<n_1;i++){
        atoms[i]=q[i];
        atom_counter++;
        if(i<h_1){
            holes[i]=qh[i];
            holes_counter++;
        }
    }


    //creo gli altri nel semiasse positivo delle x contando n_cut 
    double h= d_pv*sin(((54.74)*M_PI)/180.);
    int j=1;
    while(j<(n_l-n_cut[0])){
        int l=n_l-j;
        c[0]=j*h;

        //metto quadrati tagliati finchè ha senso
        if(l>(n_l-n_cut[2])){
            int d= l*l;
            double** bmb=cutted_square_maker(l,c,d_pv,n_cut[2]-j);
            for(int m=0;m<(n_cut[2]-j);m++) d += (-1)*4*(m+1);
            for(int i=0;i<d;i++) atoms[i+atom_counter]=bmb[i];

            atom_counter=atom_counter + d;
            j++;
        }
        else{
            for(int i=0;i<(l*l);i++) atoms[i+atom_counter]=square_maker(l,c,d_pv)[i];
            atom_counter=atom_counter + (l*l);
            j++;
        }   
    }
    
    //creo asse negativo
    int y=1;
    while(y<(n_l-n_cut[1])){
        int q=n_l-y;
        c[0]=(-1)*y*h;

         //metto quadrati tagliati finchè ha senso
        if(q>(n_l-n_cut[2])){
            int d= q*q;
            for(int m=0;m<n_cut[2]-y;m++) d-= 4*(m+1);
            for(int i=0;i<d;i++) atoms[i+atom_counter]=cutted_square_maker(q,c,d_pv,n_cut[2]-y)[i];
            atom_counter=atom_counter + d;
            y++;
        }
        else{
            for(int i=0;i<(q*q);i++) atoms[i+atom_counter]=square_maker(q,c,d_pv)[i];
            atom_counter=atom_counter + (q*q);
            y++;
        }
    }
    cout<<atom_counter<<endl;

    //PRINT SU FILE XYZ
    ofstream ofile;
    ofile.open("piramide.xyz");
    ofile<<n_atoms<<endl;
    ofile<<"Pt "<<"Pt "<<"1 "<<"1"<<endl;
    
    for(int i=0;i<n_atoms;i++){
        //cout<<i<<endl;
        ofile<<"Pt "<<atoms[i][0]<<" "<<atoms[i][1]<<" "<<atoms[i][2]<<endl;
    }






    //HOLES
    //creo gli altri nel semiasse positivo delle x contando n_cut 
    int jj=1;
    while(jj<(n_l+1-n_cut[0])){
        int l=n_l+1-jj;
        c[0]=jj*h;

        //metto quadrati tagliati finchè ha senso
        if(l>(n_l+1-n_cut[2])){
            double** bmb=edge_maker(l,c,d_pv,n_cut[2]-jj);
            int d=n_holes_square(l,n_cut[2]-jj);
            for(int i=0;i<d;i++) holes[i+holes_counter]=bmb[i];

            holes_counter=holes_counter + d;
            jj++;
        }
        else{
            for(int i=0;i<(l*l);i++) holes[i+holes_counter]=square_maker(l,c,d_pv)[i];
            holes_counter=holes_counter + (l*l);
            jj++;
        }   
    }
    
    //creo asse negativo
    int yy=1;
    while(yy<(n_l+1-n_cut[1])){
        int q=n_l+1-yy;
        c[0]=(-1)*yy*h;

         //metto quadrati tagliati finchè ha senso
        if(q>(n_l+1-n_cut[2])){
            int d=n_holes_square(q,n_cut[2]-yy);
            for(int i=0;i<d;i++) holes[i+holes_counter]=edge_maker(q,c,d_pv,n_cut[2]-yy)[i];
            holes_counter=holes_counter + d;
            yy++;
        }
        else{
            for(int i=0;i<(q*q);i++) holes[i+holes_counter]=square_maker(q,c,d_pv)[i];
            holes_counter=holes_counter + (q*q);
            yy++;
        }
    }
    cout<<holes_counter<<endl;

    //PRINT SU FILE XYZ
    ofstream ofile1;
    ofile1.open("lacune.xyz");
    ofile1<<n_holes<<endl;
    ofile1<<"Au "<<"Au "<<"1 "<<"1"<<endl;
    
    for(int i=0;i<n_holes;i++){
        //cout<<i<<endl;
        ofile1<<"Au "<<holes[i][0]<<" "<<holes[i][1]<<" "<<holes[i][2]<<endl;
    }
    
}