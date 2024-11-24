#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <unistd.h>

using namespace std;

int main(){
    //mi leggo l'energia
    ifstream ifile("medie.out");
    string str;
    getline(ifile,str);
    double a,e;
    while(ifile>>a>>a>>a>>a>>a){
        e=a;
    }
    ifile.close();

    //da termianale mi leggo l'indice del sito e prendo dal file con le coordinate la x
    int i;
    cout<<"scrivi indice del sito "<<endl;
    cin>>i;

    ifstream ifile2("intermedie.txt");
    for(int j=0;j<i;j++) getline(ifile2,str);
    double x;
    ifile2>>x;

    //prendo da terminale la cartella dove reindirizzarmi
    string c_name;
    cout<<"Digitare il nome della cartella destinazione: "<<endl;
    cin>>c_name;

    //mi sposto nella cartella
    char currentDir[1024];

    // Ottieni la directory corrente
    if (getcwd(currentDir, sizeof(currentDir)) == nullptr) {
        perror("Errore nel recuperare la directory corrente");
        return 1;
    }

    std::cout << "Directory corrente: " << currentDir << std::endl;

    // Cambia directory
    string ss="../barriere/";
    ss=ss+c_name;
    if (chdir(ss.c_str()) == 0) {
        std::cout << "Cambiata directory a: /nuova/cartella\n";
    } else {
        perror("Errore nel cambiare directory");
        return 1;
    }

    // Conferma cambiamento
    system("pwd");


    //aggiorno il file con coordinata x ed energia 
    ofstream outfile;
    outfile.open("energie.out", std::ios_base::app); // append instead of overwrite
    outfile <<setprecision(10)<<x<<" ";
    outfile<<setprecision(15)<<e<<endl; 

    //creo una nuova cartella con la posizione dell'atomo
    string s2= "mkdir pos";
    s2=s2+to_string(i);
    system(s2.c_str());

    // Torna alla directory iniziale
    if (chdir(currentDir) == 0) {
        std::cout << "Tornato alla directory iniziale: " << currentDir << std::endl;
    } else {
        perror("Errore nel tornare alla directory iniziale");
        return 1;
    }
    
    //ci muovo medie, i p.* e il file di input 
    string s3 = "mv medie.out p.* ../barriere/";
    s3=s3+c_name+"/pos"+to_string(i)+"/";
    system(s3.c_str());

    string s4 = "cp input.xyz ../barriere/";
    s4 = s4+c_name+"/pos"+to_string(i)+"/";
    system(s4.c_str());

}