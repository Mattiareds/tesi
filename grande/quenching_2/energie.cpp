#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>

using namespace std;

int main(){
    ifstream ifile("medie.out");
    string str;
    getline(ifile,str);
    double a,e;
    while(ifile>>a>>a>>a>>a>>a){
        e=a;
    }
    int i;
    cout<<"scrivi indice atomo "<<endl;
    cin>>i;

    ofstream outfile;
    outfile.open("energie.out", std::ios_base::app); // append instead of overwrite
    outfile << i<<" ";
    outfile<<setprecision(15)<<e<<endl; 

    string pino= "mkdir ";
    string pino1= "cp input.xyz ../quenching_2/";
    string pino2 = "mv p.* medie.out ../quenching_2/";
    string numero=to_string(i);
    pino=pino + numero;
    pino1=pino1 + numero;
    pino2= pino2 +numero;
    const char* command = pino.c_str();
    const char* command1 =pino1.c_str();
    const char* command2 =pino2.c_str();
    system(command);
    system(command1);
    system(command2);

    return 0;
}
