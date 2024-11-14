class Nanopiramide
{
vector<vector<double>> Piramide_fcc;
vector<vector<double>> Siti_fcc; //le prime 3 coordinate sono la posizione letta da file, l'ultima il piano di appartenenza
vector<double> Sito[3]; //il sito su cui si trova ora la particella
vector<double> Piano[4];



public:
void lettura_piramide(ifstream);
void lettura_siti(ifstream);
void scrittura_jmol(ofstream);
void associa_piano();



};