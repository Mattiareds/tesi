/*
questo mi dà il reticolo di base,
faro due elementi

siti è di larghezza 7 perché le ultime due entrate riguardano il piano e se ci sono altri piani di appartenenza (max 4)

piani ha 5 entrate perché la 5 entrata è la tipologia di piano
*/
#ifndef RETICOLO_H
#define RETICOLO_H


class Reticolo{
    std::string specie_chimica;
    std::vector<std::vector<double>> piani;
    std::vector<double> sito;
    double d_pv;
    double a,b;
    int spigolo;
    int n_cut[3];
    double barriere[5];
    int metodo;

    public:
    int N;
    int mc_step;
    int current_plane;
    std::vector<std::vector<double>> siti; 
    void lettura_file(std::ifstream& ifile);
    void lettura_dati_piramide(std::ifstream& ifile);
    void lettura_dati_moto(std::ifstream& ifile);
    void scrittura_file(int i, std::ofstream& ofile);
    void aggiunte_file(std::ofstream& ofile); //per aggiungere senza intestazione
    void aggiunte_atomo(std::ofstream& ofile);
    void seleziona_sito(int i);
    void costruttore_piani();
    void associa_piano();
    bool test_piano(int i, int j);
    std::vector <int> pv(int pos);
    void dati_guscio();
    std::ofstream output_moto(int,Reticolo&);
    int metropolis(int i, int j,double r);
    
    
};
#endif