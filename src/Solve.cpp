#include "solve.h"

using namespace std;

//Konstruktor, rs= rechte Seite
Solver::Solver(const matrix& m,const vector<double>& max, const vector<double>& rs) : 
    parameter(mal_minus_eins(m)), max_gleichung(max), rechte_seite(rs), anzahl_variablen(m.front().size()), anzahl_gleichungen(m.size())
{
    ;
}

Solver::matrix Solver::mal_minus_eins(const matrix& m){
    matrix ret_m(m);
    for (vector<double>& zeile : ret_m){
        for (double& d : zeile){
            d *=-1;
        }
    }
    return ret_m;
}

pair<double,int> Solver::chose_variable(vector<double> rueckgabe){
    double vv;
    for (int i=0; i<parameter.front().size(); i++){
        if(rueckgabe[i]==(int)rueckgabe[i]) continue;
        vv= rueckgabe[i];
        return make_pair(vv, i);
        
    } 
    return make_pair(-1,-1);
}

Solver::matrix Solver::uebergabe_zusammenfassen(const matrix& new_x, const vector<double> rs)
{
    matrix vorfaktoren(parameter);
    for (const vector<double>& element: new_x){
        vorfaktoren.push_back(element);
    }
    vorfaktoren.push_back(rs);

    return vorfaktoren;

}

int Solver::return_index(const matrix& bedingungen, int i_von_x){
    for(int i = 0; i< bedingungen.size()-1; i++){
        if(bedingungen[i][i_von_x]!=0)return i;
    }
    return -1;
}


pair<vector<double>,double> Solver::loese(){

    pair<vector<double>,double> rueckgabe = simplex.input(parameter,max_gleichung, rechte_seite);
    double aktuelles_z= rueckgabe.second;

    // Ausgabe Error, falls notwendig
    if(isinf(aktuelles_z)){
        if (rueckgabe.first[0] == -1) throw logic_error("Gleichungssystem nicht gedeckelt!\n");
        else if (rueckgabe.first[0] == -2) throw logic_error("Gleichungssystem nicht lösbar\n");
    } 

    //prüfen, ob z sowie alle zugehörigen x-Werte ganzzahlig sind und größer zi
    //Wenn ja, aktualisere optimalsten Wert z
    if (aktuelles_z==(int)aktuelles_z && zi<aktuelles_z && chose_variable(rueckgabe.first).second == -1){
        zi = aktuelles_z; 
        optimum = rueckgabe;
    }
    //Verzweigungsvariable wählen
    pair<double,int> a = chose_variable(rueckgabe.first);
    double verzweigungsvariable = a.first;
    int i_von_x = a.second;

    if(i_von_x==-1) return rueckgabe;

    //Verzweigungsvariable nach oben und unten runden
    int obere_schranke_x = ceil(verzweigungsvariable);  
    int untere_schranke_x = verzweigungsvariable;

     // Koeffizient der Verzweigungsvariable (der neuen Bedingungen) im Vector speichern 
    matrix neue_bedingung_1(1,vector<double> (parameter.front().size(),0));
    neue_bedingung_1[0][i_von_x]= -1;
    matrix neue_bedingung_2(1,vector<double> (parameter.front().size(),0));
    neue_bedingung_2[0][i_von_x]= 1;
    

    //rechte seiten der neuen Bedingungen im Vector speichern
    //in den einen Knoten muss...
    vector<double> rs_neue_bedingungen2 {-(double)obere_schranke_x};
    //in den anderen...
    vector<double> rs_neue_bedingungen1 {(double)untere_schranke_x};

    neue_bedingung_1.push_back(rs_neue_bedingungen1);
    neue_bedingung_2.push_back(rs_neue_bedingungen2);
    
    // lege neue Probleme auf Stapel
    probleme.push(neue_bedingung_1);
    probleme.push(neue_bedingung_2);
   
    //Solange Stapel an Problemen nicht leer...
    while (!probleme.empty())
    {

        // hole oberstes Problem vom Stapel
        matrix aktuelle_bedingungen = probleme.top();
        probleme.pop();
        
        matrix aktuell(parameter);
        aktuell.insert(aktuell.end(), aktuelle_bedingungen.begin(), aktuelle_bedingungen.end()-1);  //.begin()/end() gibt Speicheradresse
        
        vector<double> wert(rechte_seite);

                    // stelle   
        wert.insert(wert.end(),aktuelle_bedingungen.back().begin(), aktuelle_bedingungen.back().end());

    
        rueckgabe = simplex.input(aktuell,max_gleichung,wert);

        double aktuelles_z= rueckgabe.second;

        if(isinf(aktuelles_z)) continue;
    
        a = chose_variable(rueckgabe.first);
        verzweigungsvariable = a.first;
        i_von_x = a.second;

        //prüfen, ob z ganzzahlig und größer zi
        if (aktuelles_z==(int)aktuelles_z && zi<aktuelles_z && i_von_x==-1){
            zi = aktuelles_z; 
            optimum = rueckgabe;
        }

        

        // Abbruchbedingung: wenn alle x-werte ganzzahlig sind, werden keine Nachfolger-Probleme erzeugt
        if(i_von_x==-1) {
            continue;
        }

        //Verzweigungsvariable nach oben und unten runden
        obere_schranke_x = ceil(verzweigungsvariable);  
        untere_schranke_x = verzweigungsvariable;

        matrix aktuelle_bedingungen2(aktuelle_bedingungen);
        
        // Index der X-Stelle ermitteln, um genau dieses X einzuschränken und somti nur eine Stelle in der Matrix zu ändern
        int new_index_x = return_index(aktuelle_bedingungen, i_von_x);

        // Schranken zuweisen
        if (new_index_x!=-1){
            aktuelle_bedingungen[new_index_x][i_von_x] = -1;
            aktuelle_bedingungen2[new_index_x][i_von_x] = 1;

            aktuelle_bedingungen.back()[new_index_x] = untere_schranke_x;
            aktuelle_bedingungen2.back()[new_index_x] = -obere_schranke_x;
        }
        else{
            aktuelle_bedingungen.emplace(aktuelle_bedingungen.end()-1, vector<double> (aktuelle_bedingungen.front().size(),0));
            aktuelle_bedingungen[aktuelle_bedingungen.size()-2][i_von_x] = -1;
            aktuelle_bedingungen.back().push_back(untere_schranke_x);

            aktuelle_bedingungen2.emplace(aktuelle_bedingungen2.end()-1, vector<double> (aktuelle_bedingungen2.front().size(), 0));
            aktuelle_bedingungen2[aktuelle_bedingungen2.size()-2][i_von_x] = 1;
            aktuelle_bedingungen2.back().push_back(-obere_schranke_x);
        }
        
        probleme.push(aktuelle_bedingungen);
        probleme.push(aktuelle_bedingungen2);

        
        
    } 
    
    return optimum;
    
}




