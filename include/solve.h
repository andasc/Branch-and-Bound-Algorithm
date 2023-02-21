#pragma once
#include <iostream>
//#include <cstdlib>
#include <stack>
#include <vector>
#include <cmath>
#include <utility>

#include "simplex.h"

using namespace std;

class Solver
{   
    typedef vector<vector<double>> matrix;
    public:
    //&= Call by reference
        Solver(const matrix& m,const vector<double>& maximal, const vector<double>& rs); //Konstruktor
        pair<vector<double>,double> loese();
        

    private:
        //Objekt der Klasse Stack
        stack<matrix> probleme;
        const matrix parameter;

        //Objekt der Klasse Simplex
        Simplex simplex;
        // vector<int> xi;
        //Speicherung der bisher besten ganzzahligen Lösung
        int zi = 0;
        const vector<double> max_gleichung;
        vector<double> rechte_seite;
        const int anzahl_variablen;
        const int anzahl_gleichungen;
        pair<vector<double>,double> optimum;

        // Erleichtert dem User nur die Eingabe, damit man keine negativen Parameter setzen muss
        matrix mal_minus_eins(const matrix& m);

        // wählt nicht ganzzahliges x als Verzweigungsvariable
        pair<double,int> chose_variable(vector<double> rueckgabe);

        //fässt die Übergabe zusammen
        matrix uebergabe_zusammenfassen(const matrix& new_x, const vector<double> rs);

        // findet heraus, ob variable schon mal beschränkt wurde
        int return_index(const matrix& bedingungen, int i_von_x);

};
