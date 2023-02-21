#include <solve.h>
typedef std::vector<std::vector<double>> matrix;

int main()
{
    // Paramter der Gleichungen "auf der linken Seite"
    matrix vorfaktoren {{1,1},
                        {5,9}};

    // Parameter der Gelichung, die den maximalen Wert annehmen soll
    std::vector<double> vorfaktoren_max_gleichung {5,8};

    // Paramter der Gleichungen "auf der rechten Seite"
    std::vector<double> rechte_seiten {6,45};
    
    Solver bab(vorfaktoren,vorfaktoren_max_gleichung,rechte_seiten);

    std::pair<std::vector<double>,double> ergebnis = bab.loese();

    int k= vorfaktoren.front().size()-1;

    std::cout << "x-Werte: ["; 
    for(int i = 0; i < k; i++){
        if (ergebnis.first[i]==0) ergebnis.first[i]= abs(ergebnis.first[i]);
         std:: cout << ergebnis.first[i] << "x"<<i+1 <<", "; 
    }
    if (ergebnis.first[k]==0) ergebnis.first[k]= abs(ergebnis.first[k]);    
    std::cout << ergebnis.first[k] << "x" << k+1 << "]";
    std:: cout << "\nOptimalster Wert: " << ergebnis.second << "\n";
}