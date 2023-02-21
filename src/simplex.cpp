#include <simplex.h>

Simplex::Simplex(){}

pair<vector<double>,double> Simplex::input(matrix vorfaktoren, vector<double> max_gleichung,vector<double> rechte_seite){
    reset();
    m= vorfaktoren.size();          // Anzahl Gleichungen
    n= vorfaktoren.front().size();  // Anzahl Variablen
    for(int i=0; i< m;i++){
        for (int j=0; j< n; j++){
            A[i][j]= vorfaktoren[i][j];
        }
    }
    for (int i=0; i< m;i++){
        b[i]= rechte_seite[i];
    }
    for (int i=0; i< n;i++){
        c[i]= max_gleichung[i];
    }

    return simplex();

}

void Simplex::reset(){
    n=0;
    m=0;
    v=0;
    memset(N, 0, sizeof(N));
    memset(B, 0, sizeof(B));
    memset(c, 0, sizeof(c));
    memset(b, 0, sizeof(b));
    memset(A, 0, sizeof(A[0][0])*MAX_M*MAX_N);
}


void Simplex::pivot(int x, int y){
    
    // zuerst die x-te Zeile umordnen
    for (int j=0;j<n;j++)
    {
        if (j != y)
        {
            A[x][j] /= -A[x][y];
        }
    }
    b[x] /= -A[x][y];
    A[x][y] = 1.0 / A[x][y];
    
    // nun andere Zeilen umordnen
    for (int i=0;i<m;i++)
    {
        if (i != x)
        {
            for (int j=0;j<n;j++)
            {
                if (j != y)
                {
                    A[i][j] += A[i][y] * A[x][j];
                }
            }
            b[i] += A[i][y] * b[x];
            A[i][y] *= A[x][y];
        }
    }
    
    // nun die Zielfunktion neu anordnen
    for (int j=0;j<n;j++)
    {
        if (j != y)
        {
            c[j] += c[y] * A[x][j];
        }
    }
    v += c[y] * b[x];
    c[y] *= A[x][y];
    
    // Tausche schließlich die Basis- und Nicht-Basis-Variable
    swap(B[x], N[y]);
}

int Simplex::iterate_simplex(){

    
    int ind = -1, best_var = -1;
    for (int j=0;j<n;j++)
    {
        if (c[j] > 0)
        {
            if (best_var == -1 || N[j] < ind)
            {
                ind = N[j];
                best_var = j;
            }
        }
    }
    if (ind == -1) return 1;
    
    double max_constr = INFINITY;
    int best_constr = -1;
    for (int i=0;i<m;i++)
    {
        if (A[i][best_var] < 0)
        {
            double curr_constr = -b[i] / A[i][best_var];
            if (curr_constr < max_constr)
            {
                max_constr = curr_constr;
                best_constr = i;
            }
        }
    }
    if (isinf(max_constr)) return -1;
    else pivot(best_constr, best_var);
    
    return 0;

}

int Simplex::initialise_simplex(){
    int k = -1;
    double min_b = -1;
    for (int i=0;i<m;i++)
    {
        if (k == -1 || b[i] < min_b)
        {
            k = i;
            min_b = b[i];
        }
    }
    
    if (b[k] >= 0) // Grundlösung funktioniert!
    {
        for (int j=0;j<n;j++) N[j] = j;
        for (int i=0;i<m;i++) B[i] = n + i;
        return 0;
    }
    
    // Hilfs-LP erzeugen
    n++;
    for (int j=0;j<n;j++) N[j] = j;
    for (int i=0;i<m;i++) B[i] = n + i;
    
    // die Zielfunktion speichern
    double c_old[MAX_N];
    for (int j=0;j<n-1;j++) c_old[j] = c[j];
    double v_old = v;
    
    // Hilfszielfunktion
    c[n-1] = -1;
    for (int j=0;j<n-1;j++) c[j] = 0;
    v = 0;
    // Hilfskoeffizienten
    for (int i=0;i<m;i++) A[i][n-1] = 1;
    
    // Anfangspivot durchführen
    pivot(k, n - 1);
    
    // lösen nun aHilfszielfunktion LP
    int code;
    while (!(code = iterate_simplex()));
    
    assert(code == 1); // Hilfsfunktion LP kann nicht unbegrenzt sein!!!
    
    if (v != 0) return -1; // nicht lösbar!
    
    int z_basic = -1;
    for (int i=0;i<m;i++)
    {
        if (B[i] == n - 1)
        {
            z_basic = i;
            break;
        }
    }
    
    // wenn x_n basisch ist, einen degenerierten Pivot durchführen, um es nicht basisch zu machen
    if (z_basic != -1) pivot(z_basic, n - 1);
    
    int z_nonbasic = -1;
    for (int j=0;j<n;j++)
    {
        if (N[j] == n - 1)
        {
            z_nonbasic = j;
            break;
        }
    }
    assert(z_nonbasic != -1);
    
    for (int i=0;i<m;i++)
    {
        A[i][z_nonbasic] = A[i][n-1];
    }
    swap(N[z_nonbasic], N[n - 1]);
    
    n--;
    for (int j=0;j<n;j++) if (N[j] > n) N[j]--;
    for (int i=0;i<m;i++) if (B[i] > n) B[i]--;
    
    for (int j=0;j<n;j++) c[j] = 0;
    v = v_old;
    
    for (int j=0;j<n;j++)
    {
        bool ok = false;
        for (int jj=0;jj<n;jj++)
        {
            if (j == N[jj])
            {
                c[jj] += c_old[j];
                ok = true;
                break;
            }
        }
        if (ok) continue;
        for (int i=0;i<m;i++)
        {
            if (j == B[i])
            {
                for (int jj=0;jj<n;jj++)
                {
                    c[jj] += c_old[j] * A[i][jj];
                }
                v += c_old[j] * b[i];
                break;
            }
        }
    }
    
    return 0;
}

pair<vector<double>, double> Simplex::simplex(){
    if (initialise_simplex() == -1)
    {
        return make_pair(vector<double>(n + m, -2), INFINITY);
    }
    
    int code;
    while (!(code = iterate_simplex()));
    
    if (code == -1) return make_pair(vector<double>(n + m, -1), INFINITY);
    
    vector<double> ret;
    ret.resize(n + m);
    for (int j=0;j<n;j++)
    {
        ret[N[j]] = 0;
    }
    for (int i=0;i<m;i++)
    {
        ret[B[i]] = b[i];
    }
    
    return make_pair(ret, v);
}