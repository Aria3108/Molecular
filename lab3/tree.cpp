#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <climits>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

typedef vector<vector<double>> mat;


struct Node {
    string label;
    double height;
    Node* left = nullptr;
    Node* right = nullptr;

    Node(string l, double h = 0.0) : label(l), height(h) {}
};

// Imprimir árbol en notación Newick
void printNewick(Node* node) {
    if (!node) return;
    if (!node->left && !node->right) {
        cout << node->label;
    } else {
        cout << "(";
        printNewick(node->left);
        cout << ": " << fixed << setprecision(2) << node->height - node->left->height << ", ";
        printNewick(node->right);
        cout << ": " << fixed << setprecision(2) << node->height - node->right->height << ")";
    }
}

void writeDOT(Node* root, ofstream& out) {
    if (!root) return;

    if (root->left) {
        out << "\"" << root->label << "\" -> \"" << root->left->label << "\";\n";
        writeDOT(root->left, out);
    }
    if (root->right) {
        out << "\"" << root->label << "\" -> \"" << root->right->label << "\";\n";
        writeDOT(root->right, out);
    }
}

void exportToDOT(Node* root, const string& filename) {
    ofstream out(filename);
    out << "digraph Tree {\n";
    writeDOT(root, out);
    out << "}\n";
    out.close();
}

void printMat(const mat& M) {
    for (const auto& row : M) {
        for (double val : row) {
            cout << "\t" << fixed << setprecision(4) << val << "\t";
        }
        cout << endl;
    }
}

// Diferencias entre cadenas
int contDif(string a, string b){

    int size = min(a.size(), b.size());
    int dif = 0;
    for(int i = 0; i < size; i++){
        if(a[i] != b[i]) dif++;
    }

    return dif;
}

//Arma la matriz de distancias
mat distMat( vector<string> seq){
    int n = seq.size();
    mat dist(n, vector<double>(n,0.0));

    for(int i = 0; i < n; i++){
        for(int j = i+1 ; j < n; j++){
            int dif = contDif(seq[i], seq[j]);
            double aux = (double)dif/min(seq[i].size(), seq[j].size());
            double d = -0.75 *log(1-(4.0/3.0)*aux);
            if (d == -INFINITY) {
                d = 1e9;  // O puedes elegir otro valor muy grande para evitar distancias infinitas
            }
            dist[i][j] = d;
        }
    }

    return dist;
}

// Encuentra el par más cercano
pair<int, int> best(const mat &mdist){

    pair<int, int> bestP = { -1, -1 };
    double min = INT_MAX;
    for(int i(0); i < mdist.size(); i++){
        for(int j(i+1); j < mdist.size(); j++){
            if(mdist[i][j] < min){
                min = mdist[i][j];
                bestP = {i,j};
            }
        }
    }
    return bestP;
}

Node* UPGMA(vector<string> labels, mat dist) {
    int n = labels.size();
    vector<Node*> clusters;
    vector<int> sizes(n, 1);

    for (int i = 0; i < n; i++)
        clusters.push_back(new Node(labels[i]));

    while (clusters.size() > 1) {
        pair<int, int> b = best(dist);
        if (b.first > b.second) swap(b.first, b.second);

        double newHeight = dist[b.first][b.second] / 2.0;

        Node* newNode = new Node("(" + clusters[b.first]->label + "," + clusters[b.second]->label + ")", newHeight);
        newNode->left = clusters[b.first];
        newNode->right = clusters[b.second];

        // Crear nuevo vector de distancias
        vector<double> newRow;
        for (int k = 0; k < dist.size(); k++) {
            if (k != b.first && k != b.second) {
                double newDist = (dist[b.first][k] * sizes[b.first] + dist[b.second][k] * sizes[b.second]) / (sizes[b.first] + sizes[b.second]);
                newRow.push_back(newDist);
            }
        }

        // Actualizar distancias
        mat newDistMat;
        int idx = 0;
        for (int r = 0; r < dist.size(); r++) {
            if (r == b.first || r == b.second) continue;
            vector<double> row;
            int cidx = 0;
            for (int c = 0; c < dist.size(); c++) {
                if (c == b.first || c == b.second) continue;
                row.push_back(dist[r][c]);
                cidx++;
            }
            row.push_back(newRow[idx++]); // agregar nueva distancia
            newDistMat.push_back(row);
        }

        newRow.push_back(0.0); // diagonal
        newDistMat.push_back(newRow);

        // Actualizar estructuras
        dist = newDistMat;

        clusters.erase(clusters.begin() + b.second);
        clusters.erase(clusters.begin() + b.first);
        clusters.push_back(newNode);

        sizes.erase(sizes.begin() + b.second);
        sizes.erase(sizes.begin() + b.first);
        sizes.push_back(sizes[b.first] + sizes[b.second]);
    }

    return clusters[0];
}
 
mat buildQMatrix(const mat& D) {
    int n = D.size();
    vector<double> r(n, 0.0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            r[i] += D[i][j];

    mat Q(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j)
                Q[i][j] = (n - 2) * D[i][j] - r[i] - r[j];
            else
                Q[i][j] = 1e9; // Evitar mínimos en la diagonal
        }
    }

    return Q;
}

// Encuentra el par con el menor valor en la matriz Q
pair<int, int> findMinQ(const mat& Q) {
    double minQ = 1e9;
    pair<int, int> pairIJ = {-1, -1};
    int n = Q.size();

    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            if (Q[i][j] < minQ) {
                minQ = Q[i][j];
                pairIJ = {i, j};
            }

    return pairIJ;
}

// Algoritmo Neighbor Joining
Node* NeighborJoining(vector<string> labels, mat& D) {
    int n = labels.size();
    vector<Node*> nodes(n);

    // Inicializa los nodos
    for (int i = 0; i < n; ++i)
        nodes[i] = new Node(labels[i]);

    while (D.size() > 2) {
        mat Q = buildQMatrix(D);
        pair<int, int> minPair = findMinQ(Q);
        int i = minPair.first;
        int j = minPair.second;
        if (i > j) swap(i, j);

        int m = D.size();
        vector<double> r(m, 0.0);
        for (int k = 0; k < m; ++k)
            for (int l = 0; l < m; ++l)
                r[k] += D[k][l];

        double delta = (r[i] - r[j]) / (m - 2);
        double limbLengthI = 0.5 * D[i][j] + 0.5 * delta;
        double limbLengthJ = 0.5 * D[i][j] - 0.5 * delta;

        Node* newNode = new Node("N" + to_string(n++));
        newNode->left = nodes[i];
        newNode->right = nodes[j];
        newNode->height = max(nodes[i]->height + limbLengthI, nodes[j]->height + limbLengthJ);

        // Elimina las filas y columnas i y j de la matriz de distancias
        for (int k = 0; k < D.size(); ++k) {
            if (k != i && k != j) {
                D[k].erase(D[k].begin() + j);
                if (k > i) D[k].erase(D[k].begin() + i);
            }
        }

        D.erase(D.begin() + j);  // Elimina fila j
        D.erase(D.begin() + i);  // Elimina fila i

        // Crea la nueva fila con las distancias calculadas
        vector<double> newRow(D.size(), 0.0);
        for (int k = 0; k < D.size(); ++k) {
            if (k != i && k != j) {
                newRow[k] = 0.5 * (D[i][k] + D[j][k] - D[i][j]);
            }
        }

        // Inserta la nueva fila en la matriz
        D.push_back(newRow);
        for (int k = 0; k < D.size(); ++k) {
            D[k].push_back(newRow[k]);
        }

        // Actualiza la lista de nodos
        nodes.erase(nodes.begin() + j);
        nodes.erase(nodes.begin() + i);
        nodes.push_back(newNode);
    }

    // Agregar el último nodo raíz
    Node* root = new Node("Root");
    root->left = nodes[0];
    root->right = nodes[1];
    root->height = max(D[0][1] / 2.0 + nodes[0]->height, D[0][1] / 2.0 + nodes[1]->height);

    return root;
}


int main(){

    vector<string> seq = {
        "ATTGCCATT",
        "ATGGCCATT",
        "ATCCAATTTT",
        "ATCTTCTT",
        "ACTGACC"
    };

    /* vector<string> ap = {
        "LEON", "CARAZAS", "CACERES", "TUPAC", "COLAN"
    };*/

    mat A = distMat(seq);
    printMat(A);

    pair<int, int> b = best(A);
    cout << "mejor par: " << b.first << ", " << b.second << endl<< endl;

    
    vector<string> etiquetas = { "S1", "S2", "S3", "S4", "S5" };
    
    //UPGMA
    /*Node* root = UPGMA(etiquetas, A);
    cout << "\n Arbol en notacion Newick:\n";
    printNewick(root);
    cout << ";" << endl;

    exportToDOT(root, "arbol.dot");
    cout << "\nArchivo DOT exportado: arbol.dot" << endl;
    */
    //NEIGHBOR JOIN
    cout << endl << "-----------------------NJ ----------------------" << endl;
    cout << "nj" << endl;
    Node* rootNJ = NeighborJoining(etiquetas, A);

    cout << "RootNJ label: " << rootNJ->label << endl;
    if (rootNJ->left) cout << "Left child: " << rootNJ->left->label << endl;
    if (rootNJ->right) cout << "Right child: " << rootNJ->right->label << endl;


    cout << "\n ArbolNJ en notacion Newick:\n";
    printNewick(rootNJ);
    cout << ";" << endl;

    exportToDOT(rootNJ, "arbolNJ.dot");
    cout << "\nArchivo DOT exportado: arbolNJ.dot" << endl;

    return 0;
}
