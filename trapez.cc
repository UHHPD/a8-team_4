#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

// Testfunktionen als Funktor
class Pol1 {
public:
  double operator()(double x){ 
    return 3 * x + 2; 
  }
};

class Pol2 {
public:
  double operator()(double x){
    return -2 * pow(x, 2) + 3 * x + 1;
  }
};

class Gauss {
public:
  double operator()(double x){ 
    return 1 / (sqrt(M_PI * 2)) * exp(-x * x / 2); 
  }
};

// berechnet Werte nach Trapezformel von I_0 bis I_N

template<class Functor> vector<double> trapez(Functor f, double a, double b, int N) {
  vector<double> I(N + 1); // Feld mit N+1 Eintraegen
  const double h = b - a;
  I[0] = h / 2 * (f(a) + f(b));
  for (int k = 1; k <= N; ++k) {
    int n = pow(2, k);
    double sum = 0;
    for(int i=1; i<=(n-1); i++){
      sum += f(a+i*(h/n));
    }

    I[k] = (h/n) * ((f(a)/2) + (f(b)/2) + sum); // setze k-ten Wert im Feld
  }
  return I;
}

// berechnet die Richardsonextrapolation aus I(k-1)  und I(k)
double richardson(double Iprev, double I) { 
  return ((4*I)/3) - (Iprev/3); 
}

// berechet Naeherungen ueber das Romberg-Verfahren
// I: Ergebnis von trapez()
vector<vector<double>> romberg(vector<double> I) {
  const int N = I.size() - 1;
  vector<vector<double>> R(N + 1);

  for (int k = 0; k <= N; ++k) {
    R[k].push_back(I[k]);
  }
  for (int n = 1; n <= N+1; n++){
    for (int k = 0; k <= N-n; ++k){
     double extr = R[k+1][n-1] + (R[k+1][n-1] - R[k][n-1])/(pow(2, 2*n) - 1);
      R[k].push_back(extr);
    }
  }
  return R;
}

void testeAufgabe1() {
  Pol1 f;
  std::vector<double> I_f = trapez(f, 0, 3, 3);
  for (double tf : I_f) {
    std::cout << "A1: f:" << tf << " == " << 19.5 << ":" << (tf == 19.5 ? "ja" : "nein") << std::endl;
  }
  Pol2 g;
  std::cout << "A1: g(1) = 2 ?" << (g(1) == 2 ? "ja" : "nein") << std::endl;
  std::vector<double> I_g = trapez(g, 0, 3, 3);
  std::cout << "A1: g0:" << I_g[0] << " == " << -10.5 << ":" << (I_g[0] == -10.5 ? "ja" : "nein") << std::endl;
  std::cout << "A1: g1:" << I_g[1] << " == " << -3.75 << ":" << (I_g[1] == -3.75 ? "ja" : "nein") << std::endl;
  std::cout << "A1: g2:" << I_g[2] << " == " << -2.0625 << ":" << (I_g[2] == -2.0625 ? "ja" : "nein") << std::endl;
  double rich = richardson(I_g[0], I_g[1]);
  std::cout << "A1: Richardson : " << rich << " : " << (rich == -1.5 ? "ja " : "nein") << std::endl;
}

void testeAufgabe2() {
  Pol1 f;
  std::vector<std::vector<double>> Rf = romberg(trapez(f, 0, 3, 3));
  bool alle_richtig = true;
  int entries = 0;
  for (auto row : Rf) {
    for (double val : row) {
      alle_richtig &= val == 19.5;
      cout << val << endl;
      ++entries;
    }
    cout << endl;
  }
  std::cout << "A2: alle Eintraege für f sind 1.5:" << (alle_richtig ? "ja" : "nein") << std::endl;
  std::cout << "A2: korrekte Zahl an Einträgen:" << (entries == 10 ? "ja" : "nein") << std::endl;
  Pol2 g;
  std::vector<std::vector<double>> Rg = romberg(trapez(g, 0, 3, 3));
  std::cout << "A2: R[1][1] und R[2][1] für g gleich -1.5: " << ((Rg[1][1] == -1.5) && (Rg[2][1] == -1.5) ? " ja " : " nein") << std::endl;
}

int main() {
  // Testfunktionen:
  Pol1 f;
  Pol2 g;
  Gauss j;
  cout << "Funktion 1: f(0) = " << f(0) << endl;
  cout << "Funktion 2: g(0) = " << g(0) << endl;
  cout << "Gauss-Funktion: j(0) = " << j(0) << endl;

  // berechne Trapezformel fuer f
  vector<double> tf = trapez(f, 0., 3., 3);
  vector<double> tg = trapez(g, 0., 3., 3);
  vector<double> tj = trapez(j, 0., 3., 3);

  cout << "#############################################################\n";

  // Ausgabe:
  std::cout << "Trapez:\n";
  for (unsigned int i = 0; i < tf.size(); ++i) {
    cout << "Funktion f: I_" << i << " = " << tf[i] << endl;
  }
  for (unsigned int i = 0; i < tg.size(); ++i) {
    cout << "Funktion g: I_" << i << " = " << tg[i] << endl;
  }
  for (unsigned int i = 0; i < tj.size(); ++i) {
    cout << "Funktion j: I_" << i << " = " << tj[i] << endl;
  }
  cout << endl;

  cout << "Richardson-Extrapolation:\n";
  cout << "Funktion f:" << endl;
  for (unsigned int i = 1; i < tf.size(); ++i) {
    cout << "I_" << i-1 << "-" << "I_" << i <<  " = " << richardson(tf[i-1], tf[i]) << endl;
  }
  cout << "Funktion g:" << endl;
  for (unsigned int i = 1; i < tg.size(); ++i) {
    cout << "I_" << i-1 << "-" << "I_" << i << " = " << richardson(tg[i-1], tg[i]) << endl;
  }
  cout << "Funktion j:" << endl;
  for (unsigned int i = 1; i < tj.size(); ++i) {
    cout << "I_" << i-1 << "-" << "I_" << i <<  " = " << richardson(tj[i-1], tj[i]) << endl;
  }
  cout << endl;

  cout << "Romberg:\n";
  vector<vector<double> > R = romberg(tf);
  for(int k = 0, l = R.size() ;  k < l ; ++k) {
    for(int n = 0, m = R[k].size() ;  n < m ; ++n) {
      cout << R[k][n] << " ";
    }
    cout << endl;
  }
  //"Teste-Aufgaben" Aufrufe
  cout << endl << "Aufgabe 1" << endl;
  testeAufgabe1();
  cout << endl << "Aufgabe 2" << endl;
  testeAufgabe2();
  
}