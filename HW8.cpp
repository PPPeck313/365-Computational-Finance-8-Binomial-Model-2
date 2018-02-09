//Preston Peck
//CS 365
//November 20, 2017
//HW8

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
using namespace std;

int binomial_simple(double S, double K, double r, double q, double sigma, double T, double t0, bool call, bool American, int n, double & V);
int binomial_test();

class Derivative {
    public:
        virtual ~Derivative() {}
        virtual double TerminalPayoff(double S) { return 0; }
        virtual int ValuationTests(double S, double & V) { return 0; }
        // data
        double r;
        double q;
        double sigma;
        double T;

    protected:
        Derivative() { r = 0; q = 0; sigma = 0; T = 0; }
};

class Option : public Derivative {
    public:
        Option() { K = 0; isCall = false; isAmerican = false; }
        virtual ~Option() {}
        virtual double TerminalPayoff(double S);
        virtual int ValuationTests(double S, double &V);
        // data
        double K;
        bool isCall;
        bool isAmerican;
};

class BinomialModel {
    public:
        BinomialModel(int n);
        ~BinomialModel();
        int FairValue(int n, Derivative* p_derivative, double S, double t0, double & V);
    private:
        // methods
        void Clear();
        int Allocate(int n);
        // data
        int n_tree;
        double **stock_nodes;
        double **derivative_nodes;
};

int main() {
    double S = 0;
    double K = 100.0;
    double r = 0.05;
    double q = 0.01;
    double sigma = 0.5;
    double T = 1.0;
    double t0 = 0.0;
    double FV_Am_put = 0;
    double FV_Eur_put = 0;
    double FV_Am_call = 0;
    double FV_Eur_call = 0;
    int n = 100;
    double dS = 0.1;
    int imax = 2000;

    ofstream outfile;
    outfile.open ("binModel2.txt");
    outfile << "S" << "FV_Am_put" << " "; outfile << "FV_Eur_put" << " "; outfile << "FV_Am_call" << " "; outfile << "FV_Eur_call" << " "; outfile << std::endl;

    for (int i = 1; i <= imax; ++i) {
        S = i * dS;
        binomial_simple(S, K, r, q, sigma, T, t0, false, true,  n, FV_Am_put);
        binomial_simple(S, K, r, q, sigma, T, t0, false, false, n, FV_Eur_put);
        binomial_simple(S, K, r, q, sigma, T, t0, true,  true,  n, FV_Am_call);
        binomial_simple(S, K, r, q, sigma, T, t0, true,  false, n, FV_Eur_call);
        // print output to file outfile << S << " ";
        outfile << S << " " << FV_Am_put << " "; outfile << FV_Eur_put << " "; outfile << FV_Am_call << " "; outfile << FV_Eur_call << " "; outfile << std::endl;
    }
    outfile.close();

    binomial_test();
}

double Option::TerminalPayoff(double S) {
    if (isCall) {
        if (S > K) {
            return S - K;
        }
    }
        
    else {
        if (S < K) {
            return K - S;
        }
    }
    return 0;
}

int Option::ValuationTests(double S, double &V) {    
    // early exercise test
    if (isAmerican) {
        if (isCall) {
            V = fmax(V, fmax(S - K, 0));
        }

        else {
            V = fmax(V, fmax(K - S, 0));
        }
    }
    return 0;
}

BinomialModel::BinomialModel(int n) {
    n_tree = 0;
    stock_nodes = 0;
    derivative_nodes = 0;
    Allocate(n);
}

BinomialModel::~BinomialModel() {
    Clear(); 
}

void BinomialModel::Clear() {
    if (stock_nodes == NULL && derivative_nodes == NULL) {
        return;
    }

    else {
        for (int i = 0; i <= n_tree; ++i) {
            delete[] stock_nodes[i];
            delete[] derivative_nodes[i];
        }

        delete[] stock_nodes;
        delete[] derivative_nodes;
    }
}

int BinomialModel::Allocate(int n) {
    if (n <= n_tree) return 0;
    // deallocate old tree
    Clear();
    // allocate memory
    n_tree = n;
    stock_nodes = new double* [n + 1];
    derivative_nodes = new double* [n + 1];

    for (int i = 0; i <= n_tree; ++i) {
        stock_nodes[i] = new double[n + 1];
        derivative_nodes[i] = new double[n + 1];

        double* S_tmp = stock_nodes[i];
        double* V_tmp = derivative_nodes[i];

        for (int j = 0; j <= n_tree; ++j) {
            S_tmp[j] = 0;
            V_tmp[j] = 0;
        }
    }
    return 0;
}

int BinomialModel::FairValue(int n, Derivative* p_derivative, double S, double t0, double & V) {     
    int rc = 0;

    V = 0;

    // validation checks
    if (n < 1 || 
        S <= 0 || 
        p_derivative == NULL || 
        p_derivative->T <= t0 ||
        p_derivative->sigma <= 0.0) {
        return 0;
    }
    
    // declaration of local variables (I use S_tmp and V_tmp)
    double* S_tmp = NULL;
    double* V_tmp = NULL;

    double r = p_derivative->r;
    double q = p_derivative->q;
    double T = p_derivative->T;
    double sigma = p_derivative->sigma;

    // calculate parameters
    double dt = (T - t0) / double(n);
    double df = exp(-r * dt);
    double growth = exp((r - q) * dt);
    double u = exp(sigma * sqrt(dt));
    double d = 1.0 / u;
    double p_prob = (growth - d)/(u - d);
    double q_prob = 1.0 - p_prob;

    // more validation checks
    if (p_prob < 0.0 || p_prob > 1.0) {
        return 1;
    }

    // allocate memory if required (call Allocate(n))
    Allocate(n);

    // set up stock prices in tree
    S_tmp = stock_nodes[0];
    S_tmp[0] = S;

    for (int i = 1; i <= n; ++i) {
        double* prev = stock_nodes[i - 1];
        S_tmp = stock_nodes[i];
        S_tmp[0] = prev[0] * d;

        for (int j = 1; j <= n; ++j) {
            S_tmp[j] = S_tmp[j - 1] * u * u;
        }
    }

    // set terminal payoff  (call virtual function in derivative class to calculate payoff)
    int i = n;//n is current scope, n_tree represent maximum tree size
    S_tmp = stock_nodes[i];
    V_tmp = derivative_nodes[i];
    
    for (int j = 0; j <= n; ++j) {
        V_tmp[j] = p_derivative->TerminalPayoff(S_tmp[j]);
    }

    // valuation loop  (call virtual function in derivative class for valuation tests)
    for (int i = n-1; i >= 0; --i) {
        S_tmp = stock_nodes[i];
        V_tmp = derivative_nodes[i];
        double* V_next = derivative_nodes[i + 1];

        for (int j = 0; j <= i; ++j) {
            V_tmp[j] = df * (p_prob * V_next[j + 1] + q_prob * V_next[j]);
            p_derivative->ValuationTests(S_tmp[j], V_tmp[j]);  // VALUATION TESTS
        }
    }
    // option fair value
    V_tmp = derivative_nodes[0];
    V = V_tmp[0];
    return 0; 
}

int binomial_test() {
    int rc = 0;
    // output file
    std::ofstream ofs("output.txt");
    double S = 100;
    double K = 100;
    double r = 0.05;
    double q = 0.01;
    double sigma = 0.5;
    double T = 1.0;
    double t0 = 0;

    Option Eur_put;
    Eur_put.r = r;
    Eur_put.q = q;
    Eur_put.sigma = sigma;
    Eur_put.T = T;
    Eur_put.K = K;
    Eur_put.isCall = false;
    Eur_put.isAmerican = false;

    Option Am_put;
    Am_put.r = r;
    Am_put.q = q;
    Am_put.sigma = sigma;
    Am_put.T = T;
    Am_put.K = K;
    Am_put.isCall = false;
    Am_put.isAmerican = true;

    Option Eur_call;
    Eur_call.r = r;
    Eur_call.q = q;
    Eur_call.sigma = sigma;
    Eur_call.T = T;
    Eur_call.K = K;
    Eur_call.isCall = true;
    Eur_call.isAmerican = false;

    Option Am_call;
    Am_call.r = r;
    Am_call.q = q;
    Am_call.sigma = sigma;
    Am_call.T = T;
    Am_call.K = K;
    Am_call.isCall = true;
    Am_call.isAmerican = true;

    double FV_Am_put = 0;
    double FV_Eur_put = 0;
    double FV_Am_call = 0;
    double FV_Eur_call = 0;
    int n = 100;
    BinomialModel binom(n);
    double dS = 0.1;
    int imax = 2000;
    int i;

    ofs << std::setw(16) << "S" << "  ";
    ofs << std::setw(16) << "FV_Am_put" << "  ";
    ofs << std::setw(16) << "FV_Eur_put" << "  ";
    ofs << std::setw(16) << "FV_Am_call" << "  ";
    ofs << std::setw(16) << "FV_Eur_call" << "  ";
    ofs << std::endl;

    for (i = 1; i <= imax; ++i) {
        S = i * dS;
        rc = binom.FairValue(n, &Am_put, S, t0, FV_Am_put);
        rc = binom.FairValue(n, &Eur_put, S, t0, FV_Eur_put);
        rc = binom.FairValue(n, &Am_call, S, t0, FV_Am_call);
        rc = binom.FairValue(n, &Eur_call, S, t0, FV_Eur_call);
        ofs << std::setw(16) << S << "  ";
        ofs << std::setw(16) << FV_Am_put << "  ";
        ofs << std::setw(16) << FV_Eur_put << "  ";
        ofs << std::setw(16) << FV_Am_call << "  ";
        ofs << std::setw(16) << FV_Eur_call << "  ";
        ofs << std::endl;
    }
    ofs.close();
    return 0; 
}

//7.1 Function signature
int binomial_simple(double S, double K, double r, double q, double sigma, double T, double t0, bool call, bool American, int n, double & V) {
    //7.2 Validation tests
    if (n < 1 || S <= 0 || T <= t0 || sigma <= 0.0) {
        return 1;
    }

    //7.3 Parameters
    double dt = (T - t0) / double(n);
    double df = exp(-r * dt);
    double growth = exp((r - q) * dt);

    double u = exp(sigma * sqrt(dt));
    double d = 1.0 / u;
    double p_prob = (growth - d) / (u - d);
    double q_prob = 1.0 - p_prob;

    if (p_prob < 0.0 || p_prob > 1.0) {
        return 1;
    }

    //7.4 Allocate memory/ set up arrays
    //array
    //allocate memory
    double** stock_nodes = new double* [n + 1];//2D
    double** option_nodes = new double* [n + 1];//2D
    double* S_tmp = NULL;
    double* V_tmp = NULL;

    for (int i = 0; i <= n; ++i) {
        stock_nodes[i] = new double[n + 1];
        option_nodes[i] = new double[n + 1];

        S_tmp = stock_nodes[i];
        V_tmp = option_nodes[i];

        for (int j = 0; j <= n; ++j) {
            S_tmp[j] = 0;
            V_tmp[j] = 0;
        }
    }

    //7.5 Set up stock prices in nodes
    S_tmp = stock_nodes[0];
    S_tmp[0] = S;

    for (int i = 1; i <= n; ++i) {
        double* prev = stock_nodes[i - 1];
        S_tmp = stock_nodes[i];
        S_tmp[0] = prev[0] * d;

        for (int j = 1; j <= n; ++j) {
            S_tmp[j] = S_tmp[j - 1] * u * u;
        }
    }
    
    //7.6 Terminal payoff
    int i = n;
    S_tmp = stock_nodes[i];
    V_tmp = option_nodes[i];

    for (int j = 0; j <= n; ++j) {
        double intrinsic = 0;
        
        if (call) {
            if (S_tmp[j] > K) {
                intrinsic = S_tmp[j] - K;
            }
        }
        
        else {
            if (S_tmp[j] < K) {
                intrinsic = K - S_tmp[j];
            }
        }

        V_tmp[j] = intrinsic;
    }
    
    //7.7 Main valuation loop
    for (int i = n - 1; i >= 0; --i) {
        S_tmp = stock_nodes[i];
        V_tmp = option_nodes[i];
        double* V_next = option_nodes[i + 1];

        for (int j = 0; j <= i; ++j) {
            V_tmp[j] = df * (p_prob * V_next[j+1] + q_prob * V_next[j]);
            
            // early exercise test
            if (American) {
                if (call) {
                    V_tmp[j] = fmax(V_tmp[j], fmax(S_tmp[j] - K, 0));
                }

                else {
                    V_tmp[j] = fmax(V_tmp[j], fmax(K - S_tmp[j], 0));
                }
            }
        } 
    }

    //7.8 option fair value
    i = 0;
    V_tmp = option_nodes[i];
    V = V_tmp[0];
        
    //7.9 Memory deallocation
    for (int i = 0; i <= n; ++i) {
        delete[] stock_nodes[i];
        delete[] option_nodes[i];
    }

    delete[] stock_nodes;
    delete[] option_nodes;
    return 0;
}