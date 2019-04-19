#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<map>
#include<vector>
#include<algorithm>

#include<Normalize.hpp>
#include<GaussianFunction.hpp>

using namespace std;

int main(){

    double X0=31.226,Y0=2891;
    double dt=0.005,L=7.5,H=7.5,sigma=0.8493218;

    vector<double> A;
    for (size_t i=0;i<L/dt;++i){
        double x=X0-L/2+dt*i;
        A.push_back(GaussianFunction(x,sigma,X0));
    }

    Normalize(A);

    ofstream fpout("Gaussian");
    for (size_t i=0;i<A.size();++i){
        fpout << X0-L/2+i*dt << " " << Y0-H*A[i] << '\n';
    }
    fpout << X0+L/2 << " " << 2891 << '\n';
    fpout << X0-L/2 << " " << 2891 << '\n';
    fpout.close();

	return 0;
}
