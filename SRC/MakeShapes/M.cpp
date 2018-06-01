#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<map>
#include<vector>
#include<algorithm>

#include<Normalize.hpp>
#include<MollifierFunction.hpp>

using namespace std;

int main(){

    ofstream fpout("Mollifier");
    double X0=34.1265,Y0=2891;
    double dt=0.005,L=2,H=30,delta=0.5;

    vector<double> A;
    for (size_t i=0;i<L/dt;++i){
        double x=fabs(dt*i-L/2)-(L/2-delta);
        if (x>=0) A.push_back(MollifierFunction(x,delta));
        else A.push_back(MollifierFunction(0,delta));
    }
    Normalize(A);

    for (size_t i=0;i<A.size();++i)
        fpout << X0-L/2+i*dt << " " << Y0-H*A[i] << '\n';
    fpout << X0+L/2 << " " << 2891 << '\n';
    fpout.close();

	return 0;
}
