#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<map>
#include<vector>
#include<algorithm>

using namespace std;

int main(){

    ofstream fpout("Ellipse");
    double A=1,B=30,X0=34.1265,Y0=2891;
    for (double x=-A;x<=A;x+=0.005){
        fpout << x+X0 << " " << Y0-sqrt((1-x*x/A/A)*B*B) << '\n';
    }
    fpout << X0+A << " " << 2891 << '\n';
    fpout.close();

	return 0;
}
