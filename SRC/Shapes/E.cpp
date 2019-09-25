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
    double A=2.5,B=30,X0=31,Y0=2891;
    for (double x=-A;x<=A;x+=0.005){
        fpout << x+X0 << " " << Y0-sqrt((1-x*x/A/A)*B*B) << '\n';
        printf("%.7lf %.7lf\n",x+X0,Y0-sqrt((1-x*x/A/A)*B*B));
    }
    fpout << X0+A << " " << 2891 << '\n';
    printf("%.7lf %.7lf\n",X0+A,2891.0);
    fpout.close();

	return 0;
}
