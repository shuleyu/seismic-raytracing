#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<map>
#include<vector>
#include<algorithm>

#include<Normalize.hpp>
#include<RampFunction.hpp>

using namespace std;

int main(){

    double centerX = 30, dx = 0.005, halfWidth = 2;
    double buttonY = 2891, height = 30;
    double delta = 0.5;


    ofstream fpout("Trapzoid.txt");
    vector<double> A;
    for (size_t i = 0; i < halfWidth * 2 / dx; ++i){

        double x = fabs(dx * i - halfWidth) - (halfWidth - delta);

        if (x >= 0) A.push_back(RampFunction(delta - x, 0, delta));
        else A.push_back(1);
    }

    Normalize(A);

    for (size_t i = 0; i < A.size(); ++i) {

        fpout << centerX - halfWidth + i * dx << " " << buttonY - height * A[i] << '\n';
    }
    fpout << centerX + halfWidth << " " << buttonY << '\n';
    fpout.close();

	return 0;
}
