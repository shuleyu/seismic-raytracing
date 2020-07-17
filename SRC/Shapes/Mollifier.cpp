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

    double centerX = 30, dx = 0.005, halfWidth = 1;
    double height = 30, bottomY = 2891;
    double delta=0.5;




    ofstream fpout("Mollifier.txt");
    vector<double> A;

    for (size_t i = 0; i < halfWidth * 2 / dx; ++i){
        
        double x = fabs(dx * i - halfWidth) - (halfWidth - delta);

        if (x >= 0) A.push_back(MollifierFunction(x, delta));
        else A.push_back(MollifierFunction(0, delta));
    }

    Normalize(A);

    for (size_t i = 0; i < A.size(); ++i) {

        fpout << centerX - halfWidth + i * dx << " " << bottomY - height * A[i] << '\n';
    }
    fpout << centerX + halfWidth << " " << bottomY << '\n';
    fpout.close();

	return 0;
}
