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


    double centerX = 30, dx = 0.005, halfWidth = 1;
    double height = 30, bottomY = 2891;
    double sigma = 0.8493218;



    ofstream fpout("Gaussian.txt");
    vector<double> A;
    for (double x = -halfWidth; x <= halfWidth; x += dx) {

        A.push_back(GaussianFunction(x, sigma, 0));
    }

    Normalize(A);

    for (size_t i = 0; i < A.size(); ++i){

        fpout << centerX - halfWidth + i * dx << " " << bottomY - height * A[i] << '\n';
    }
    fpout << centerX + halfWidth << " " << bottomY << '\n';
    fpout << centerX - halfWidth << " " << bottomY << '\n';
    fpout.close();

	return 0;
}
