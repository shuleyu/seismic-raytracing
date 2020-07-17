#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<map>
#include<vector>
#include<algorithm>

using namespace std;

int main(){


    double halfWidth = 1, centerX = 30, dx = 0.005;
    double bottomY = 2891, height = 30;



    ofstream fpout("Ellipse.txt");
    for (double x = -halfWidth; x <= halfWidth; x += dx){

        fpout << centerX + x << " " << bottomY - sqrt((1 - x * x / halfWidth / halfWidth) * height * height) << '\n';
    }
    fpout << centerX + halfWidth << " " << bottomY << '\n';
    fpout.close();

	return 0;
}
