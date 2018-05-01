#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>

#include<Interpolate.hpp>

using namespace std;

// Simple program, usage:
// ./Interpolate.out infile(oldX,oldY) newX

int main(int argc, char **argv){

    if (argc!=3) return 1;

    string infile(argv[1]),infile2(argv[2]);

    vector<double> oldX,oldY,newX;

    ifstream fpin(infile);
    double x,y;
    while (fpin >> x >> y){
        oldX.push_back(x);
        oldY.push_back(y);
    }
    fpin.close();


    fpin.open(infile2);
    while (fpin >> x ) newX.push_back(x);
    fpin.close();

    auto newY=Interpolate(oldX,oldY,newX,false);

    for (size_t i=0;i<newY.size();++i)
        printf("%.8lf %.8lf\n",newX[i],newY[i]);

    return 0;
}
