#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<string>
#include<iterator>

#include<CreateGrid.hpp>
#include<PREM.hpp>
#include<RayPath.hpp>

using namespace std;

int main(int argc, char **argv){

    enum PIenum{FLAG1};
    enum PSenum{Layers,Ray,OutFile,FLAG2};
    enum Penum{TakeOffAngle,EVDE,CriticalAngle,FLAG3};

    /****************************************************************

				Deal with inputs. (Store them in PI,PS,P)

    ****************************************************************/

	if (argc!=4){
		cerr << "In C++: Argument Error!" << endl;
		return 1;
	}

    int int_num,string_num,double_num;

    vector<int> PI;
    vector<string> PS;
    vector<double> P;

    int_num=atoi(argv[1]);
    string_num=atoi(argv[2]);
    double_num=atoi(argv[3]);

	if (FLAG1!=int_num){
		cerr << "In C++: Ints name count error!" << endl;
		return 1;
	}
	if (FLAG2!=string_num){
		cerr << "In C++: Strings name count error!" << endl;
		return 1;
	}
	if (FLAG3!=double_num){
		cerr << "In C++: Doubles name count error!" << endl;
		return 1;
	}

	string InputTmpStr;
	int InputTmpInt,InputTmpCnt;
	double InputTmpDouble;

	InputTmpCnt=0;
	while (getline(cin,InputTmpStr)){
		++InputTmpCnt;
		stringstream ss(InputTmpStr);
		if (InputTmpCnt<=int_num){
			if (ss >> InputTmpInt && ss.eof()){
				PI.push_back(InputTmpInt);
			}
			else{
				cerr << "In C++: Ints reading Error !" << endl;
				return 1;
			}
		}
		else if (InputTmpCnt<=int_num+string_num){
			PS.push_back(InputTmpStr);
		}
		else if (InputTmpCnt<=int_num+string_num+double_num){
			if (ss >> InputTmpDouble && ss.eof()){
				P.push_back(InputTmpDouble);
			}
			else{
				cerr << "In C++: Doubles reading Error !" << endl;
				return 1;
			}
		}
		else{
			cerr << "In C++: Redundant inputs !" << endl;
			return 1;
		}
	}
	if (InputTmpCnt!=int_num+string_num+double_num){
		cerr << "In C++: Not enough inputs !" << endl;
		return 1;
	}

    /****************************************************************

                              Job begin.

    ****************************************************************/

	const double RE=6371.0;

 	// Set up Layers.
	ifstream fpin;
	vector<double> r,vs,vp,tmpr;
	double depth,inc,next_inc,prev_depth,min_inc;

	fpin.open(PS[Layers]);
	fpin >> prev_depth >> inc;
	min_inc=inc;
	while (fpin >> depth >> next_inc ){
		if (!r.empty()) r.pop_back();
		tmpr=CreateGrid(RE-depth,RE-prev_depth,inc,2);
		r.insert(r.end(),make_move_iterator(tmpr.rbegin()),make_move_iterator(tmpr.rend()));
		if (next_inc<=0) break;
		inc=next_inc;
		min_inc=min(min_inc,inc);
		prev_depth=depth;
	}
	fpin.close();

	// Set up S and P velcoities.
	for (const auto &item:r) {
		vp.push_back(Rvp(item));
		vs.push_back(Rvs(item));
	}

	// Set up Raypath.
	char by,type;
	vector<char> ChangeBy,ChangeToWaveType;
	vector<double> ChangeDepth;
	fpin.open(PS[Ray]);
	while (fpin >> depth >> by >> type) {
		if (ChangeDepth.empty()) depth=P[EVDE];
		ChangeDepth.push_back(depth);
		ChangeBy.push_back(by);
		ChangeToWaveType.push_back(type);
	}
	fpin.close();

	// Calculate ray parameter.
	double vsAtSource=Dvs(P[EVDE]),vpAtSource=Dvp(P[EVDE]);
	double rayp=M_PI/180*(RE-P[EVDE])*sin(P[TakeOffAngle]/180*M_PI)/(ChangeToWaveType[0]=='P'?vpAtSource:vsAtSource);

	// Do each ray path section.

	ofstream fpout(PS[OutFile]);
	double Distance=0,TravelDist=0,TravelTime=0;
	bool GoUp=false;
	for (size_t i=0;i+1<ChangeDepth.size();++i){

		if (ChangeBy[i]=='R') GoUp=!GoUp;

		vector<double> degree,radius;
		const auto &v=(ChangeToWaveType[i]=='P'?vp:vs);
		double Top=min(ChangeDepth[i],ChangeDepth[i+1]),Bot=max(ChangeDepth[i],ChangeDepth[i+1]);
		auto ans=RayPath(r,v,rayp,Top,Bot,degree,radius,P[CriticalAngle]);

		if (!GoUp)
			for (size_t j=0;j<degree.size();++j)
				fpout << Distance+degree[j] << " " << radius[j] << '\n';
		else
			for (size_t j=degree.size();j>0;--j)
				fpout << Distance+degree.back()-degree[j-1] << " " << radius[j-1] << '\n';

		Distance+=degree.back();
		TravelTime+=ans.first;
		TravelDist+=ans.second;

		// if ray turns.
		if (!GoUp && radius.back()-(RE-Bot)>min_inc) GoUp=true;

	}
	fpout.close();

	cout << "TakeOff    : " << P[TakeOffAngle] << " deg."     << endl;
	cout << "Depth      : " << P[EVDE]         << " km."      << endl;
	cout << "Rayp       : " << rayp            << " sec/deg." << endl << endl;
	cout << "Distance   : " << Distance        << " deg."     << endl;
	cout << "TravelTime : " << TravelTime      << " sec."     << endl;
	cout << "TravelDist: " << TravelDist     << " km."      << endl;

    return 0;
}
