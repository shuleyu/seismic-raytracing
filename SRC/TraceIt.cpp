#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<string>
#include<iterator>
#include<algorithm>
#include<iomanip>

#include<CreateGrid.hpp>
#include<PointInPolygon.hpp>
#include<PREM.hpp>
#include<RayPath.hpp>

using namespace std;

class Ray {
	public:
		int InRegion;
		double Pt,Pr,TravelTime,TravelDist,TakeOffAngle,RayP;
		bool IsP;
		bool GoUp () {
			return fabs(this->TakeOffAngle)>=90;
		}
		bool GoLeft () {
			return this->TakeOffAngle<0;
		}
		Ray (int i,double th, double r, double t, double d, double a, double rp, bool p) :
		     InRegion(i),Pt(th),Pr(r),TravelTime(t),TravelDist(d), TakeOffAngle(a), RayP(rp), IsP(p) {}
};

int main(int argc, char **argv){

    enum PIenum{TS,TD,RD,MaximumLegs,FLAG1};
    enum PSenum{Layers,Depths,Polygons,OutFilePrefix,OutInfoFile,StartWith,FLAG2};
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

	// Calculate ray parameter.
	const double RE=6371.0;
	double vsAtSource=Dvs(P[EVDE]),vpAtSource=Dvp(P[EVDE]);
	double RaypS=M_PI/180*(RE-P[EVDE])*sin(fabs(P[TakeOffAngle])/180*M_PI)/vsAtSource,RaypP=RaypS*vsAtSource/vpAtSource;

 	// Set up Layers.
	vector<vector<double>> R{vector<double> ()},Vs=R,Vp=R;
	double prev_depth,depth,inc,next_inc,MinInc;

	ifstream fpin(PS[Layers]);
	fpin >> prev_depth >> inc;
	MinInc=inc;

	while (fpin >> depth >> next_inc ){
		if (!R[0].empty()) R[0].pop_back();
		auto tmpr=CreateGrid(RE-depth,RE-prev_depth,inc,2);
		R[0].insert(R[0].end(),make_move_iterator(tmpr.rbegin()),make_move_iterator(tmpr.rend()));
		if (next_inc<=0) break;
		inc=next_inc;
		MinInc=min(MinInc,inc);
		prev_depth=depth;
	}
	fpin.close();

	// Set up S and P velcoities for PREM.
	Vp[0].resize(R[0].size());
	Vs[0].resize(R[0].size());
	for (size_t i=0;i<R[0].size();++i){
		Vp[0][i]=Rvp(R[0][i]);
		Vs[0][i]=Rvs(R[0][i]);
	}

	// Read in special depths.
	vector<double> SpecialDepths;
	fpin.open(PS[Depths]);
	while (fpin >> depth) SpecialDepths.push_back(depth);
	sort(SpecialDepths.begin(),SpecialDepths.end());
	fpin.close();

	// Read in 2D structure polygons.
	// Set up special radius and velocities for each regions.
	// PREM is region 0.
	string tmpstr;
	double dvp,dvs;
	vector<pair<double,double>> region;
	vector<vector<pair<double,double>>> Regions{region};

	fpin.open(PS[Polygons]);
	while (getline(fpin,tmpstr)){

		char x;
		stringstream ss(tmpstr);
		double theta,radius,minr,maxr;

		if (tmpstr[0]=='>') {
			if (!region.empty()) {

				Regions.push_back(region);
				R.push_back(vector<double> ());

				auto it1=lower_bound(R.begin()->rbegin(),R.begin()->rend(),minr);
				auto it2=upper_bound(R.begin()->rbegin(),R.begin()->rend(),maxr);

				if (*prev(it2)!=maxr) R.back().push_back(maxr);
				R.back().insert(R.back().end(),it2.base(),it1.base());
				MinInc=min(MinInc,R.back()[0]-R.back()[1]);
				if (*it1!=minr) {
					MinInc=min(MinInc,R.back().back());
					R.back().push_back(minr);
				}

				Vs.push_back(vector<double> ());
				Vp.push_back(vector<double> ());
				for (auto &item:R.back()) {
					Vs.back().push_back((1.0-dvs/100)*Rvs(item));
					Vp.back().push_back((1.0-dvp/100)*Rvp(item));
				}
			}
			ss >> x >> dvp >> dvs;

			minr=RE+1,maxr=-1;
			region.clear();
		}
		else {
			ss >> theta >> radius;
			minr=min(minr,radius);
			maxr=max(maxr,radius);
			region.push_back({theta,radius});
		}
	}
	fpin.close();

cout << R[1][0] << " " << R[1].back() << endl;


	// Start ray tracing.

	vector<Ray> RayHeads;
	if (PS[StartWith]=="S") RayHeads.push_back(Ray(0,0,RE-P[EVDE],0,0,P[TakeOffAngle],RaypS,false));
	else {
		RayHeads.push_back(Ray(0,0,RE-P[EVDE],0,0,P[TakeOffAngle],RaypP,true));
		if (PS[StartWith]!="P") RayHeads.push_back(Ray(0,0,RE-P[EVDE],0,0,P[TakeOffAngle],RaypS,false));
	}

	for (int Leg=1;Leg<=PI[MaximumLegs];++Leg){

		size_t N=RayHeads.size();


		for (size_t i=0;i<N;++i){

cout << endl;
cout << "Leg #, RayID#, theta, radius, GoUp, PorS" << endl;
cout << Leg << " " << i << " " << RayHeads[i].Pt << " " << RayHeads[i].Pr << " " << (RayHeads[i].GoUp()?"Up ":"Down ") <<  (RayHeads[i].IsP?"P ":"S ") << endl;

			// Find the begin and end depth for the next leg (among special depths).
			double nextDepth,MinDiff=RE+1;
			size_t Cloest=0;
			for (size_t j=0;j<SpecialDepths.size();++j) {
				double CurDiff=fabs(SpecialDepths[j]-(RE-RayHeads[i].Pr));
				if (MinDiff>CurDiff) {
					MinDiff=CurDiff;
					Cloest=j;
				}
			}
			if (RayHeads[i].GoUp()){
				if (fabs(SpecialDepths[Cloest]-(RE-RayHeads[i].Pr))<MinInc) nextDepth=SpecialDepths[Cloest-1];
				else if (SpecialDepths[Cloest]>RE-RayHeads[i].Pr) nextDepth=SpecialDepths[Cloest-1];
				else nextDepth=SpecialDepths[Cloest];
			}
			else {
				if (fabs(SpecialDepths[Cloest]-(RE-RayHeads[i].Pr))<MinInc) nextDepth=SpecialDepths[Cloest+1];
				else if (SpecialDepths[Cloest]<RE-RayHeads[i].Pr) nextDepth=SpecialDepths[Cloest+1];
				else nextDepth=SpecialDepths[Cloest];
			}

			double Top=min(RE-RayHeads[i].Pr,nextDepth),Bot=max(RE-RayHeads[i].Pr,nextDepth);

cout << "Top/Bot" << endl;
cout << Top << " " << Bot << endl;

			// Do ray-tracing.
			int CurRegion=RayHeads[i].InRegion;

			const auto &v=(RayHeads[i].IsP?Vp:Vs);
			vector<double> degree;
			size_t radius;
			auto ans=RayPath(R[CurRegion],v[CurRegion],RayHeads[i].RayP,Top,Bot,degree,radius,P[CriticalAngle]);

			// Reverse the direction of the new leg if currently going upward.
			if (RayHeads[i].GoUp()) {
				double totalDist=degree.back();
				for (auto &item:degree) item=totalDist-item;
				reverse(degree.begin(),degree.end());
				radius=radius+1-degree.size();
			}

			// See if the ray path enter another region.
			int RayEnd=-1,NextRegion=-1;
			int M1=(RayHeads[i].GoUp()?-1:1),M2=(RayHeads[i].GoLeft()?-1:1);

			for (size_t j=0;j<degree.size();++j){

				pair<double,double> p;
				size_t rIndex=radius+M1*(j+1-degree.size());

				p=make_pair(RayHeads[i].Pt+M2*degree[j],R[CurRegion][rIndex]);

				if (CurRegion && PointInPolygon(Regions[CurRegion],p,-1)) continue;
				if (CurRegion) RayEnd=j; // Ray enters PREM.

				for (size_t k=1;k<Regions.size();++k){ // k=0 means PREM.

					if (PointInPolygon(Regions[k],p,-1)) { // Ray enters another region.
						RayEnd=j;
						NextRegion=k;
						break;
					}
				}
				if (RayEnd!=-1){ // Ray changes regions.
					if (NextRegion==-1) NextRegion=0; // Ray re-enters PREM.
					break;
				}
			}
printf ("%.16lf\n",R[CurRegion][radius+M1*(RayEnd+1-degree.size())]);

			// Re-calculate travel distance and travel time if the ray enter another region.
			if (RayEnd!=-1){
				ans.first.first=ans.first.second=0;
				for (int j=0;j<RayEnd-1;++j){
					size_t rIndex=radius+M1*(j+1-degree.size());
					double dist=sqrt(R[CurRegion][rIndex]*R[CurRegion][rIndex]+R[CurRegion][rIndex+1]*R[CurRegion][rIndex+1]-2*R[CurRegion][rIndex]*R[CurRegion][rIndex+1]*cos(M_PI/180*(degree[j+1]-degree[j])));
					ans.first.second+=dist;
					ans.first.first+=dist/v[CurRegion][rIndex];
				}
			}

			// Output valid part ray paths.
			ofstream fpout;
			fpout.open(PS[OutFilePrefix]+to_string(i+1),ofstream::app);
			fpout << "> " << (RayHeads[i].IsP?"P":"S") << endl;
			if (RayEnd==-1) {
				RayEnd=degree.size();
				NextRegion=CurRegion;
			}
			for (int j=0;j<RayEnd;++j) {
				size_t rIndex=radius+M1*(j+1-degree.size());
				fpout << RayHeads[i].Pt+M2*degree[j] << " " << R[CurRegion][rIndex] << '\n';
			}
			fpout.close();

			// Update current RayHead and add extra rays.
			bool ts=(PI[TS]==1),td=(PI[TD]==1),rd=(PI[RD]==1);
			size_t rIndex=radius+M1*(RayEnd-degree.size());
			double incidentAngle=asin(RayHeads[i].RayP*180/M_PI*v[CurRegion][rIndex]/R[CurRegion][rIndex])*180/M_PI;

			RayHeads[i].Pr=R[CurRegion][rIndex];
			RayHeads[i].Pt+=M2*degree[RayEnd-1];
			RayHeads[i].TravelTime+=ans.first.first;
			RayHeads[i].TravelDist+=ans.first.second;

// 			if (CurRegion==NextRegion) {

				if (!RayHeads[i].GoUp() && ans.second) ts=td=rd=false;
				if (fabs(RE-R[CurRegion][rIndex])<MinInc) ts=td=false;
				if ((!RayHeads[i].GoUp() && fabs(RayHeads[i].Pr-3480)<MinInc) || (RayHeads[i].GoUp() && fabs(RayHeads[i].Pr-1221.5)<MinInc)) {
					ts&=RayHeads[i].IsP;
					td&=(!RayHeads[i].IsP);
				}
				if (!RayHeads[i].GoUp() && RayHeads[i].IsP && fabs(RayHeads[i].Pr-1221.5)<MinInc) rd=false;
				if (RayHeads[i].GoUp() && RayHeads[i].IsP && fabs(RayHeads[i].Pr-3480)<MinInc) rd=false;

				if (ts) {
					RayHeads.push_back(RayHeads[i]);
					if (RayHeads[i].GoUp()) RayHeads.back().TakeOffAngle=M2*91;
					else RayHeads.back().TakeOffAngle=M2;

					if (CurRegion!=NextRegion) {
						RayHeads.back().InRegion=NextRegion;
					}
				}

				if (td) {
					RayHeads.push_back(RayHeads[i]);
					if (RayHeads[i].GoUp()) RayHeads.back().TakeOffAngle=M2*91;
					else RayHeads.back().TakeOffAngle=M2;
					RayHeads.back().IsP=!RayHeads.back().IsP;

					if (CurRegion!=NextRegion) {
						RayHeads.back().InRegion=NextRegion;
					}
				}

				if (rd) {
					RayHeads.push_back(RayHeads[i]);
					if (RayHeads[i].GoUp()) RayHeads.back().TakeOffAngle=M2;
					else RayHeads.back().TakeOffAngle=M2*91;
					RayHeads.back().IsP=!RayHeads.back().IsP;
				}

				// Do RS part.

				if (RayHeads[i].GoUp()) RayHeads[i].TakeOffAngle=M2;
				else {
					if (ans.second) RayHeads[i].TakeOffAngle=M2*(180-incidentAngle);
					else RayHeads[i].TakeOffAngle=M2*91;
				}
// 			}
// 			else {

// 			}
		}
	}

	ofstream fpout(PS[OutInfoFile]);
	for (size_t i=0;i<RayHeads.size();++i)
		fpout << i+1 << " " << RayHeads[i].TravelTime << " sec. " << RayHeads[i].Pt << " deg. " << RayHeads[i].RayP << " sec/deg. " << RayHeads[i].TravelDist << "km." << endl;
	fpout.close();

	return 0;
}
