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

    enum PIenum{TS,TD,RD,CalculationStep,FLAG1};
    enum PSenum{InputRays,Layers,Depths,Polygons,OutFilePrefix,OutInfoFile,FLAG2};
    enum Penum{CriticalAngle,FLAG3};

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

 	// Set up Layers.
	const double RE=6371.0;
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
	vector<pair<double,double>> region,RegionRMinMax{{0,6371}};
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

				auto rit1=lower_bound(R.begin()->rbegin(),R.begin()->rend(),minr);
				auto rit2=upper_bound(R.begin()->rbegin(),R.begin()->rend(),maxr);
                auto it1=rit1.base(),it2=rit2.base();
                if (it1!=R.begin()->end()) it1=next(it1);
                if (it2!=R.begin()->begin()) it2=prev(it2);

				R.back().insert(R.back().end(),it2,it1);
                RegionRMinMax.push_back({minr,maxr});

				Vs.push_back(vector<double> ());
				Vp.push_back(vector<double> ());
				for (auto &item:R.back()) {
					Vs.back().push_back((1.0+dvs/100)*Rvs(item));
					Vp.back().push_back((1.0+dvp/100)*Rvp(item));
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

    // Read in first legs.
	fpin.open(PS[InputRays]);
	vector<Ray> RayHeads;

    double theta,takeoff;
    char phase;
    while (fpin >> theta >> depth >> takeoff >> phase){

        // ray parameter.
        double rayp=M_PI/180*(RE-depth)*sin(fabs(takeoff)/180*M_PI)/(phase=='P'?Dvp(depth):Dvs(depth));
        
        // starts in which region.
        int rid=0;
        for (size_t i=1;i<Regions.size();++i) {
            if (PointInPolygon(Regions[i],make_pair(theta,RE-depth),1)) {
                rid=i;
                break;
            }
        }
        RayHeads.push_back(Ray(rid,theta,RE-depth,0,0,takeoff,rayp,phase=='P'));
    }
    fpin.close();

	// Start ray tracing.

    size_t PrevN=0;
	for (int Leg=1;Leg<=PI[CalculationStep];++Leg){

		size_t N=RayHeads.size();
        PrevN=N;

        // Marks of trivia rays for delete.
        vector<size_t> D;

		for (size_t i=0;i<N;++i){

cout << endl;
cout << "Leg #, RayID#, Region#, theta, radius, GoUp, PorS" << endl;
cout << Leg << " " << i << " " << RayHeads[i].InRegion << " " << RayHeads[i].Pt << " " << RayHeads[i].Pr << " " << (RayHeads[i].GoUp()?"Up ":"Down ") <<  (RayHeads[i].IsP?"P ":"S ") << endl;

			// Find the begin and end depth for the next leg (among special depths).
			int CurRegion=RayHeads[i].InRegion;
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
            Top=max(Top,RE-RegionRMinMax[CurRegion].second);
            Bot=min(Bot,RE-RegionRMinMax[CurRegion].first);

printf ("Top, Bot: %.16lf %.16lf\n",Top,Bot);
printf ("Search between: %.16lf %.16lf\n",RE-R[CurRegion][0],RE-R[CurRegion].back());

			// Do ray-tracing.

			const auto &v=(RayHeads[i].IsP?Vp:Vs);
			vector<double> degree;
			size_t radius;
			auto ans=RayPath(R[CurRegion],v[CurRegion],RayHeads[i].RayP,Top,Bot,degree,radius,P[CriticalAngle]);

if (Leg==2 && i==2) printf("TravelTime: %.15lf\n",ans.first.first);
if (Leg==3 && i==2) printf("TravelTime: %.15lf\n",ans.first.first);
if (Leg==4 && i==2) printf("TravelTime: %.15lf\n",ans.first.first);

cout << R[CurRegion].size() << " " << degree.size() << " " << radius << endl;

            // If the new leg is trivia, mark it for deleting.
            if (degree.size()==1) {
                D.push_back(i);
                continue;
            }

			// Reverse the direction of the next leg if next leg is going upward.
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

				if (CurRegion && PointInPolygon(Regions[CurRegion],p,1)) continue;
				if (CurRegion) RayEnd=j; // Ray enters PREM.

                // Current ray in PREM, search for the region it enters next.
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

cout << "RayEnd at step: " << RayEnd << " / " << degree.size()  << endl;
// printf ("Last point radius: %.16lf\n",R[CurRegion][radius+M1*(RayEnd-degree.size())]);

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
		}

        // Delete trivia rays.
        for (auto rit=D.rbegin();rit!=D.rend();++rit) {
            RayHeads.erase(RayHeads.begin()+*rit);
            PrevN-=D.size();
        }
	}

	ofstream fpout(PS[OutInfoFile]);
	for (size_t i=0;i<PrevN;++i)
		fpout << i+1 << " " << RayHeads[i].TravelTime << " sec. " << RayHeads[i].Pt << " deg. " << RayHeads[i].RayP << " sec/deg. " << RayHeads[i].TravelDist << "km." << endl;
	fpout.close();

	return 0;
}
