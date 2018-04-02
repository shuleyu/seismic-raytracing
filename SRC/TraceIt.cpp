#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<string>
#include<iterator>
#include<algorithm>
#include<iomanip>
#include<map>

#include<CreateGrid.hpp>
#include<LineJunction.hpp>
#include<LocDist.hpp>
#include<Lon2180.hpp>
#include<Lon2360.hpp>
#include<PointInPolygon.hpp>
#include<PREM.hpp>
#include<RayPath.hpp>
#include<SegmentJunction.hpp>

using namespace std;

class Ray {
	public:
		bool IsP,GoUp,GoLeft,Valid;
        string Color;
		int InRegion,Prev;
		double Pt,Pr,TravelTime,TravelDist,RayP;

		Ray(bool p, bool g, bool l, bool v, string c, int i,double th, double r, double t, double d, double rp) :
		    IsP(p), GoUp(g), GoLeft(l), Valid(v), Color(c), InRegion(i), Prev(-1), Pt(th),Pr(r),TravelTime(t),TravelDist(d), RayP(rp) {}
};

int main(int argc, char **argv){

    enum PIenum{TS,TD,RD,CalculationStep,FLAG1};
    enum PSenum{InputRays,Layers,Depths,Polygons,OutFilePrefix,OutInfoFile,ReceiverFile,PolygonOutPrefix,FLAG2};
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

 	// Create PREM layers.
	const double RE=6371.0;
	vector<vector<double>> R{vector<double> ()},Vs=R,Vp=R;
	double prev_depth,depth,inc,next_inc,MinInc;

	ifstream fpin(PS[Layers]);
	fpin >> prev_depth >> inc;
	MinInc=inc; // Record the minimum layer increment for later use. (Solve round-off error issues)

	while (fpin >> depth >> next_inc ){

		auto tmpr=CreateGrid(RE-depth+1e-6,RE-prev_depth+1e-6,inc,2);

		if (!R[0].empty()) R[0].pop_back();
		R[0].insert(R[0].end(),make_move_iterator(tmpr.rbegin()),make_move_iterator(tmpr.rend()));

		if (next_inc<=0) break;

		prev_depth=depth;
		inc=next_inc;
		MinInc=min(MinInc,inc);
	}
	fpin.close();


 	// Create PREM velocities.
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
	// Create layers and velocities for each polygons.
	// PREM's region index is 0.
    //
    // The layers of these polygons are derived from the layers of PREM.
    // The boundaries of these derived layers are slightly larger than the polygon boudnary (in r direction).
	string tmpstr;
	double dvp,dvs;
	vector<pair<double,double>> tmpregion,tmpregion2;
	vector<vector<pair<double,double>>> Regions{tmpregion};
    vector<pair<double,double>> dV{make_pair(1,1)};

	fpin.open(PS[Polygons]);
	while (getline(fpin,tmpstr)){

        if (tmpstr.empty()) continue;

		char x;
		stringstream ss(tmpstr);
		double theta,d,r,minr,maxr;

		if (tmpstr[0]=='>') {
			if (!tmpregion.empty()) {

                // Rectify input region (tmpregion --> tmpregion2).
                // Record current region (tmpregion2) into "Regions".
                tmpregion2.clear();
                for (size_t i=0;i<tmpregion.size();++i){
                    size_t j=(i+1)%tmpregion.size();
                    double Rdist=(tmpregion[j].second-tmpregion[i].second);
                    double Tdist=(tmpregion[j].first-tmpregion[i].first);

                    size_t maxIndex=(tmpregion[i].second>tmpregion[j].second?i:j),NPTS=2;
                    double dL=RE,dR,dT;
                    while (dL>1){ // Set the target (maximum) rectified segment length.
                        NPTS*=2;
                        dR=Rdist/(NPTS-1),dT=Tdist/(NPTS-1);
                        dL=LocDist(tmpregion[maxIndex].first,0,tmpregion[maxIndex].second,tmpregion[maxIndex].first+dT,0,tmpregion[maxIndex].second+dR);
                    }
                    for (size_t k=0;k<NPTS;++k) tmpregion2.push_back(make_pair(tmpregion[i].first+k*dT,tmpregion[i].second+k*dR));
                    tmpregion2.pop_back();

                    ofstream fpout(PS[PolygonOutPrefix]+to_string(Regions.size()));
                    for (auto &item:tmpregion2) fpout << item.first << " " << item.second << '\n';
                    fpout.close();
                }
				Regions.push_back(tmpregion2);

                // Use the reigon boundaries (minr,maxr) to derive layers around the region.
                // we will use these new layers for altered velocities.
				R.push_back(vector<double> ());
				auto rit1=lower_bound(R.begin()->rbegin(),R.begin()->rend(),minr);
				auto rit2=upper_bound(R.begin()->rbegin(),R.begin()->rend(),maxr);
                auto it1=rit1.base(),it2=rit2.base();

                // new layers are slightly larger than the region.
                if (it1!=R.begin()->end()) it1=next(it1);
                if (it1!=R.begin()->end()) it1=next(it1);
                if (it2!=R.begin()->begin()) it2=prev(it2);
                if (it2!=R.begin()->begin()) it2=prev(it2);

				R.back().insert(R.back().end(),it2,it1);


                // Use the new layers for this region to get altered velocities.
                dV.push_back(make_pair(1.0+dvp/100,1.0+dvs/100));
				Vp.push_back(vector<double> ());
				Vs.push_back(vector<double> ());
				for (const auto &item:R.back()) {
					Vp.back().push_back(dV.back().first*Rvp(item));
					Vs.back().push_back(dV.back().second*Rvs(item));
				}
			}

            // Read next region. Initalizations.
			ss >> x >> dvp >> dvs;
			minr=RE+1,maxr=-1;
			tmpregion.clear();

            // The last polygon is not recorded, therefore put a ">" at the end of input.
		}
		else {
            // Read next position of current polygon.
			ss >> theta >> d;
            r=RE-d;
			minr=min(minr,r);
			maxr=max(maxr,r);
			tmpregion.push_back({theta,r});
		}
	}
	fpin.close();


    // Read in input rays.
	fpin.open(PS[InputRays]);

    char phase;
    string color;
    double theta,takeoff;
	vector<Ray> RayHeads;

    while (fpin >> theta >> depth >> takeoff >> phase >> color){

        // Starts in any polygons?
        int rid=0;
        for (size_t i=1;i<Regions.size();++i)
            if (PointInPolygon(Regions[i],make_pair(theta,RE-depth),1)) {rid=i;break;}

        // Calculate the ray parameter.
        // Ray parameter doesn't change only when reflection/refraction surface is horizontal.
        double rayp=M_PI/180*(RE-depth)*sin(fabs(takeoff)/180*M_PI)/(phase=='P'?Dvp(depth)*dV[rid].first:Dvs(depth)*dV[rid].second);

        // Record this input ray (into "RayHeads").
        RayHeads.push_back(Ray(phase=='P',fabs(takeoff)>=90,takeoff<0,true,color,rid,theta,RE-depth,0,0,rayp));
    }
    fpin.close();


	// Start ray tracing.
    double CurPosition=1;
    map<string,double> Position;
    size_t QueueL,QueueR=0; // For info output: only output valid rays.

	for (int Step=0;Step<PI[CalculationStep];++Step){

        QueueL=QueueR;
        QueueR=RayHeads.size();

		for (size_t i=QueueL;i<QueueR;++i){ // for each ray head, find next leg.

            if (!RayHeads[i].Valid) continue;

			// Find the begin and end depth for the next leg.

            /// ... among special depths.

            //// Which special depth is cloest to ray head depth?
			double RayHeadDepth=RE-RayHeads[i].Pr,NextDepth,MinDiff=RE+1;
			size_t Cloest=0;
			for (size_t j=0;j<SpecialDepths.size();++j) {
				double CurDiff=fabs(SpecialDepths[j]-RayHeadDepth);
				if (MinDiff>CurDiff) {MinDiff=CurDiff;Cloest=j;}
			}

            //// Next depth should be the cloest special depth, unless ...
            NextDepth=SpecialDepths[Cloest];
			if ( RayHeads[i].GoUp && (SpecialDepths[Cloest]>RayHeadDepth || RayHeadDepth-SpecialDepths[Cloest]<MinInc) )
                NextDepth=SpecialDepths[Cloest-1];
			else if ( !RayHeads[i].GoUp && (SpecialDepths[Cloest]<RayHeadDepth || SpecialDepths[Cloest]-RayHeadDepth<MinInc) )
				NextDepth=SpecialDepths[Cloest+1];

			double Top=min(RayHeadDepth,NextDepth),Bot=max(RayHeadDepth,NextDepth);

            /// ... among region "layering" boundaries.
			int CurRegion=RayHeads[i].InRegion;
            Top=max(Top,RE-R[CurRegion][0]);
            Bot=min(Bot,RE-R[CurRegion].back());

cout << endl << "Calculating   : ";
int x=i;
vector<int> xx;
while (x!=-1) {
    xx.push_back(x);
    x=RayHeads[x].Prev;
}
for (auto rit=xx.rbegin();rit!=xx.rend();++rit) cout << (1+*rit) << " --> ";
printf ("\nLocation      : %.15lf deg, %.15lf km (In region: %d)\n", RayHeads[i].Pt, RE-RayHeads[i].Pr,CurRegion);
cout << "Will go as    : " << (RayHeads[i].IsP?"P, ":"S, ") << (RayHeads[i].GoUp?"Up, ":"Down, ") << (RayHeads[i].GoLeft?"Left":"Right") << endl;
printf ("Start, End (D): %.16lf --> %.16lf km\n",RayHeadDepth,NextDepth);
printf ("Actual s/e (D): %.16lf --> %.16lf km\n",Top,Bot);
printf ("Search between: %.16lf ~ %.16lf km with rayp  : %.16lf sec/deg\n",RE-R[CurRegion][0],RE-R[CurRegion].back(),RayHeads[i].RayP);
printf ("\n...\n\n");
cout << flush;


			// Do ray-tracing.
			const auto &v=(RayHeads[i].IsP?Vp:Vs);
			size_t radius;
			vector<double> degree;
			auto ans=RayPath(R[CurRegion],v[CurRegion],RayHeads[i].RayP,Top,Bot,degree,radius,P[CriticalAngle]);

            // If the new leg is found to be trivia, mark it as invalid.
            size_t RayLength=degree.size();
            if (RayLength==1) {
                RayHeads[i].Valid=false;
                continue;
            }

            // If the new leg is a reflection from S(going downward) to P(going upward), and the new leg turns, also mark it as invalid.
            int PrevID=RayHeads[i].Prev;
            if (PrevID!=-1 && !RayHeads[PrevID].GoUp && !RayHeads[PrevID].IsP && RayHeads[i].GoUp && RayHeads[i].IsP && ans.second) {
                RayHeads[i].Valid=false;
                continue;
            }


			// Reverse the ray-tracing result of the next leg if next leg is going upward.
			if (RayHeads[i].GoUp) {
				double totalDist=degree.back();
				for (auto &item:degree) item=totalDist-item;
				reverse(degree.begin(),degree.end());
			}


            // Create a projection between layer index and ray index.
            int up=(RayHeads[i].GoUp?1:-1);
            auto rIndex = [RayLength,radius,up](size_t j){
                if (up==1) return radius-j;
                else return radius-RayLength+1+j;
            };


			// See if the ray path enter another region.
			int RayEnd=-1,NextRegion=-1;
			int M=(RayHeads[i].GoLeft?-1:1);

			for (size_t j=0;j<degree.size();++j){  // Let's say the ray head can't be a ray end point.

				pair<double,double> p={RayHeads[i].Pt+M*degree[j],R[CurRegion][rIndex(j)]};

                if (CurRegion){
                    // Current ray already in some polygon ...
                    if (PointInPolygon(Regions[CurRegion],p)) continue; // ... and stays in that polygon.
                    else { // ... but enters PREM.
                        RayEnd=j;
                        NextRegion=0;
                        break;
                    }
                }
                else {
                    // Current ray in PREM, search for the region it enters next.
                    for (size_t k=1;k<Regions.size();++k){
                        if (PointInPolygon(Regions[k],p)) { // If ray enters another region.
                            RayEnd=j;
                            NextRegion=k;
                            break;
                        }
                    }
                    if (RayEnd!=-1) break;
                }
			}

if (NextRegion==-1) cout << "No region change." <<  endl;
else cout << "Next Region    (T): " << NextRegion <<  endl;

			// For rays entering different regions, we need to calculate some new angles/ray parameters.
			bool ts=(PI[TS]==1),td=(PI[TD]==1),rd=(PI[RD]==1);
            double Incident=0,Rayp_td,Rayp_ts,Rayp_rd,Rayp_rs,NextPt_T,NextPr_T,NextPt_R,NextPr_R,JuncPt,JuncPr,Takeoff_ts,Takeoff_td,Takeoff_rd,Takeoff_rs;
			if (RayEnd!=-1){ // If ray ends pre-maturely.


                // For reflection, the new ray starts from the last point in the original region.
                NextPt_R=RayHeads[i].Pt+M*degree[RayEnd-1];
                NextPr_R=R[CurRegion][rIndex(RayEnd-1)];


                // For transmission/refraction, the new ray starts from the first point in the next region.
                NextPt_T=RayHeads[i].Pt+M*degree[RayEnd];
                NextPr_T=R[CurRegion][rIndex(RayEnd)];


                // Re-calculate travel distance and travel time.
				ans.first.first=ans.first.second=0;
				for (int j=0;j<RayEnd-1;++j){
					double dist=sqrt( pow(R[CurRegion][rIndex(j)],2) + pow(R[CurRegion][rIndex(j+1)],2)
                                      -2*R[CurRegion][rIndex(j)]*R[CurRegion][rIndex(j+1)]*cos(M_PI/180*(degree[j+1]-degree[j])) );
					ans.first.second+=dist;
					ans.first.first+=dist/v[CurRegion][rIndex(j+1)];
				}


                // Find the interception between the ray (RayEnd-1 ~ RayEnd) and the some polygon boundary.
                size_t L1,L2,SearchRegion=(NextRegion==0?CurRegion:NextRegion);
                pair<double,double> p2={NextPt_R,NextPr_R},q2={NextPt_T,NextPr_T};

                for (L1=0;L1<Regions[SearchRegion].size();++L1){
                    L2=(L1+1)%Regions[SearchRegion].size();
                    auto res=SegmentJunction(Regions[SearchRegion][L1],Regions[SearchRegion][L2],p2,q2);
                    if (res.first) break;
                }
                if (L1==Regions[SearchRegion].size()) cout << "!!!!!!!!!!! Can't find junction! Bugs HERE !!!!!!!!!!!" << endl;


                // Find the junction point between ray and polygon boundary.
                const pair<double,double> &p1=Regions[SearchRegion][L1],&q1=Regions[SearchRegion][L2];

                double s1=(q1.second-p1.second)/(q1.first-p1.first),s2=(q2.second-p2.second)/(q2.first-p2.first);
                auto junc=LineJunction(p1,s1,p2,s2);
                JuncPt=junc.second.first;
                JuncPr=junc.second.second;

// printf("Junction at: %.15lf,%.15lf.\n",JuncPt,JuncPr);


                // Get the geometry of the last section of the ray and the geometry boundary.
                double Rayd=180/M_PI*atan2(q2.second-p2.second,(q2.first-p2.first)*M_PI/180*JuncPr);
                double TiltAngle=180/M_PI*atan2(q1.second-p1.second,(q1.first-p1.first)*M_PI/180*JuncPr);


                // OCD on travel times and travel distance.
                double dlx=(p2.first-JuncPt)*M_PI*JuncPr/180,dly=p2.second-JuncPr;
                double dl=sqrt(dlx*dlx+dly*dly);
                ans.first.second+=dl;
                ans.first.first+=dl/v[CurRegion][rIndex(RayEnd-1)];


                // Calculate incident angle. (the "acute angle" between ray and the Normal of region boundary)
                double Rayd_Hor=Lon2360(Rayd-TiltAngle);
                Incident=fabs(Lon2180(Rayd_Hor));
                Incident=(Incident>90?Incident-90:90-Incident);

printf("Incident angle    : %.15lf deg\n",Incident);
printf("Ray angle         : %.15lf deg\n",Rayd);
printf("Tilt angle        : %.15lf deg\n",TiltAngle);


                /// A. transmission to the same wave type.
                double c1,c2;
                if (RayHeads[i].IsP) {
                    c1=Vp[CurRegion][rIndex(RayEnd)];
                    c2=c1/dV[CurRegion].first*dV[NextRegion].first;
                }
                else {
                    c1=Vs[CurRegion][rIndex(RayEnd)];
                    c2=c1/dV[CurRegion].second*dV[NextRegion].second;
                }

                //// how much difference between incident angle and takeoff angle?
                Takeoff_ts=asin(sin(Incident*M_PI/180)*c2/c1)*180/M_PI-Incident;
                ts&=!std::isnan(Takeoff_ts);

                //// take 2D structure shape into consideration, calculate the takeoff angle in earth's reference.
                if ((0<Rayd_Hor && Rayd_Hor<=90) || (180<Rayd_Hor && Rayd_Hor<=270)) Takeoff_ts=Lon2180(Rayd_Hor-Takeoff_ts+TiltAngle+90);
                else Takeoff_ts=Lon2180(Rayd_Hor+Takeoff_ts+TiltAngle+90);
                Rayp_ts=M_PI/180*NextPr_T*sin(fabs(Takeoff_ts)*M_PI/180)/c2;


                /// B. transmission to the different wave type.
                if (RayHeads[i].IsP) {
                    c1=Vp[CurRegion][rIndex(RayEnd)];
                    c2=Vs[CurRegion][rIndex(RayEnd)]/dV[CurRegion].second*dV[NextRegion].second;
                }
                else {
                    c1=Vs[CurRegion][rIndex(RayEnd)];
                    c2=Vp[CurRegion][rIndex(RayEnd)]/dV[CurRegion].first*dV[NextRegion].first;
                }

                //// how much difference between incident angle and takeoff angle?
                Takeoff_td=asin(sin(Incident*M_PI/180)*c2/c1)*180/M_PI-Incident;
                td&=!std::isnan(Takeoff_td);

                //// take 2D structure shape into consideration, calculate the takeoff angle in earth's reference.
                if ((0<Rayd_Hor && Rayd_Hor<=90) || (180<Rayd_Hor && Rayd_Hor<=270)) Takeoff_td=Lon2180(Rayd_Hor-Takeoff_td+TiltAngle+90);
                else Takeoff_td=Lon2180(Rayd_Hor+Takeoff_td+TiltAngle+90);
                Rayp_td=M_PI/180*NextPr_T*sin(fabs(Takeoff_td)*M_PI/180)/c2;


                /// C. reflection to a different wave type.
                c1=Vs[CurRegion][rIndex(RayEnd-1)];
                c2=Vp[CurRegion][rIndex(RayEnd-1)];
                if (RayHeads[i].IsP) swap(c1,c2);

                //// how much difference between incident angle and takeoff angle?
                Takeoff_rd=asin(sin(Incident*M_PI/180)*c2/c1)*180/M_PI-Incident;
                rd&=!std::isnan(Takeoff_rd);

                //// take 2D structure shape into consideration, calculate the takeoff angle in earth's reference.
                double x=Lon2360(-Rayd_Hor);
                if ((0<x && x<=90) || (180<x && x<=270)) Takeoff_rd=Lon2180(x-Takeoff_rd+TiltAngle+90);
                else Takeoff_rd=Lon2180(x+Takeoff_rd+TiltAngle+90);
                Rayp_rd=M_PI/180*NextPr_R*sin(fabs(Takeoff_rd)*M_PI/180)/c2;


                /// D. reflection to a same wave type.
                Takeoff_rs=Lon2180(-Rayd_Hor+TiltAngle+90);
                Rayp_rs=M_PI/180*NextPr_R*sin(fabs(Takeoff_rs)*M_PI/180)/(RayHeads[i].IsP?Vp[CurRegion][rIndex(RayEnd-1)]:Vs[CurRegion][rIndex(RayEnd-1)]);

			}
            else { // If ray doesn't end pre-maturelly (stays in the same region and reflect/refract on horizontal intervals)
				RayEnd=degree.size();
				NextRegion=CurRegion;
                NextPt_T=NextPt_R=RayHeads[i].Pt+M*degree[RayEnd-1];
                NextPr_T=NextPr_R=R[CurRegion][rIndex(RayEnd-1)];

                // ray parameters shouldn't change.
                Rayp_td=Rayp_ts=Rayp_rd=Rayp_rs=RayHeads[i].RayP;

                // Calculate incident angles.
                Incident=180/M_PI*asin(RayHeads[i].RayP*v[CurRegion][rIndex(RayEnd-1)]*180/M_PI/R[CurRegion][rIndex(RayEnd-1)]);

                // angles are dummies.
                Takeoff_ts=sin(Incident/180*M_PI)/v[CurRegion][rIndex(RayEnd-1)]
                Takeoff_td=Takeoff_ts;
                Takeoff_rd=(RayHeads[i].GoUp?M:M*91);
                Takeoff_rs=Takeoff_rd;
//                 Takeoff_rs=(RayHeads[i].GoUp?:90+Incident);


            } // End of dealing with rays entering another region.

printf("\nTakeOff angle (TS): %.15lf deg\n",Takeoff_ts);
printf("RayP          (TS): %.15lf\n",Rayp_ts);
printf("TakeOff angle (TD): %.15lf deg\n",Takeoff_td);
printf("RayP          (TD): %.15lf deg\n",Rayp_td);
printf("TakeOff angle (RD): %.15lf deg\n",Takeoff_rd);
printf("RayP          (RD): %.15lf deg\n",Rayp_rd);
printf("TakeOff angle (RS): %.15lf deg\n",Takeoff_rs);
printf("RayP          (RS): %.15lf deg\n\n",Rayp_rs);

cout << "RayEnd at Index  :" << RayEnd-1 << " (inclusive) / " << RayLength << endl;
printf ("RayEnd at    (R) :%.15lf deg, %.15lf (inclusive) km\n",NextPt_R,RE-NextPr_R);
printf ("RayEnd at    (T) :%.15lf deg, %.15lf (inclusive) km\n\n",NextPt_T,RE-NextPr_T);


			// Update current RayHead.
			RayHeads[i].TravelTime=ans.first.first;
			RayHeads[i].TravelDist=ans.first.second;


			// Output valid part ray paths.
			ofstream fpout;
			fpout.open(PS[OutFilePrefix]+to_string(i+1),ofstream::app);
			fpout << "> " << (RayHeads[i].IsP?"P ":"S ") << RayHeads[i].TravelTime << " sec." << endl;
			for (int j=0;j<RayEnd;++j)
				fpout << RayHeads[i].Pt+M*degree[j] << " " << R[CurRegion][rIndex(j)] << '\n';
			fpout.close();

            if (fabs(NextPr_R-RE)<MinInc) { // Ray reaches surface, output receiver info.

                fpout.open(PS[ReceiverFile],ofstream::app);

                fpout << RayHeads[i].Color << " ";

                if (Position.find(RayHeads[i].Color)==Position.end()) {
                    Position[RayHeads[i].Color]=CurPosition++;
                }
                fpout << Position[RayHeads[i].Color] << " " ;

                int x=i;
                vector<int> xx;
                while (x!=-1) {
                    xx.push_back(x);
                    x=RayHeads[x].Prev;
                }
                double tt=0;
                for (auto rit=xx.rbegin();rit!=xx.rend();++rit) {
                    tt+=RayHeads[*rit].TravelTime;
                }
                fpout << NextPt_R << " deg. " << tt << " sec. " << endl;
                fpout.close();
            }


            // Add rules of: (t)ransmission/refrection and (r)eflection to (s)ame or (d)ifferent way type.
            // Notice reflection with the same wave type is always allowed. ("rs" is always true)

            /// if ray going down and turns.
            if (!RayHeads[i].GoUp && ans.second) ts=td=rd=false;

            /// if ray ends at the surface.
            if (fabs(RE-R[CurRegion][rIndex(RayEnd-1)])<MinInc) ts=td=false;

            /// if ray goes down to CMB, no transmission S; if ray goes up to ICB, not transmission S.
            if ((!RayHeads[i].GoUp && fabs(NextPr_T-3480)<MinInc) || (RayHeads[i].GoUp && fabs(NextPr_T-1221.5)<MinInc)) {
                ts&=RayHeads[i].IsP;
                td&=(!RayHeads[i].IsP);
            }

            /// if ray goes down to ICB as P, no reflection S.
            if (!RayHeads[i].GoUp && RayHeads[i].IsP && fabs(NextPr_R-1221.5)<MinInc) rd=false;

            /// if ray goes up to CMB as P, no reflection S.
            if (RayHeads[i].GoUp && RayHeads[i].IsP && fabs(NextPr_R-3480)<MinInc) rd=false;


            // Add new rays according to the rules decided above.

            if (ts) {
                RayHeads.push_back(RayHeads[i]);
                RayHeads.back().Prev=i;
                RayHeads.back().Pt=NextPt_T;
                RayHeads.back().Pr=NextPr_T;
                RayHeads.back().RayP=Rayp_ts;
                RayHeads.back().GoUp=(fabs(Takeoff_ts)>90);
                RayHeads.back().GoLeft=(Takeoff_ts<0);
                RayHeads.back().InRegion=NextRegion;
            }

            if (td) {
                RayHeads.push_back(RayHeads[i]);
                RayHeads.back().IsP=!RayHeads.back().IsP;
                RayHeads.back().Prev=i;
                RayHeads.back().Pt=NextPt_T;
                RayHeads.back().Pr=NextPr_T;
                RayHeads.back().RayP=Rayp_td;
                RayHeads.back().GoUp=(fabs(Takeoff_td)>90);
                RayHeads.back().GoLeft=(Takeoff_td<0);
                RayHeads.back().InRegion=NextRegion;
            }

            if (rd) {
                RayHeads.push_back(RayHeads[i]);
                RayHeads.back().IsP=!RayHeads.back().IsP;
                RayHeads.back().Prev=i;
                RayHeads.back().Pt=NextPt_R;
                RayHeads.back().Pr=NextPr_R;
                RayHeads.back().RayP=Rayp_rd;
                RayHeads.back().GoUp=(fabs(Takeoff_rd)>90);
                RayHeads.back().GoLeft=(Takeoff_rd<0);
            }

            // rs is always possible.
            RayHeads.push_back(RayHeads[i]);
            RayHeads.back().Prev=i;
            RayHeads.back().Pt=NextPt_R;
            RayHeads.back().Pr=NextPr_R;
            RayHeads.back().RayP=Rayp_rs;
            RayHeads.back().GoUp=(fabs(Takeoff_rs)>90);
            RayHeads.back().GoLeft=(Takeoff_rs<0);

		} // End of rayhead loop.

    } // End of step loop.


    // Output Infos.
	ofstream fpout(PS[OutInfoFile]);
	for (size_t i=QueueL;i<QueueR;++i){
        if (RayHeads[i].Valid) {

            fpout << RayHeads[i].Color << " ";

            int x=i;
            vector<int> xx;
            while (x!=-1) {
                xx.push_back(x);
                x=RayHeads[x].Prev;
            }
            for (auto rit=xx.rbegin();rit!=xx.rend();++rit) fpout << (1+*rit) << ((*rit)==*xx.begin()?'\n':' ');
        }
    }

	fpout.close();

	return 0;
}
