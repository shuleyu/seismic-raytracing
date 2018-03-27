#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<string>
#include<iterator>
#include<algorithm>
#include<iomanip>

#include<CreateGrid.hpp>
#include<LineJunction.hpp>
#include<LocDist.hpp>
#include<PointInPolygon.hpp>
#include<PREM.hpp>
#include<RayPath.hpp>
#include<SegmentJunction.hpp>

using namespace std;

class Ray {
	public:
		bool IsP,GoUp,GoLeft;
        string Color;
		int InRegion;
		double Pt,Pr,TravelTime,TravelDist,RayP;

		Ray(bool p, bool g, bool l, string c, int i,double th, double r, double t, double d, double rp) :
		    IsP(p), GoUp(g), GoLeft(l), Color(c), InRegion(i),Pt(th),Pr(r),TravelTime(t),TravelDist(d), RayP(rp) {}
};

int main(int argc, char **argv){

    enum PIenum{TS,TD,RD,CalculationStep,FLAG1};
    enum PSenum{InputRays,Layers,Depths,Polygons,OutFilePrefix,OutInfoFile,PolygonOutPrefix,FLAG2};
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

		auto tmpr=CreateGrid(RE-depth,RE-prev_depth,inc,2);

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
        RayHeads.push_back(Ray(phase=='P',fabs(takeoff)>=90,takeoff<0,color,rid,theta,RE-depth,0,0,rayp));
    }
    fpin.close();


	// Start ray tracing.
    size_t ValidN=0; // For info output: only output valid rays.

	for (int Leg=1;Leg<=PI[CalculationStep];++Leg){

		size_t N=RayHeads.size();
        ValidN=N;

        // Indices of trivia rays(no propagation due to past-critical angle, etc.).
        // delete them when this leg is complete.
        vector<size_t> D;

		for (size_t i=0;i<N;++i){ // for each ray head, find next leg.

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

cout << endl;
cout << "Leg and RayID : " << Leg << " " << i << endl;
cout << "In region     : " << RayHeads[i].InRegion << endl;
printf ("Location      : %.15lf deg, %.15lf km\n", RayHeads[i].Pt, RE-RayHeads[i].Pr);
cout << "Go up or down : " << (RayHeads[i].GoUp?"Up ":"Down ") << endl;
cout << "will be P or S: " << (RayHeads[i].IsP?"P ":"S ") << endl;
printf ("RayHD,TargetD : %.16lf --> %.16lf km\n",RayHeadDepth,NextDepth);
printf ("SrhTop,SrhBot : %.16lf --> %.16lf km\n",Top,Bot);
printf ("Search between: %.16lf ~ %.16lf km\n",RE-R[CurRegion][0],RE-R[CurRegion].back());
printf ("...\n");
cout << flush;


			// Do ray-tracing.
			const auto &v=(RayHeads[i].IsP?Vp:Vs);
			size_t radius;
			vector<double> degree;
			auto ans=RayPath(R[CurRegion],v[CurRegion],RayHeads[i].RayP,Top,Bot,degree,radius,P[CriticalAngle]);

// cout << RayHeads[i].RayP << endl;

            // If the new leg is found to be trivia, mark it for deleting.
            size_t RayLength=degree.size();
            if (RayLength==1) {D.push_back(i);continue;}


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
                        NextRegion=0;
                        RayEnd=j;
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

cout << "Next Region  (T) :" << NextRegion <<  endl;

			// For rays entering different regions, we need to calculate some new angles/ray parameters.
            double Incident=0,Rayp_td,Rayp_ts,Rayp_rd,Rayp_rs,NextPt_T,NextPr_T,NextPt_R,NextPr_R,JuncPt,JuncPr;
			bool ts=(PI[TS]==1),td=(PI[TD]==1),rd=(PI[RD]==1);
			if (RayEnd!=-1){

                NextPt_R=RayHeads[i].Pt+M*degree[RayEnd-1];
                NextPr_R=R[CurRegion][rIndex(RayEnd-1)];
                NextPt_T=RayHeads[i].Pt+M*degree[RayEnd];
                NextPr_T=R[CurRegion][rIndex(RayEnd)];

                // Re-calculate travel distance and travel time if the ray ends pre-maturely.
				ans.first.first=ans.first.second=0;
				for (int j=0;j<RayEnd-1;++j){
					double dist=sqrt( pow(R[CurRegion][rIndex(j)],2) + pow(R[CurRegion][rIndex(j+1)],2)
                                      -2*R[CurRegion][rIndex(j)]*R[CurRegion][rIndex(j+1)]*cos(M_PI/180*(degree[j+1]-degree[j])) );
					ans.first.second+=dist;
					ans.first.first+=dist/v[CurRegion][rIndex(j)];
				}

                // Find the interception between the original ray (RayEnd-1 ~ RayEnd) and the some polygon boundary.
                pair<double,double> p2={RayHeads[i].Pt+M*degree[RayEnd-1],R[CurRegion][rIndex(RayEnd-1)]};
                pair<double,double> q2={RayHeads[i].Pt+M*degree[RayEnd],R[CurRegion][rIndex(RayEnd)]};
// printf("%.15lf,%.15lf,%.15lf,%.15lf\n",p2.first,p2.second,q2.first,q2.second);

                size_t L1,L2;
                size_t SearchRegion=(NextRegion==0?CurRegion:NextRegion);

                for (L1=0;L1<Regions[SearchRegion].size();++L1){
                    L2=(L1+1)%Regions[SearchRegion].size();
                    auto res=SegmentJunction(Regions[SearchRegion][L1],Regions[SearchRegion][L2],p2,q2);
                    if (res.first) break;
                }
                if (L1==Regions[SearchRegion].size()) cout << "Error Found!" << endl;

                // Find the junction point between ray and polygon boundary.
                const pair<double,double> &p1=Regions[SearchRegion][L1],&q1=Regions[SearchRegion][L2];
// printf("%.15lf,%.15lf,%.15lf,%.15lf\n",p1.first,p1.second,q1.first,q1.second);

                double s1=(q1.second-p1.second)/(q1.first-p1.first),s2=(q2.second-p2.second)/(q2.first-p2.first);
                auto junc=LineJunction(p1,s1,p2,s2);
                JuncPt=junc.second.first;
                JuncPr=junc.second.second;
                double TiltAngle=180/M_PI*atan2(q1.second-p1.second,(q1.first-p1.first)*M_PI/180*JuncPr);

// printf("Juncture at: %.15lf,%.15lf.\n",JuncPt,JuncPr);

                // Calculate incident angle on the polygon boundary.
                double dx1_p=p1.first-JuncPt,dy1_p=p1.second-JuncPr;
                dx1_p*=2*M_PI*JuncPr/360;
                double dl1_p=sqrt(dx1_p*dx1_p+dy1_p*dy1_p);

                double dx1_q=q1.first-JuncPt,dy1_q=q1.second-JuncPr;
                dx1_q*=2*M_PI*JuncPr/360;
                double dl1_q=sqrt(dx1_q*dx1_q+dy1_q*dy1_q);

                double dx1,dy1,dl1;
                if (dl1_p>=dl1_q) dx1=dx1_p,dy1=dy1_p,dl1=dl1_p;
                else dx1=dx1_q,dy1=dy1_q,dl1=dl1_q;


                double dx2_p=p2.first-JuncPt,dy2_p=p2.second-JuncPr;
                dx2_p*=2*M_PI*JuncPr/360;
                double dl2_p=sqrt(dx2_p*dx2_p+dy2_p*dy2_p);

                double dx2_q=q2.first-JuncPt,dy2_q=q2.second-JuncPr;
                dx2_q*=2*M_PI*JuncPr/360;
                double dl2_q=sqrt(dx2_q*dx2_q+dy2_q*dy2_q);

                // OCD on travel times.
                ans.first.second+=dl2_p;
                ans.first.first+=dl2_p/v[CurRegion][rIndex(RayEnd-1)];

                double dx2,dy2,dl2;
                if (dl2_p>=dl2_q) dx2=dx2_p,dy2=dy2_p,dl2=dl2_p;
                else dx2=dx2_q,dy2=dy2_q,dl2=dl2_q;

printf("d1: %.15lf,%.15lf.\n",dx1,dy1);
printf("d2: %.15lf,%.15lf.\n",dx2,dy2);
                Incident=180/M_PI*acos((dx1*dx2+dy1*dy2)/dl1/dl2);
                if (Incident>90) Incident=180-Incident;
                Incident=90-Incident;

printf ("Incident angle :%.15lf deg\n",Incident);

                // Find the new takeoff angle and ray parameter for rays entering another region.

                double c1,c2;

                /// A. transmission to a same wave type.
                if (RayHeads[i].IsP) {
                    c1=Vp[CurRegion][rIndex(RayEnd)];
                    c2=c1/dV[CurRegion].first*dV[NextRegion].first;
                }
                else {
// cout << "CurRegion: "  << CurRegion << " " << "NextRegion " << NextRegion << endl;
                    c1=Vs[CurRegion][rIndex(RayEnd)];
                    c2=c1/dV[CurRegion].second*dV[NextRegion].second;
                }
// printf("c1/c2     : %.15lf,%.15lf\n",c1,c2);
                double Takeoff_ts=asin(sin(Incident*M_PI/180)*c2/c1)*180/M_PI;
// printf("TakeOff angle1: %.15lf deg\n",Takeoff_ts);

// printf("RayP before      : %.15lf\n",RayHeads[i].RayP);
//                         double Rayp_rb=M_PI/180*JuncPr*sin(fabs(Incident)*M_PI/180)/Rvs(JuncPr);
// printf("RayP right before: %.15lf\n",Rayp_rb);

                //// take 2D structure shape into consideration.

// printf("Tilt angle    : %.15lf deg\n",TiltAngle);
                if (0<TiltAngle && TiltAngle<90) Takeoff_ts+=M*TiltAngle;
                else if (90<=TiltAngle && TiltAngle<=180) Takeoff_ts-=M*(180-TiltAngle);
                else if (-90<=TiltAngle && TiltAngle<=0) Takeoff_ts+=M*TiltAngle;
                else Takeoff_ts+=M*(180+TiltAngle);

                if (::isnan(Takeoff_ts)) ts=false;
// printf("TakeOff angle2: %.15lf deg\n",Takeoff_ts);
// printf("c2            : %.15lf km/sec.\n",c2);
                Rayp_ts=M_PI/180*NextPr_T*sin(fabs(Takeoff_ts)*M_PI/180)/c2;
// printf("RayP after       : %.15lf\n",Rayp_ts);


                /// D. reflection to a same wave type.
                double Takeoff_rs=(RayHeads[i].GoUp?Incident:180-Incident);
printf("TakeOff angle1 (R): %.15lf deg\n",Takeoff_rs);
printf("Tilt angle     (R): %.15lf deg\n",TiltAngle);
                if (0<TiltAngle && TiltAngle<90) Takeoff_rs+=M*TiltAngle;
                else if (90<=TiltAngle && TiltAngle<=180) Takeoff_rs-=M*(180-TiltAngle);
                else if (-90<=TiltAngle && TiltAngle<=0) Takeoff_rs+=M*TiltAngle;
                else Takeoff_rs+=M*(180+TiltAngle);
printf("TakeOff angle2 (R): %.15lf deg\n",Takeoff_rs);
                Rayp_rs=M_PI/180*NextPr_R*sin(fabs(Takeoff_rs)*M_PI/180)/(RayHeads[i].IsP?Vp[CurRegion][rIndex(RayEnd-1)]:Vs[CurRegion][rIndex(RayEnd-1)]);
printf("RayP           (R): %.15lf deg\n",Rayp_rs);

			}
            else { // If ray doesn't end pre-maturelly, ray parameters shouldn't change.
				RayEnd=degree.size();
				NextRegion=CurRegion;
                NextPt_T=NextPt_R=RayHeads[i].Pt+M*degree[RayEnd-1];
                NextPr_T=NextPr_R=R[CurRegion][rIndex(RayEnd-1)];
                Rayp_td=Rayp_ts=Rayp_rd=Rayp_rs=RayHeads[i].RayP;
                Incident=asin(RayHeads[i].RayP*180/M_PI*v[CurRegion][rIndex(RayEnd-1)]/R[CurRegion][rIndex(RayEnd-1)])*180/M_PI;
// printf ("Incident angle :%.15lf deg\n",Incident);
            }

cout << "RayEnd at Index  :" << RayEnd-1 << " (inclusive) / " << RayLength << endl;
printf ("RayEnd at (R)    :%.15lf deg, %.15lf (inclusive) km\n",NextPt_R,RE-NextPr_R);
printf ("RayEnd at (T)    :%.15lf deg, %.15lf (inclusive) km\n",NextPt_T,RE-NextPr_T);


			// Output valid part ray paths.
			ofstream fpout;
			fpout.open(PS[OutFilePrefix]+to_string(i+1),ofstream::app);
			fpout << "> " << (RayHeads[i].IsP?"P":"S") << endl;
			for (int j=0;j<RayEnd;++j)
				fpout << RayHeads[i].Pt+M*degree[j] << " " << R[CurRegion][rIndex(j)] << '\n';
			fpout.close();


			// Update current RayHead.
			RayHeads[i].TravelTime+=ans.first.first;
			RayHeads[i].TravelDist+=ans.first.second;


            // Add rules of: (t)ransmission/refrection and (r)eflection to (s)ame or (d)ifferent way type.
            // Notice reflection with the same wave type is always allowed("rs" is always true; actually, updated ray head will be rs).

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
                RayHeads.back().Pt=NextPt_T;
                RayHeads.back().Pr=NextPr_T;
                RayHeads.back().RayP=Rayp_ts;
                RayHeads.back().InRegion=NextRegion;
            }

            if (td) {
                RayHeads.push_back(RayHeads[i]);
                RayHeads.back().IsP=!RayHeads.back().IsP;
                RayHeads[i].RayP=Rayp_td;
                RayHeads.back().InRegion=NextRegion;
            }

            if (rd) {
                RayHeads.push_back(RayHeads[i]);
                RayHeads.back().IsP=!RayHeads.back().IsP;
                RayHeads[i].RayP=Rayp_rd;
                RayHeads.back().InRegion=CurRegion;
            }

            // rs
            RayHeads[i].RayP=Rayp_rs;
            RayHeads[i].GoUp=!RayHeads[i].GoUp;
            RayHeads[i].Pt=NextPt_R;
            RayHeads[i].Pr=NextPr_R;

		} // End of ray loop.

        // Delete trivia rays and update number of valid rays.
        for (auto rit=D.rbegin();rit!=D.rend();++rit)
            RayHeads.erase(RayHeads.begin()+*rit);
        ValidN-=D.size();

    } // End of leg loop.

	ofstream fpout(PS[OutInfoFile]);
	for (size_t i=0;i<ValidN;++i)
		fpout << i+1 << " " << RayHeads[i].Color << " " << RayHeads[i].TravelTime << " sec. " << RayHeads[i].Pt << " deg. "
              << RayHeads[i].TravelDist << "km." << endl;
	fpout.close();

	return 0;
}
