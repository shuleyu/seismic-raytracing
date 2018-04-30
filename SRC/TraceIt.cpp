#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<string>
#include<iterator>
#include<algorithm>
#include<iomanip>
#include<map>
#include<complex>

#include<CreateGrid.hpp>
#include<LineJunction.hpp>
#include<LocDist.hpp>
#include<Lon2180.hpp>
#include<Lon2360.hpp>
#include<PointInPolygon.hpp>
#include<PREM.hpp>
#include<RayPath.hpp>
#include<SegmentJunction.hpp>
#include<PlaneWaveCoefficients.hpp>

using namespace std;

// Define the ray head structure.
class Ray {
    public:
        bool IsP,GoUp,GoLeft;
        string Color,Comp,Debug;
        int InRegion,Prev,RemainingLegs,Surfacing;
        double Pt,Pr,TravelTime,TravelDist,RayP,Amp,Inc,Takeoff;

        Ray(bool p, bool g, bool l, string c, string cmp, int i,int rl, double th, double r, double t, double d, double rp,double to) :
            IsP(p), GoUp(g), GoLeft(l), Color(c), Comp(cmp), Debug(""), InRegion(i), Prev(-1), RemainingLegs(rl), Surfacing(0), Pt(th),Pr(r),TravelTime(t),TravelDist(d), RayP(rp), Amp(1),Inc(0), Takeoff(to) {}
};


// Utilities for 1D-altering the PREM model.
vector<double> MakeRef(const double &depth,const vector<vector<double>> &dev){
    double rho=Drho(depth),vs=Dvs(depth),vp=Dvp(depth);
    for (size_t i=0;i<dev.size();++i){
        if (depth>=dev[i][0]) {
            vp*=(1+dev[i][1]/100);
            vs*=(1+dev[i][2]/100);
            rho*=(1+dev[i][3]/100);
            break;
        }
    }
    return {vp,vs,rho};
}


int main(int argc, char **argv){

    enum PIenum{DebugInfo,TS,TD,RD,FLAG1};
    enum PSenum{InputRays,Layers,Depths,Ref,Polygons,OutFilePrefix,ReceiverFile,PolygonOutPrefix,FLAG2};
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


    // Create 1D reference layers. (R[0]. 0 means 1D reference model)
    const double RE=6371.0;
    vector<vector<double>> R{vector<double> ()};
    double prev_depth,depth,inc,next_inc,MinInc;

    ifstream fpin(PS[Layers]);
    fpin >> prev_depth >> inc;
    MinInc=inc; // Record the minimum layer increment for later use. (Solve round-off error related issues)

    while (fpin >> depth >> next_inc){

        auto tmpr=CreateGrid(RE-depth+1e-6,RE-prev_depth+1e-6,inc,2); // +1e-6 here to solve round-off error related issues;
                                                                      // this will also affects our choice when calculating coefficients.
                                                                      // See the definition of variable "si".

        if (!R[0].empty()) R[0].pop_back();
        R[0].insert(R[0].end(),make_move_iterator(tmpr.rbegin()),make_move_iterator(tmpr.rend()));
        if (next_inc<=0) break;

        prev_depth=depth;
        inc=next_inc;
        MinInc=min(MinInc,inc);
    }
    fpin.close();


    // Read in 1D reference deviation.
    double dvp,dvs,drho;
    vector<vector<double>> Deviation;

    fpin.open(PS[Ref]);
    while (fpin >> depth >> dvp >> dvs >> drho) Deviation.push_back({depth,dvp,dvs,drho});
    auto cmp=[](const vector<double> &p1, const vector<double> &p2){return p1[0]>p2[0];};
    sort(Deviation.begin(),Deviation.end(),cmp); // sort it to descending order (this property is used in MakeRef).
    fpin.close();


    // Create 1D reference model properties. (Vp[0],Vs[0] and Rho[0])
    vector<vector<double>> Vp{vector<double> ()},Vs=Vp,Rho=Vp;
    Vp[0].resize(R[0].size());
    Vs[0].resize(R[0].size());
    Rho[0].resize(R[0].size());
    for (size_t i=0;i<R[0].size();++i){
        auto ans=MakeRef(RE-R[0][i],Deviation);
        Vp[0][i]=ans[0];
        Vs[0][i]=ans[1];
        Rho[0][i]=ans[2];
    }


    // Read in special depths (force reflection & refraction at these depths).
    vector<double> SpecialDepths;
    fpin.open(PS[Depths]);
    while (fpin >> depth) SpecialDepths.push_back(depth);
    sort(SpecialDepths.begin(),SpecialDepths.end()); // sort it ascending.
    fpin.close();


    // Read in 2D structure polygons.
    // Create layers and properties for each polygons.
    //
    // The layers of these polygons are derived from the layers of the 1D refernce.
    // The boundaries of these derived layers are slightly larger than the polygon vertical limit, for convenience.
    string tmpstr;
    vector<pair<double,double>> tmpregion,tmpregion2;
    vector<vector<pair<double,double>>> Regions{tmpregion};
    vector<double> dVp{1},dVs{1},dRho{1};

    fpin.open(PS[Polygons]);
    while (getline(fpin,tmpstr)){

        if (tmpstr.empty()) continue;

        char c;
        stringstream ss(tmpstr);
        double theta,depth,r,minr,maxr;

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

                    ofstream fpout(PS[PolygonOutPrefix]+to_string(Regions.size()-1));
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
                dVp.push_back(1.0+dvp/100);
                dVs.push_back(1.0+dvs/100);
                dRho.push_back(1.0+drho/100);

                Rho.push_back(vector<double> ());
                Vp.push_back(vector<double> ());
                Vs.push_back(vector<double> ());
                for (const auto &item:R.back()) {
                    auto ans=MakeRef(RE-item,Deviation);
                    Vp.back().push_back(dVp.back()*ans[0]);
                    Vs.back().push_back(dVs.back()*ans[1]);
                    Rho.back().push_back(dRho.back()*ans[2]);
                }
            }

            // Read next region. Initalizations.
            ss >> c >> dvp >> dvs >> drho;
            minr=RE+1,maxr=-1;
            tmpregion.clear();

            // The last polygon is not recorded, therefore a ">" is needed at the end of polygon inputs (see INFILE).
        }
        else {
            // Read next position of current polygon.
            ss >> theta >> depth;
            r=RE-depth;
            minr=min(minr,r);
            maxr=max(maxr,r);
            tmpregion.push_back({theta,r});
        }
    }
    fpin.close();


    // For plotting: 1D reference property deviation depths.
    ofstream fpout;
    if (!Deviation.empty()){
        fpout.open(PS[PolygonOutPrefix]+"0");
        for (auto &item:Deviation) {
            fpout << ">\n";
            double d=RE-item[0];
            for (double t=0;t<360;t=t+0.1) fpout << t << " " << d << '\n';
        }
        fpout.close();
    }


    // Read in initial rays.
    int steps;
    string color,phase;
    double theta,takeoff;
    vector<Ray> RayHeads;

    fpin.open(PS[InputRays]);
    while (fpin >> theta >> depth >> takeoff >> phase >> color >> steps){

        // Source in any polygons?
        int rid=0;
        for (size_t i=1;i<Regions.size();++i)
            if (PointInPolygon(Regions[i],make_pair(theta,RE-depth),1)) {rid=i;break;}

        // Calculate ray parameter.
        auto ans=MakeRef(depth,Deviation);
        double rayp=M_PI/180*(RE-depth)*sin(fabs(takeoff)/180*M_PI)/(phase[0]=='P'?ans[0]*dVp[rid]:ans[1]*dVs[rid]);

        // Push this ray into "RayHeads" for future processing.
        RayHeads.push_back(Ray(phase[0]=='P',fabs(takeoff)>=90,takeoff<0,color,phase,rid,steps,theta,RE-depth,0,0,rayp,takeoff));
    }
    fpin.close();

    fpout.open(PS[ReceiverFile]);
    fpout << "<TextHeight> <Takeoff> <Dist> <TravelTime> <DispAmp> <Incident> <RemainingLegs> <WaveTypeTrain> <RayTrain>" << '\n';
    fpout.close();


    // Start ray tracing. (Woohoo!)
    //
    // Will process each "Ray" leg in "RayHeads".
    // For future legs generated by reflction/refraction, create new "Ray" and push to the end of "RayHeads" queue.
    double PlotPosition=1;
    map<string,double> PlotColorPosition;

    size_t QueueL=0,QueueR=RayHeads.size(); // Rays between [QueueL,QueueR) needs processing.

    while (QueueL!=QueueR) {

        for (size_t i=QueueL;i<QueueR;++i){

            if (RayHeads[i].RemainingLegs==0) continue;

            // Locate the begining and ending depths for the next leg.

            /// ... among special depths.

            //// Which special depth is cloest to ray head depth?
            double RayHeadDepth=RE-RayHeads[i].Pr,NextDepth,MinDiff=RE+1;
            size_t Cloest=0;
            for (size_t j=0;j<SpecialDepths.size();++j) {
                double CurDiff=fabs(SpecialDepths[j]-RayHeadDepth);
                if (MinDiff>CurDiff) {MinDiff=CurDiff;Cloest=j;}
            }

            //// Next depth should be the cloest special depth, but more things needs to be considered:
            //// Is the ray going up or down? Is the ray head already at the cloest special depth?
            NextDepth=SpecialDepths[Cloest];
            if ( RayHeads[i].GoUp && (SpecialDepths[Cloest]>RayHeadDepth || RayHeadDepth-SpecialDepths[Cloest]<MinInc) )
                NextDepth=SpecialDepths[Cloest-1];
            else if ( !RayHeads[i].GoUp && (SpecialDepths[Cloest]<RayHeadDepth || SpecialDepths[Cloest]-RayHeadDepth<MinInc) )
                NextDepth=SpecialDepths[Cloest+1];

            double Top=min(RayHeadDepth,NextDepth),Bot=max(RayHeadDepth,NextDepth);

            /// ... among current 2D "Regions" vertical limits.
            int CurRegion=RayHeads[i].InRegion;
            Top=max(Top,RE-R[CurRegion][0]);
            Bot=min(Bot,RE-R[CurRegion].back());


            // Print some debug info.
            if (PI[DebugInfo]) {
                RayHeads[i].Debug+=to_string(1+i)+" --> ";
                cout << '\n' << "Calculating   : " << RayHeads[i].Debug;
                printf ("\nLocation      : %.15lf deg, %.15lf km (In region: %d)\n", RayHeads[i].Pt, RE-RayHeads[i].Pr,CurRegion);
                cout << "Will go as    : " << (RayHeads[i].IsP?"P, ":"S, ") << (RayHeads[i].GoUp?"Up, ":"Down, ")
                     << (RayHeads[i].GoLeft?"Left":"Right") << endl;
                printf ("Start, End (D): %.16lf --> %.16lf km\n",RayHeadDepth,NextDepth);
                printf ("Actual s/e (D): %.16lf --> %.16lf km\n",Top,Bot);
                printf ("Search between: %.16lf ~ %.16lf km with rayp  : %.16lf sec/deg\n",RE-R[CurRegion][0],RE-R[CurRegion].back(),RayHeads[i].RayP);
                printf ("\n...\n\n");
                cout << flush;
            }


            // Use ray-tracing code "RayPath".
            size_t radius;
            vector<double> degree;
            const auto &v=(RayHeads[i].IsP?Vp:Vs);
            auto ans=RayPath(R[CurRegion],v[CurRegion],RayHeads[i].RayP,Top,Bot,degree,radius,P[CriticalAngle]);


            // If the new leg is trivia, no further operation needed.
            size_t RayLength=degree.size();
            if (RayLength==1) {RayHeads[i].RemainingLegs=0;continue;}


            // If the new leg is a reflection of down-going S to up-going P, and also the new leg turns, also mark it as invalid.
            int PrevID=RayHeads[i].Prev;
            if (PrevID!=-1 && !RayHeads[PrevID].GoUp && !RayHeads[PrevID].IsP && RayHeads[i].GoUp && RayHeads[i].IsP && ans.second) {
                RayHeads[i].RemainingLegs=0;
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


            // Follow the new ray path to see if the new leg enters another region.
            int RayEnd=-1,NextRegion=-1,M=(RayHeads[i].GoLeft?-1:1);

            for (size_t j=0;j<degree.size();++j){

                pair<double,double> p={RayHeads[i].Pt+M*degree[j],R[CurRegion][rIndex(j)]}; // point location on new leg.

                if (CurRegion!=0){ // New leg starts in some 2D polygon ...

                    if (PointInPolygon(Regions[CurRegion],p)) continue; // ... and this point stays in that polygon.
                    else { // ... but this point enters another polygon.

                        RayEnd=j;

                        // which region is the new leg entering?
                        for (size_t k=1;k<Regions.size();++k){
                            if (PointInPolygon(Regions[k],p)) {
                                NextRegion=k;
                                break;
                            }
                        }
                        if (NextRegion==-1) NextRegion=0; // if can't find next 2D polygons, it must return to the 1D reference region.
                        break;
                    }
                }
                else { // New leg staats in 1D reference region. Search for the region it enters.
                    for (size_t k=1;k<Regions.size();++k){
                        if (PointInPolygon(Regions[k],p)) { // If ray enters another region.
                            RayEnd=j;
                            NextRegion=k;
                            break;
                        }
                    }
                    if (RayEnd!=-1) break; // If ray enters another region.
                }
            }


            // Print some debug info.
            if (PI[DebugInfo]) {
                if (NextRegion==-1) cout << "No region change." <<  endl;
                else cout << "Next Region is    : " << NextRegion << " (for refraction)" << endl;
            }


            // Prepare reflection/refraction flags. (notice "rs" [r]eflection to [s]ame wave type is always possible)
            bool ts=(PI[TS]==1),td=(PI[TD]==1 && RayHeads[i].Comp!="SH"),rd=(PI[RD]==1 && RayHeads[i].Comp!="SH");


            // Locate the end of new leg, which is needed to calculate incident angle, coefficients, next ray parameter, etc.
            // (if interface is not horizontal("TiltAngle"), ray parameter will change.)
            //
            // Decision made: If the last line segment of the new leg crosses interface(at "JuncPt/JuncPr"),
            // for reflection: the end point outside of new region ("NextP?_R") is the next ray starting point;
            // for refraction: the end point inside of the new region ("NextP?_T") is the next ray starting piont.
            double NextPt_R,NextPr_R,NextPt_T,NextPr_T,JuncPt,JuncPr,TiltAngle,Rayp_td=-1,Rayp_ts=-1,Rayp_rd=-1,Rayp_rs=-1;

            pair<double,double> p2,q2; // Two end points of the last line segment of the new leg.
                                       // Notice, for normal rays hit the horizontal interface, one end point is on the interface.

            if (RayEnd!=-1){ // If ray ends pre-maturely (last line segment crossing the interface).

                // For reflection, the future rays start from the last point in the current region (index: RayEnd-1).
                NextPt_R=RayHeads[i].Pt+M*degree[RayEnd-1];
                NextPr_R=R[CurRegion][rIndex(RayEnd-1)];


                // For transmission/refraction, the futuer rays start from the first point in the next region (index: RayEnd).
                NextPt_T=RayHeads[i].Pt+M*degree[RayEnd];
                NextPr_T=R[CurRegion][rIndex(RayEnd)];


                // Re-calculate travel distance and travel time till the last point in the current region.
                ans.first.first=ans.first.second=0;
                for (int j=0;j<RayEnd-1;++j){
                    double dist=sqrt( pow(R[CurRegion][rIndex(j)],2) + pow(R[CurRegion][rIndex(j+1)],2)
                                      -2*R[CurRegion][rIndex(j)]*R[CurRegion][rIndex(j+1)]*cos(M_PI/180*(degree[j+1]-degree[j])) );
                    ans.first.second+=dist;
                    ans.first.first+=dist/v[CurRegion][rIndex(j+1)];
                }


                // Find the junction between the last line segment (index: RayEnd-1 ~ RayEnd) and polygon boundary segment (index: L1 ~ L2).
                size_t L1,L2,SearchRegion=(NextRegion==0?CurRegion:NextRegion);
                p2={NextPt_R,NextPr_R};
                q2={NextPt_T,NextPr_T};

                for (L1=0;L1<Regions[SearchRegion].size();++L1){ // Search around the polygon.
                    L2=(L1+1)%Regions[SearchRegion].size();
                    auto res=SegmentJunction(Regions[SearchRegion][L1],Regions[SearchRegion][L2],p2,q2);
                    if (res.first) break;
                }
                if (L1==Regions[SearchRegion].size()) {
                    cerr << "!!!!!!!!!!! Can't find junction! Bugs here !!!!!!!!!!!" << endl;
                    return 1;
                }


                // Find the junction point between ray and polygon boundary.
                const pair<double,double> &p1=Regions[SearchRegion][L1],&q1=Regions[SearchRegion][L2];

                double s1=(q1.second-p1.second)/(q1.first-p1.first),s2=(q2.second-p2.second)/(q2.first-p2.first);
                auto junc=LineJunction(p1,s1,p2,s2);
                JuncPt=junc.second.first;
                JuncPr=junc.second.second;


                // Print some debug info.
                if (PI[DebugInfo]) printf("Junction at       : %.15lf,%.15lf.\n",JuncPt,JuncPr);


                // Twick travel times and travel distance, compensate for the lost part.
                double dlx=(p2.first-JuncPt)*M_PI*JuncPr/180,dly=p2.second-JuncPr;
                double dl=sqrt(dlx*dlx+dly*dly);
				ans.first.second+=dl;
				ans.first.first+=dl/v[CurRegion][rIndex(RayEnd-1)]; // Use the velocit within current region to avoid possible "inf" travel time.


                // Get the geometry of the boundary.
                TiltAngle=180/M_PI*atan2(q1.second-p1.second,(q1.first-p1.first)*M_PI/180*JuncPr);

            }
            else { // If ray doesn't end pre-maturelly (stays in the same region and reflect/refract on horizontal intervals)
                   // (one end point of the last line segment (index: RayEnd-1) is on the interface)

                RayEnd=degree.size();
                NextRegion=CurRegion;

                // Get futuer rays starting point.
                NextPt_T=NextPt_R=RayHeads[i].Pt+M*degree[RayEnd-1];
                NextPr_T=NextPr_R=R[CurRegion][rIndex(RayEnd-1)];


                // Notice ray parameter doesn't change if reflection/refraction interface is horizontal.
                Rayp_td=Rayp_ts=Rayp_rd=Rayp_rs=RayHeads[i].RayP;


                // Get the last segment of the new leg.
                p2={RayHeads[i].Pt+M*degree[RayEnd-2],R[CurRegion][rIndex(RayEnd-2)]};
                q2={NextPt_T,NextPr_T};

                // Get the geometry of the boundary.
                TiltAngle=0;
                JuncPt=NextPt_T;
                JuncPr=NextPr_T;

            } // End of dealing with rays entering another region.


            // Get the geometry of the last section. (Ray direction: "Rayd" [-180 ~ 180])
            double Rayd=180/M_PI*atan2(q2.second-p2.second,(q2.first-p2.first)*M_PI/180*JuncPr);


            // Calculate incident angle. (the acute angle between ray direction and the normal of interface)
            double Rayd_Hor=Lon2360(Rayd-TiltAngle),Incident=fabs(Lon2180(Rayd_Hor));
            Incident=(Incident>90?Incident-90:90-Incident);


            // Print some debug info.
            if (PI[DebugInfo]) {
                printf("Incident angle    : %.15lf deg\n",Incident);
                printf("Ray direction     : %.15lf deg\n",Rayd);
                printf("Interface tilt    : %.15lf deg\n",TiltAngle);
            }


            // Prepare to calculate reflection/refractoin(transmission) coefficients.
            string Mode,AdaptedMode,Polarity=(RayHeads[i].Comp=="SH"?"SH":"PSV");
            if (fabs(RE-NextPr_R)<MinInc) Mode="SA"; // At the surface.
            else if (fabs(3480-NextPr_R)<MinInc) Mode=(RayHeads[i].GoUp?"LS":"SL"); // At the CMB.
            else if (fabs(1221.5-NextPr_R)<MinInc) Mode=(RayHeads[i].GoUp?"SL":"LS"); // At the ICB.
            else if (1221.5<NextPr_R && NextPr_R<3480) Mode="LL";
            else Mode="SS"; // AdaptedMode is for 2D regions with extreme S velocity drop.

            double rho1,vp1,vs1,rho2,vp2,vs2,c1,c2;;

            if (CurRegion!=NextRegion) { // if ray enters a different region.
                rho1=Rho[CurRegion][rIndex(RayEnd-1)];
                rho2=Rho[CurRegion][rIndex(RayEnd)]/dRho[CurRegion]*dRho[NextRegion];
                vp1=Vp[CurRegion][rIndex(RayEnd-1)];
                vp2=Vp[CurRegion][rIndex(RayEnd)]/dVp[CurRegion]*dVp[NextRegion];
                vs1=Vs[CurRegion][rIndex(RayEnd-1)];
                vs2=Vs[CurRegion][rIndex(RayEnd)]/dVs[CurRegion]*dVs[NextRegion];
            }
            else { // if ray stays in the same region. ray ends normally.
                int si=RayEnd;
                if (RayHeads[i].GoUp) --si;

                rho1=Rho[CurRegion][rIndex(si-1)];
                rho2=Rho[CurRegion][rIndex(si)];
                vp1=Vp[CurRegion][rIndex(si-1)];
                vp2=Vp[CurRegion][rIndex(si)];
                vs1=Vs[CurRegion][rIndex(si-1)];
                vs2=Vs[CurRegion][rIndex(si)];
            }
            if (vs1<0.01) Mode[0]='L';
            if (vs2<0.01 && Mode[1]=='S') Mode[1]='L';

            // Print some debug info.
            if (PI[DebugInfo]) printf("rho1/vp1/vs1      : %.15lf g/cm3, %.15lf km/sec, %.15lf km/sec\n",rho1,vp1,vs1);
            if (PI[DebugInfo]) printf("rho2/vp2/vs2      : %.15lf g/cm3, %.15lf km/sec, %.15lf km/sec\n",rho2,vp2,vs2);

            complex<double> R_PP,R_PS,R_SP,R_SS,T_PP,T_PS,T_SP,T_SS;
            R_PP=R_PS=R_SP=R_SS=T_PP=T_PS=T_SP=T_SS=0;

            /// A. Refractions/Transmissions to the same wave type.

            //// Coefficients. (T_PP,T_SS)
            auto Coef=PlaneWaveCoefficients(rho1,vp1,vs1,rho2,vp2,vs2,Incident,Polarity,Mode);
            if (Mode=="SS") {
                if (RayHeads[i].Comp=="SH") T_SS=Coef[1];
                else {T_PP=Coef[4];T_SS=Coef[7];}
            }
            if (Mode=="SL" && RayHeads[i].Comp=="P") T_PP=Coef[4];
            if (Mode=="LS" && RayHeads[i].Comp=="P") T_PP=Coef[1];
            if (Mode=="LL" && RayHeads[i].Comp=="P") T_PP=Coef[1];

            if (RayHeads[i].IsP) ts&=(T_PP.imag()==0);
            else ts&=(T_SS.imag()==0);

            //// take-off angles.
            if (RayHeads[i].IsP) {c1=vp1;c2=vp2;}
            else {c1=vs1;c2=vs2;}

            //// difference between incident angle and takeoff angle.
            double Takeoff_ts=asin(sin(Incident*M_PI/180)*c2/c1)*180/M_PI-Incident;
            ts&=!std::isnan(Takeoff_ts);

            //// take 2D structure shape into consideration.
            if ((0<Rayd_Hor && Rayd_Hor<=90) || (180<Rayd_Hor && Rayd_Hor<=270)) Takeoff_ts=Lon2180(Rayd_Hor-Takeoff_ts+TiltAngle+90);
            else Takeoff_ts=Lon2180(Rayd_Hor+Takeoff_ts+TiltAngle+90);

            //// new ray paramter.
            if (CurRegion!=NextRegion) Rayp_ts=M_PI/180*NextPr_T*sin(fabs(Takeoff_ts)*M_PI/180)/c2;
            ts&=!std::isnan(Rayp_ts);


            /// B. Refractions/Transmissions to different wave type.

            //// Coefficients. (T_PS,T_SP)
            if (Mode=="SS" && RayHeads[i].Comp!="SH") {T_PS=Coef[1];T_SP=Coef[6];}
            if (Mode=="SL" && RayHeads[i].Comp=="SV") T_SP=Coef[5];
            if (Mode=="LS" && RayHeads[i].Comp=="P") T_PS=Coef[2];

            if (RayHeads[i].IsP) td&=(T_PS.imag()==0);
            else td&=(T_SP.imag()==0);

            //// take-off angles.
            if (RayHeads[i].IsP) {c1=vp1;c2=vs2;}
            else {c1=vs1;c2=vp2;}

            //// difference between incident angle and takeoff angle.
            double Takeoff_td=asin(sin(Incident*M_PI/180)*c2/c1)*180/M_PI-Incident;
            td&=!std::isnan(Takeoff_td);

            //// take 2D structure shape into consideration.
            if ((0<Rayd_Hor && Rayd_Hor<=90) || (180<Rayd_Hor && Rayd_Hor<=270)) Takeoff_td=Lon2180(Rayd_Hor-Takeoff_td+TiltAngle+90);
            else Takeoff_td=Lon2180(Rayd_Hor+Takeoff_td+TiltAngle+90);

            //// new ray paramter.
            if (CurRegion!=NextRegion) Rayp_td=M_PI/180*NextPr_T*sin(fabs(Takeoff_td)*M_PI/180)/c2;
            td&=!std::isnan(Rayp_td);


            /// C. Reflection to a different wave type.

            //// Coefficients. (R_PS,R_SP)
            if (Mode=="SS" && RayHeads[i].Comp!="SH") {R_PS=Coef[1];R_SP=Coef[2];}
            if (Mode=="SL" && RayHeads[i].Comp!="SH") {R_PS=Coef[1];R_SP=Coef[2];}
            if (Mode=="SA" && RayHeads[i].Comp!="SH") {R_PS=Coef[1];R_SP=Coef[2];}

            //// take-off angles.
            c1=vs1;c2=vp1;
            if (RayHeads[i].IsP) swap(c1,c2);

            //// difference between incident angle and takeoff angle.
            double Takeoff_rd=asin(sin(Incident*M_PI/180)*c2/c1)*180/M_PI-Incident;
            rd&=!std::isnan(Takeoff_rd);

            //// take 2D structure shape into consideration.
            double x=Lon2360(-Rayd_Hor);
            if ((0<x && x<=90) || (180<x && x<=270)) Takeoff_rd=Lon2180(x-Takeoff_rd+TiltAngle+90);
            else Takeoff_rd=Lon2180(x+Takeoff_rd+TiltAngle+90);

            //// new ray paramter.
            if (CurRegion!=NextRegion) Rayp_rd=M_PI/180*NextPr_R*sin(fabs(Takeoff_rd)*M_PI/180)/c2;
            rd&=!std::isnan(Rayp_rd);


            /// D. reflection to a same wave type.
            //// Coefficients. (R_PP,R_SS)
            if (Mode=="SS") {
                if (RayHeads[i].Comp=="SH") R_SS=Coef[0];
                else {R_PP=Coef[0];R_SS=Coef[3];}
            }
            if (Mode=="SL") {
                if (RayHeads[i].Comp=="SH") R_SS=1.0;
                else {R_PP=Coef[0];R_SS=Coef[3];}
            }
            if (Mode=="SA") {
                if (RayHeads[i].Comp=="SH") R_SS=1.0;
                else {R_PP=Coef[0];R_SS=Coef[3];}
            }
            if ((Mode=="LS" || Mode=="LL") && RayHeads[i].Comp=="P") R_PP=Coef[0];
            if (ans.second) R_SS=R_PP=1;

            //// new ray paramter.
            double Takeoff_rs=Lon2180(-Rayd_Hor+TiltAngle+90);
            if (CurRegion!=NextRegion) Rayp_rs=M_PI/180*NextPr_R*sin(fabs(Takeoff_rs)*M_PI/180)/c1;


            // Print some debug info.
            if (PI[DebugInfo]) {
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
            }


            // Update current RayHead.
            RayHeads[i].TravelTime=ans.first.first;
            RayHeads[i].TravelDist=ans.first.second;
            --RayHeads[i].RemainingLegs;
            RayHeads[i].Inc=Incident;


            // Output valid part ray paths.
            ofstream fpout;
            fpout.open(PS[OutFilePrefix]+to_string(i+1),ofstream::app);
            fpout << "> " << (RayHeads[i].Color=="black"?(RayHeads[i].IsP?"blue":"red"):RayHeads[i].Color) << " "
                          << (RayHeads[i].IsP?"P ":"S ") << RayHeads[i].TravelTime << " sec. " << RayHeads[i].Inc << " IncDeg. "
                          << RayHeads[i].Amp << " DispAmp. " << RayHeads[i].TravelDist << " km. " << endl;
            for (int j=0;j<RayEnd;++j)
                fpout << RayHeads[i].Pt+M*degree[j] << " " << R[CurRegion][rIndex(j)] << '\n';
            fpout.close();

            // If ray reaches surface, output info at the surface.
            if (fabs(NextPr_R-RE)<MinInc) ++RayHeads[i].Surfacing;
            if (fabs(NextPr_R-RE)<MinInc && RayHeads[i].Surfacing<2) {

                fpout.open(PS[ReceiverFile],ofstream::app);

                if (PlotColorPosition.find(RayHeads[i].Color)==PlotColorPosition.end())
                    PlotColorPosition[RayHeads[i].Color]=PlotPosition++;
                fpout << PlotColorPosition[RayHeads[i].Color] << " ";

                // Accumulate the travel-time.
                int I=i;
                double tt=0;
                vector<int> hh;
                while (I!=-1) {
                    hh.push_back(I);
                    tt+=RayHeads[I].TravelTime;
                    I=RayHeads[I].Prev;
                }
                fpout << RayHeads[hh.back()].Takeoff << " " << NextPt_R << " "
                      << tt << " " << RayHeads[i].Amp << " " << RayHeads[i].Inc << " " << RayHeads[i].RemainingLegs << " ";
                for (auto rit=hh.rbegin();rit!=hh.rend();++rit) fpout << (RayHeads[*rit].IsP?"P":"S") << ((*rit)==*hh.begin()?" ":"->");
                for (auto rit=hh.rbegin();rit!=hh.rend();++rit) fpout << (1+*rit) << ((*rit)==*hh.begin()?"\n":"->");
                fpout.close();
                continue;
            }


            // Add rules of: (t)ransmission/refrection and (r)eflection to (s)ame or (d)ifferent way type.
            // Notice reflection with the same wave type is always allowed. ("rs" is always true)

            /// if ray going down and turns.
            if (!RayHeads[i].GoUp && ans.second) ts=td=rd=false;

            /// if ray ends at the surface.
            if (fabs(RE-NextPr_R)<MinInc) ts=td=false;

            /// if ray goes down to CMB, no transmission S; if ray goes up to ICB, not transmission S.
            if ((!RayHeads[i].GoUp && fabs(NextPr_T-3480)<MinInc) || (RayHeads[i].GoUp && fabs(NextPr_T-1221.5)<MinInc)) {
                ts&=RayHeads[i].IsP;
                td&=(!RayHeads[i].IsP);
            }

            /// if ray goes down to ICB as P, no reflection S.
            if (!RayHeads[i].GoUp && RayHeads[i].IsP && fabs(NextPr_R-1221.5)<MinInc) rd=false;

            /// if ray goes up to CMB as P, no reflection S.
            if (RayHeads[i].GoUp && RayHeads[i].IsP && fabs(NextPr_R-3480)<MinInc) rd=false;


            // Add new ray heads to "RayHeads" according to the rules ans reflection/refraction angle calculation results.

            if (ts) {
                RayHeads.push_back(RayHeads[i]);
                RayHeads.back().Prev=i;
                RayHeads.back().Pt=NextPt_T;
                RayHeads.back().Pr=NextPr_T;
                RayHeads.back().RayP=Rayp_ts;
                RayHeads.back().GoUp=(fabs(Takeoff_ts)>90);
                RayHeads.back().GoLeft=(Takeoff_ts<0);
                RayHeads.back().InRegion=NextRegion;
                RayHeads.back().Amp*=(RayHeads.back().IsP?T_PP.real():T_SS.real());
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
                RayHeads.back().Amp*=(RayHeads.back().IsP?T_PS.real():T_SP.real());
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
                double sign1=(R_PS.imag()==0?(R_PS.real()<0?-1:1):1);
                double sign2=(R_SP.imag()==0?(R_SP.real()<0?-1:1):1);
                RayHeads.back().Amp*=(RayHeads.back().IsP?(sign1*abs(R_PS)):(sign2*abs(R_SP)));
            }

            // rs is always possible.
            RayHeads.push_back(RayHeads[i]);
            RayHeads.back().Prev=i;
            RayHeads.back().Pt=NextPt_R;
            RayHeads.back().Pr=NextPr_R;
            RayHeads.back().RayP=Rayp_rs;
            RayHeads.back().GoUp=(fabs(Takeoff_rs)>90);
            RayHeads.back().GoLeft=(Takeoff_rs<0);
            double sign1=(R_PP.imag()==0?(R_PP.real()<0?-1:1):1);
            double sign2=(R_SS.imag()==0?(R_SS.real()<0?-1:1):1);
            RayHeads.back().Amp*=(RayHeads.back().IsP?(sign1*abs(R_PP)):(sign2*abs(R_SS)));

        } // End of processing Rays between [QueueL,QueueR)

        QueueL=QueueR;
        QueueR=RayHeads.size();

    } // End of ray tracing.

    return 0;
}
