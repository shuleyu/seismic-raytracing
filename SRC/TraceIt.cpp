#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<string>
#include<iterator>
#include<algorithm>
#include<iomanip>
#include<set>
#include<map>
#include<complex>
#include<thread>
#include<atomic>
#include<unistd.h>

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
#include<ReadParameters.hpp>

#define _TURNINGANGLE 89.99

using namespace std;

enum PI{DebugInfo,TS,TD,RS,RD,StopAtSurface,nThread,FLAG1};
enum PS{InputRays,Layers,Depths,Ref,Polygons,ReceiverFileName,PolygonFilePrefix,RayFilePrefix,FLAG2};
enum PF{RectifyLimit,FLAG3};
ReadParameters<PI,PS,PF> P;

// Define the ray head structure.
class Ray {
    public:
        bool IsP,GoUp,GoLeft;
        string Color,Comp,Debug;
        int InRegion,Prev,RemainingLegs,Surfacing;
        double Pt,Pr,TravelTime,TravelDist,RayP,Amp,Inc,Takeoff;

        Ray()=default;
        Ray(bool p, bool g, bool l, string c, string cmp, 
            int i,int rl, double th, double r, double t, double d, double rp,double to) :
            IsP(p), GoUp(g), GoLeft(l), Color(c), Comp(cmp), Debug(""),
            InRegion(i), Prev(-1), RemainingLegs(rl), Surfacing(0),
            Pt(th),Pr(r),TravelTime(t),TravelDist(d), RayP(rp), Amp(1),Inc(0), Takeoff(to) {}
};


// Utilities for 1D-altering the PREM model.
vector<double> MakeRef(const double &depth,const vector<vector<double>> &dev){
    double rho=Drho(depth),vs=Dvs(depth),vp=Dvp(depth);
    for (const auto &item: dev) {
        if (item[0]<=depth && depth<=item[1]) {
            vp*=(1+item[2]/100);
            vs*=(1+item[3]/100);
            rho*=(1+item[4]/100);
            break;
        }
    }
    return {vp,vs,rho};
}

// Utilities for finding the index in an array that is cloeset to a given radius.
// array is sorted descending.
size_t findClosetLayer(const vector<double> &R, const double &r){
    auto rit=upper_bound(R.rbegin(),R.rend(),r);
    if (rit==R.rbegin()) return R.size()-1;
    else if (rit==R.rend()) return 0;
    else {
        auto it=rit.base();
        if (fabs(*it-r)<fabs(*rit-r)) return distance(R.begin(),it);
        else return distance(rit,R.rend())-1;
    }
}

// Utilities for finding the index in an array that is cloeset to a given depth.
// array is sorted ascending.
size_t findClosetDepth(const vector<double> &D, const double &d){
    auto it=upper_bound(D.begin(),D.end(),d);
    if (it==D.begin()) return 0;
    else if (it==D.end()) return D.size()-1;
    else {
        if (fabs(*it-d)<fabs(*prev(it)-d)) return distance(D.begin(),it);
        else return distance(D.begin(),it)-1;
    }
}

// generating rays born from RayHeads[i]
void followThisRay(size_t i,  atomic<size_t> &Cnt, atomic<int> &Estimation, atomic<size_t> &Running,
                   vector<string> &ReachSurfaces, vector<Ray> &RayHeads,
                   double RE, int branches, const vector<double> &specialDepths, const vector<vector<double>> &R,
                   const vector<vector<double>> &Vp,const vector<vector<double>> &Vs,const vector<vector<double>> &Rho,
                   double MinInc, const vector<vector<pair<double,double>>> &Regions, const vector<vector<double>> &RegionBounds,
                   const vector<double> &dVp, const vector<double> &dVs,const vector<double> &dRho){

    if (RayHeads[i].RemainingLegs==0) {
        Running.fetch_sub(1);
        return;
    }


    // Locate the begining and ending depths for the next leg.

    /// ... among special depths.

    //// Which special depth is cloest to ray head depth?
    double RayHeadDepth=RE-RayHeads[i].Pr;
    size_t Cloest=findClosetDepth(specialDepths,RayHeadDepth);

    //// Next depth should be the cloest special depth at the correct side (ray is going up/down).
    //// Is the ray going up or down? Is the ray already at the cloest special depth? If yes, adjust the next depth.
    double NextDepth=specialDepths[Cloest];
    if (RayHeads[i].GoUp && (specialDepths[Cloest]>RayHeadDepth || RayHeadDepth==specialDepths[Cloest]))
        NextDepth=specialDepths[Cloest-1];
    else if (!RayHeads[i].GoUp && (specialDepths[Cloest]<RayHeadDepth || specialDepths[Cloest]==RayHeadDepth))
        NextDepth=specialDepths[Cloest+1];

    double Top=min(RayHeadDepth,NextDepth),Bot=max(RayHeadDepth,NextDepth);

    /// ... among current 2D "Regions" vertical limits.
    int CurRegion=RayHeads[i].InRegion;
    Top=max(Top,RE-RegionBounds[CurRegion][3]);
    Bot=min(Bot,RE-RegionBounds[CurRegion][2]);

    // Print some debug info.
    if (P[DebugInfo]) {
        RayHeads[i].Debug+=to_string(1+i)+" --> ";
        cout << '\n' << "----------------------" ;
        cout << '\n' << "Calculating    : " << RayHeads[i].Debug;
        cout << "\nStart in region       : " << CurRegion;
        printf ("\nStart Location        : %.15lf deg, %.15lf km\n", RayHeads[i].Pt, RE-RayHeads[i].Pr);
        cout << "Will go as            : " << (RayHeads[i].IsP?"P, ":"S, ") << (RayHeads[i].GoUp?"Up, ":"Down, ")
             << (RayHeads[i].GoLeft?"Left":"Right") << endl;
        printf ("Ray tracing start, end: %.16lf --> %.16lf km with rayp: %.16lf sec/deg\n",Top,Bot,RayHeads[i].RayP);
        printf ("Search region bounds  : %.16lf ~ %.16lf\n", RE-R[CurRegion][0],RE-R[CurRegion].back());
        printf ("\nStart ray tracing ...\n\n");
        cout << flush;
    }


    // Use ray-tracing code "RayPath".
    size_t lastRadiusIndex;
    vector<double> degree;
    const auto &v=(RayHeads[i].IsP?Vp:Vs);
    auto ans=RayPath(R[CurRegion],v[CurRegion],RayHeads[i].RayP,Top,Bot,degree,lastRadiusIndex,_TURNINGANGLE);


    // If the new leg is trivia, no further operation needed.
    size_t RayLength=degree.size();
    if (RayLength==1) {
        RayHeads[i].RemainingLegs=0;
        Running.fetch_sub(1);
        return;
    }


    // This should never happen.
    // If the new leg is a reflection of down-going S to up-going P, and also the new leg turns, also mark it as invalid.
//     int PrevID=RayHeads[i].Prev;
//     if (PrevID!=-1 && !RayHeads[PrevID].GoUp && !RayHeads[PrevID].IsP && RayHeads[i].GoUp && RayHeads[i].IsP && ans.second) {
//         RayHeads[i].RemainingLegs=0;
//         Running.fetch_sub(1);
//         return;
//     }


    // Reverse the ray-tracing result if new leg is going upward.
    if (RayHeads[i].GoUp) {
        double totalDist=degree.back();
        for (auto &item:degree) item=totalDist-item;
        reverse(degree.begin(),degree.end());
    }


    // Create a projection from ray index to layer index.
    bool uP=RayHeads[i].GoUp;
    auto rIndex = [RayLength,lastRadiusIndex,uP](size_t j){
        if (uP) return (int)lastRadiusIndex-(int)j;
        else return (int)j+(int)lastRadiusIndex-(int)RayLength+1;
    };


    // Follow the new ray path to see if the new leg enters another region.
    int RayEnd=-1,NextRegion=-1,M=(RayHeads[i].GoLeft?-1:1);

    for (size_t j=0;j<degree.size();++j){

        pair<double,double> p={RayHeads[i].Pt+M*degree[j],R[CurRegion][rIndex(j)]}; // point on the newly calculated ray.

        if (CurRegion!=0){ // starts in some 2D polygon ...

            if (PointInPolygon(Regions[CurRegion],p,-1,RegionBounds[CurRegion])) continue; // ... and this point stays in that polygon.
            else { // ... but this point enters another polygon.

                RayEnd=j;

                // which region is the new leg entering?
                for (size_t k=1;k<Regions.size();++k){
                    if ((int)k==CurRegion) continue;
                    if (PointInPolygon(Regions[k],p,1,RegionBounds[k])) {
                        NextRegion=k;
                        break;
                    }
                }
                if (NextRegion==-1) NextRegion=0; // if can't find next 2D polygons, it must had return to the 1D reference region.
                break;
            }
        }
        else { // New leg starts in 1D reference region. Search for the region it enters.
            for (size_t k=1;k<Regions.size();++k){
                if (PointInPolygon(Regions[k],p,-1,RegionBounds[k])) { // If ray enters another region.
                    RayEnd=j;
                    NextRegion=k;
                    break;
                }
            }
            if (RayEnd!=-1) break; // If ray enters another region.
            else NextRegion=0;
        }
    }


    // Print some debug info.
    if (P[DebugInfo]) {
        cout << "Next refraction region is    : " << NextRegion << endl;
    }


    // Prepare reflection/refraction flags. (notice "rs" [r]eflection to [s]ame wave type is always possible)
    bool ts=(P[TS]==1),td=(P[TD]==1 && RayHeads[i].Comp!="SH"),rd=(P[RD]==1 && RayHeads[i].Comp!="SH");


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
        if (L1==Regions[SearchRegion].size())
            throw runtime_error("!!!!!!!!!!! Can't find junction! Bugs here !!!!!!!!!!!");


        // Find the junction point between ray and polygon boundary.
        const pair<double,double> &p1=Regions[SearchRegion][L1],&q1=Regions[SearchRegion][L2];

        double s1=(q1.second-p1.second)/(q1.first-p1.first),s2=(q2.second-p2.second)/(q2.first-p2.first);
        auto junc=LineJunction(p1,s1,p2,s2);
        JuncPt=junc.second.first;
        JuncPr=junc.second.second;


        // Print some debug info.
        if (P[DebugInfo]) printf("Junction at       : %.15lf,%.15lf.\n",JuncPt,JuncPr);


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
    if (P[DebugInfo]) {
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
    if (P[DebugInfo]) printf("rho1/vp1/vs1      : %.15lf g/cm3, %.15lf km/sec, %.15lf km/sec\n",rho1,vp1,vs1);
    if (P[DebugInfo]) printf("rho2/vp2/vs2      : %.15lf g/cm3, %.15lf km/sec, %.15lf km/sec\n",rho2,vp2,vs2);

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
    if (P[DebugInfo]) {
        printf("\nTakeOff angle (TS): %.15lf deg\n",Takeoff_ts);
        printf("RayP          (TS): %.15lf\n",Rayp_ts);
        printf("TakeOff angle (TD): %.15lf deg\n",Takeoff_td);
        printf("RayP          (TD): %.15lf deg\n",Rayp_td);
        printf("TakeOff angle (RD): %.15lf deg\n",Takeoff_rd);
        printf("RayP          (RD): %.15lf deg\n",Rayp_rd);
        printf("TakeOff angle (RS): %.15lf deg\n",Takeoff_rs);
        printf("RayP          (RS): %.15lf deg\n\n",Rayp_rs);

        cout << "RayEnd at Index  :" << RayEnd-1 << " (inclusive) / " << RayLength << endl;
        printf ("RayEnd at (for reflection)   :%.15lf deg, %.15lf (inclusive) km\n",NextPt_R,RE-NextPr_R);
        printf ("RayEnd at (for transmission) :%.15lf deg, %.15lf (inclusive) km\n\n",NextPt_T,RE-NextPr_T);
        cout << endl;
    }

    // Update current RayHead.
    RayHeads[i].TravelTime=ans.first.first;
    RayHeads[i].TravelDist=ans.first.second;
    --RayHeads[i].RemainingLegs;
    RayHeads[i].Inc=Incident;


    // Output valid part ray paths.
    if (P[RayFilePrefix]!="NONE") {
        ofstream fpout(P[RayFilePrefix]+to_string(i+1));
        fpout << "> " << (RayHeads[i].Color=="black"?(RayHeads[i].IsP?"blue":"red"):RayHeads[i].Color) << " "
                      << (RayHeads[i].IsP?"P ":"S ") << RayHeads[i].TravelTime << " sec. " << RayHeads[i].Inc << " IncDeg. "
                      << RayHeads[i].Amp << " DispAmp. " << RayHeads[i].TravelDist << " km. " << endl;
        for (int j=0;j<RayEnd;++j)
            fpout << RayHeads[i].Pt+M*degree[j] << " " << R[CurRegion][rIndex(j)] << '\n';
        fpout.close();
    }

    // If ray reaches surface, output info at the surface.
    if (fabs(NextPr_R-RE)<MinInc) ++RayHeads[i].Surfacing;
    if (fabs(NextPr_R-RE)<MinInc && (P[StopAtSurface]==0 || RayHeads[i].Surfacing<2)) {

        // Accumulate the travel-time.
        int I=i;
        double tt=0;
        vector<int> hh;
        while (I!=-1) {
            hh.push_back(I);
            tt+=RayHeads[I].TravelTime;
            I=RayHeads[I].Prev;
        }

        stringstream ss;
        ss << RayHeads[hh.back()].Takeoff << " " << RayHeads[i].RayP << " " << RayHeads[i].Inc << " " << NextPt_R << " "
           << tt << " " << RayHeads[i].Amp << " " << RayHeads[i].RemainingLegs << " ";
        for (auto rit=hh.rbegin();rit!=hh.rend();++rit)
            ss << (RayHeads[*rit].IsP?(RayHeads[*rit].GoUp?"p":"P"):(RayHeads[*rit].GoUp?"s":"S")) << ((*rit)==*hh.begin()?" ":"->");
        for (auto rit=hh.rbegin();rit!=hh.rend();++rit)
            ss << (1+*rit) << ((*rit)==*hh.begin()?"\n":"->");
        ReachSurfaces[i]=ss.str();
        if (P[StopAtSurface]==1) {
            int z=RayHeads[i].RemainingLegs;
            if (branches>1) z=(1-pow(branches,RayHeads[i].RemainingLegs))/(1-branches);
            Estimation.fetch_sub(branches*z);
            Running.fetch_sub(1);
            return;
        }
    }

    if (RayHeads[i].RemainingLegs==0) {
        Running.fetch_sub(1);
        return;
    }


    // Add rules of: (t)ransmission/refrection and (r)eflection to (s)ame or (d)ifferent way type.
    // Notice reflection with the same wave type is always allowed. ("rs" is always possible)

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
        Ray newRay=RayHeads[i];
        newRay.Prev=i;
        newRay.Pt=NextPt_T;
        newRay.Pr=NextPr_T;
        newRay.RayP=Rayp_ts;
        newRay.GoUp=(fabs(Takeoff_ts)>90);
        newRay.GoLeft=(Takeoff_ts<0);
        newRay.InRegion=NextRegion;
        double sign1=(T_PP.imag()==0?(T_PP.real()<0?-1:1):1);
        double sign2=(T_SS.imag()==0?(T_SS.real()<0?-1:1):1);
        newRay.Amp*=(newRay.IsP?(sign1*abs(T_PP)):(sign2*abs(T_SS)));
        RayHeads[Cnt.fetch_add(1)]=newRay;
    }

    if (td) {
        Ray newRay=RayHeads[i];
        newRay.IsP=!newRay.IsP;
        newRay.Prev=i;
        newRay.Pt=NextPt_T;
        newRay.Pr=NextPr_T;
        newRay.RayP=Rayp_td;
        newRay.GoUp=(fabs(Takeoff_td)>90);
        newRay.GoLeft=(Takeoff_td<0);
        newRay.InRegion=NextRegion;
        double sign1=(T_PS.imag()==0?(T_PS.real()<0?-1:1):1);
        double sign2=(T_SP.imag()==0?(T_SP.real()<0?-1:1):1);
        newRay.Amp*=(newRay.IsP?(sign1*abs(T_PS)):(sign2*abs(T_SP)));
        RayHeads[Cnt.fetch_add(1)]=newRay;
    }

    if (rd) {
        Ray newRay=RayHeads[i];
        newRay.IsP=!newRay.IsP;
        newRay.Prev=i;
        newRay.Pt=NextPt_R;
        newRay.Pr=NextPr_R;
        newRay.RayP=Rayp_rd;
        newRay.GoUp=(fabs(Takeoff_rd)>90);
        newRay.GoLeft=(Takeoff_rd<0);
        double sign1=(R_PS.imag()==0?(R_PS.real()<0?-1:1):1);
        double sign2=(R_SP.imag()==0?(R_SP.real()<0?-1:1):1);
        newRay.Amp*=(newRay.IsP?(sign1*abs(R_PS)):(sign2*abs(R_SP)));
        RayHeads[Cnt.fetch_add(1)]=newRay;
    }

    int y=(P[TS]!=0 && !ts)+(P[TD]!=0 && !td)+(P[RD]!=0 && !rd);
    int z=RayHeads[i].RemainingLegs;
    if (branches>1) z=(1-pow(branches,RayHeads[i].RemainingLegs))/(1-branches);
    Estimation.fetch_sub(y*z);

    // rs is always possible.
    if (P[RS]) {
        Ray newRay=RayHeads[i];
        newRay.Prev=i;
        newRay.Pt=NextPt_R;
        newRay.Pr=NextPr_R;
        newRay.RayP=Rayp_rs;
        newRay.GoUp=(fabs(Takeoff_rs)>90);
        newRay.GoLeft=(Takeoff_rs<0);
        double sign1=(R_PP.imag()==0?(R_PP.real()<0?-1:1):1);
        double sign2=(R_SS.imag()==0?(R_SS.real()<0?-1:1):1);
        newRay.Amp*=(newRay.IsP?(sign1*abs(R_PP)):(sign2*abs(R_SS)));
        RayHeads[Cnt.fetch_add(1)]=newRay;
    }

    Running.fetch_sub(1);
    return;
}


// dealt with inputs.
int main(int argc, char **argv){


    P=ReadParameters<PI,PS,PF> (argc,argv,cin,FLAG1,FLAG2,FLAG3);

    // Read in source settings.
    ifstream fpin;
    int steps,color;
    string comp;
    double theta,depth,takeoff;
    vector<int> initRaySteps,initRayComp,initRayColor;
    vector<double> initRayTheta,initRayDepth,initRayTakeoff;
    
    fpin.open(P[InputRays]);
    while (fpin >> theta >> depth >> takeoff >> comp >> color >> steps){
        // check.
        if (depth<0 || depth>6371)
            throw runtime_error("Source depth error @ line " + to_string(initRaySteps.size()+1) + "...");
        if (comp!="P" && comp!="SV" && comp!="SH")
            throw runtime_error("Source component error @ line " + to_string(initRaySteps.size()+1) + "...");
        if (steps<=0)
            throw runtime_error("Source step error: step<=0 @ line " + to_string(initRaySteps.size()+1) + "...");

        initRaySteps.push_back(steps);
        initRayTheta.push_back(Lon2360(theta));
        initRayDepth.push_back(depth);
        initRayTakeoff.push_back(Lon2180(takeoff));
        initRayComp.push_back(comp=="P"?0:(comp=="SV"?1:2));
        initRayColor.push_back(color);
    }
    fpin.close();
    // check.
    if (initRaySteps.empty()) throw runtime_error("No valid source, input format error?");


    // Read in grid setting.
    double depth1,depth2,inc,prev_depth=-numeric_limits<double>::max();
    vector<double> gridDepth1,gridDepth2,gridInc;
    fpin.open(P[Layers]);
    while (fpin >> depth1 >> depth2 >> inc){
        // check.
        if (prev_depth<0 && depth1!=0)
            throw runtime_error("Grid input error: not starting from depth=0 km ...");
        if (prev_depth>=0 && prev_depth!=depth1)
            throw runtime_error("Grid input error: not continuous ...");
        if (inc<=0)
            throw runtime_error("Grid input error: increment<=0 km ...");
        if (inc>depth2-depth1)
            throw runtime_error("Grid input error: increment .. Ha? ...");

        prev_depth=depth2;
        gridDepth1.push_back(depth1);
        gridDepth2.push_back(depth2);
        gridInc.push_back(inc);
    }
    fpin.close();
    // check.
    if (depth2!=6371) throw runtime_error("Grid input error: not ending at depth=6371 km ...");


    // Read in 1D reference deviation.
    double dvp,dvs,drho;
    vector<vector<double>> Deviation;
    map<double,int> refCheck;
    fpin.open(P[Ref]);
    while (fpin >> depth1 >> depth2 >> dvp >> dvs >> drho) {
        // check.
        if (depth1<0 || depth1>6371 || depth2<0 || depth2>6371 || depth2<=depth1)
            throw runtime_error("1D reference input depth error @ line "+ to_string(Deviation.size()+1) +" ...");
        if (dvp<-100 || dvs<-100 || drho<-100)
            throw runtime_error("1D reference input property error @ line "+ to_string(Deviation.size()+1) +" ...");

        ++refCheck[depth1];
        --refCheck[depth2];
        Deviation.push_back({depth1,depth2,dvp,dvs,drho});
    }
    fpin.close();
    //check.
    int refCnt=0;
    for (const auto &item: refCheck) {
        refCnt+=item.second;
        if (refCnt>1) throw runtime_error("1D reference input depth error: interval overlapping ...");
    }


    // Read in special depths (force reflection & refraction at these depths).
    set<double> specialSet{0,2891,5149.5,6371};
    fpin.open(P[Depths]);
    while (fpin >> depth) {
        //check.
        if (depth<0 || depth>6371)
            throw runtime_error("Additional reflection/transmission depth error @ line "+ to_string(specialSet.size()-3) +" ...");

        specialSet.insert(depth);
    }
    fpin.close();
    // sorted.
    vector<double> specialDepths(specialSet.begin(),specialSet.end());


    // Read in polygons (2D features).
    string tmpstr;
    char c;
    vector<vector<double>> regionProperties,regionPolygonsTheta,regionPolygonsDepth;
    fpin.open(P[Polygons]);
    while (getline(fpin,tmpstr)){
        if (tmpstr.empty()) continue;
        stringstream ss(tmpstr);

        if (tmpstr[0]=='>') { // a new region.
            ss >> c >> dvp >> dvs >> drho;
            if (dvp<-100 || dvs<-100 || drho<-100)
                throw runtime_error("2D region property error @ line "+ to_string(regionPolygonsTheta.size()+1) + " ...");
            regionProperties.push_back({dvp,dvs,drho});
            regionPolygonsTheta.push_back(vector<double> ());
            regionPolygonsDepth.push_back(vector<double> ());
        }
        else {
            ss >> theta >> depth;
            if (depth<0 || depth>6371)
                throw runtime_error("2D region depth error @ polygon "+ to_string(regionProperties.size()) +
                                    ", line: "+ to_string(regionPolygonsTheta.back().size()+1) + " ...");
            regionPolygonsTheta.back().push_back(Lon2360(theta));
            regionPolygonsDepth.back().push_back(depth);
        }
    }
    fpin.close();


    // Create ray tracing inputs.
    //
    // Currently we have these variables ------ :
    // ReadParameters<PI,PS,PF> P;
    // vector<int> initRaySteps,initRayComp,initRayColor;
    // vector<double> initRayTheta,initRayDepth,initRayTakeoff,gridDepth1,gridDepth2,gridInc,specialDepths;
    // vector<vector<double>> Deviation,regionProperties,regionPolygonsTheta,regionPolygonsDepth;
    //
    // For future I/O modification, you can stops here.

    const double RE=6371;

    // For plotting: 1D reference property deviation depths.
    if (!Deviation.empty() && P[PolygonFilePrefix]!="NONE"){
        ofstream fpout(P[PolygonFilePrefix]+"0");
        for (const auto &item:Deviation) {
            fpout << ">\n";
            for (double t=0;t<360;t=t+0.1) fpout << t << " " << RE-item[0] << '\n';
            fpout << ">\n";
            for (double t=0;t<360;t=t+0.1) fpout << t << " " << RE-item[1] << '\n';
        }
        fpout.close();
    }


    // Create 1D reference layers. (R[0]. 0 means 1D reference model)
    vector<vector<double>> R{vector<double> ()};
    for (size_t i=0;i<gridDepth1.size();++i){
        auto tmpr=CreateGrid(RE-gridDepth2[i],RE-gridDepth1[i],gridInc[i],2);
        if (!R[0].empty()) R[0].pop_back();
        R[0].insert(R[0].end(),tmpr.rbegin(),tmpr.rend());
    }


    // Fix round-off-errors: adding exact values in special depths and modefied 1D model to R[0].
    set<double> depthToCorrect(specialDepths.begin(),specialDepths.end());
    for (const auto &item:Deviation) {
        depthToCorrect.insert(item[0]);
        depthToCorrect.insert(item[1]);
    }
    vector<double> tmpArray;
    swap(R[0],tmpArray);
    tmpArray[0]=RE;tmpArray.back()=0;
    auto it=depthToCorrect.begin();
    for (int i=0;i<(int)tmpArray.size();++i) {
        if (it==depthToCorrect.end())
            R[0].push_back(tmpArray[i]);
        else if (tmpArray[i]==RE-*it) {
            R[0].push_back(tmpArray[i]);
            ++it;
        }
        else if (tmpArray[i]>RE-*it) {
            R[0].push_back(tmpArray[i]);
        }
        else {
            R[0].push_back(RE-*it);
            ++it;
            --i;
        }
    }
    while (it!=depthToCorrect.end()){
        R[0].push_back(RE-*it);
        ++it;
    }


    // ?? What's this for?
    // Find the minimum layer increment for later use. (To solve floating point number round-off-error-related issues)
    double MinInc=numeric_limits<double>::max();
    for (size_t i=1;i<R[0].size();++i) MinInc=min(MinInc,R[0][i-1]-R[0][i]);



    // Find the bounds of input polygons.
    vector<vector<double>> RegionBounds{{-numeric_limits<double>::max(),numeric_limits<double>::max(),
                                         -numeric_limits<double>::max(),numeric_limits<double>::max()}};
                                         // the 1D reference bounds is as large as possible.

    for (size_t i=0;i<regionPolygonsTheta.size();++i){
        double Xmin=numeric_limits<double>::max(),Xmax=-Xmin,Ymin=Xmin,Ymax=-Ymin;
        for (size_t j=0;j<regionPolygonsTheta[i].size();++j){
            size_t k=(j+1)%regionPolygonsTheta[i].size();
            double theta1=regionPolygonsTheta[i][j],theta2=regionPolygonsTheta[i][k];
            double radius1=RE-regionPolygonsDepth[i][j],radius2=RE-regionPolygonsDepth[i][k];

            Xmin=min(Xmin,theta1);Xmax=max(Xmax,theta2);
            Ymin=min(Ymin,radius1);Ymax=max(Ymax,radius2);
        }
        RegionBounds.push_back({Xmin,Xmax,Ymin,Ymax});
    }


    // Find the closest layer value in R[0] to Ymin/Ymax of each polygon.
    vector<size_t> adjustedYmin{0},adjustedYmax{0};
    for (size_t i=0;i<regionPolygonsTheta.size();++i){
        adjustedYmin.push_back(findClosetLayer(R[0],RegionBounds[i+1][2]));
        adjustedYmax.push_back(findClosetLayer(R[0],RegionBounds[i+1][3]));
    }


    // Rectify input polygons. And derived the polygon layers from the layers of the 1D reference:
    vector<pair<double,double>> tmpRegion;
    vector<vector<pair<double,double>>> Regions{tmpRegion}; // place holder for Region[0], which is the 1D reference.

    for (size_t i=0;i<regionPolygonsTheta.size();++i){

        tmpRegion.clear();
        for (size_t j=0;j<regionPolygonsTheta[i].size();++j){

            // Find the fine enough rectify for this section.
            size_t k=(j+1)%regionPolygonsTheta[i].size();
            double radius1=RE-regionPolygonsDepth[i][j],radius2=RE-regionPolygonsDepth[i][k];
            double theta1=regionPolygonsTheta[i][j],theta2=regionPolygonsTheta[i][k];

            if (radius1==RegionBounds[i+1][2]) radius1=R[0][adjustedYmin[i+1]];
            if (radius1==RegionBounds[i+1][3]) radius1=R[0][adjustedYmax[i+1]];
            if (radius2==RegionBounds[i+1][2]) radius2=R[0][adjustedYmin[i+1]];
            if (radius2==RegionBounds[i+1][3]) radius2=R[0][adjustedYmax[i+1]];

            double Tdist=theta2-theta1,Rdist=radius2-radius1;

            size_t NPTS=2;
            double dL=RE,dR,dT;
            while (dL>P[RectifyLimit]){
                NPTS*=2;
                dR=Rdist/(NPTS-1),dT=Tdist/(NPTS-1);
                dL=LocDist(theta1,0,radius1,theta1+dT,0,radius1+dR);
            }

            // Add rectified section to this polygon.
            for (size_t l=0;l+1<NPTS;++l)
                tmpRegion.push_back(make_pair(theta1+l*dT,radius1+l*dR));
        }

        // Add this rectified polygon to region array.
        Regions.push_back(tmpRegion);

    }

    // Output regions.
    if (P[PolygonFilePrefix]!="NONE"){
        for (size_t i=1;i<Regions.size();++i) {
            ofstream fpout(P[PolygonFilePrefix]+to_string(i));
            for (const auto &item:Regions[i]) fpout << item.first << " " << item.second << '\n';
            fpout.close();
        }
    }


    // adjust bounds to the values in R[0].
    for (size_t i=0;i<regionPolygonsTheta.size();++i){
        RegionBounds[i+1][2]=R[0][adjustedYmin[i+1]];
        RegionBounds[i+1][3]=R[0][adjustedYmax[i+1]];
    }

    // derive layers for this polygon.
    for (size_t i=0;i<regionPolygonsTheta.size();++i)
        R.push_back(vector<double> (R[0].begin()+adjustedYmax[i+1],R[0].begin()+adjustedYmin[i+1]+1));

    // properties for these polygon.
    vector<double> dVp{1},dVs{1},dRho{1}; // Region 0 has dVp=1 ,...
    for (size_t i=0;i<regionPolygonsTheta.size();++i) {
        dVp.push_back(1.0+regionProperties[i][0]/100);
        dVs.push_back(1.0+regionProperties[i][1]/100);
        dRho.push_back(1.0+regionProperties[i][2]/100);
    }

    // derive properties for all regions.
    vector<vector<double>> Vp(R.size(),vector<double> ()),Vs=Vp,Rho=Vp;
    for (size_t i=0;i<R.size();++i) {
        for (const auto &item:R[i]) {
            auto ans=MakeRef(RE-item,Deviation);
            Vp[i].push_back(dVp[i]*ans[0]);
            Vs[i].push_back(dVs[i]*ans[1]);
            Rho[i].push_back(dRho[i]*ans[2]);
        }
    }


    // Create initial rays.
    vector<Ray> RayHeads;
    size_t potentialSize=0;
    int branches=(P[TS]!=0)+(P[TD]!=0)+(P[RS]!=0)+(P[RD]!=0);

    for (size_t i=0;i<initRaySteps.size();++i){
        if (branches<=1) potentialSize+=initRaySteps[i];
        else potentialSize+=(1-pow(branches,initRaySteps[i]))/(1-branches);

        // Source in any polygons?
        int rid=0;
        for (size_t i=1;i<Regions.size();++i)
            if (PointInPolygon(Regions[i],make_pair(initRayTheta[i],RE-initRayDepth[i]),1,RegionBounds[i])) {rid=i;break;}

        // Calculate ray parameter.
        auto ans=MakeRef(initRayDepth[i],Deviation);
        double v=(initRayComp[i]==0?ans[0]*dVp[rid]:ans[1]*dVs[rid]);
        double rayp=M_PI/180*(RE-initRayDepth[i])*sin(fabs(initRayTakeoff[i])/180*M_PI)/v;

        // Push this ray into "RayHeads" for future processing.
        RayHeads.push_back(Ray(initRayComp[i]==0,fabs(initRayTakeoff[i])>=90,initRayTakeoff[i]<0,
                               (initRayColor[i]==0?"black":to_string(initRayColor[i])),
                               (initRayComp[i]==0?"P":(initRayComp[i]==1?"SV":"SH")),
                               rid,initRaySteps[i],initRayTheta[i],RE-initRayDepth[i],0,0,rayp,initRayTakeoff[i]));

    }

    if (potentialSize>RayHeads.max_size())
        throw runtime_error("Too many rays to handle: decrease the number of legs or the number of input rays...");

    atomic<size_t> Cnt;
    Cnt.store(RayHeads.size());
    atomic<int> Estimation;
    Estimation.store(potentialSize);
    RayHeads.resize(potentialSize);

    // Start ray tracing. (Finally!)
    //
    // Process each "Ray" leg in "RayHeads".
    // For future legs generated by reflction/refraction, create new "Ray" and assign it to the proper position in "RayHeads" vector.
    vector<thread> allThread(RayHeads.capacity());
    vector<string> ReachSurfaces(RayHeads.capacity());

    size_t Doing=0,Done=0;
    atomic<size_t> Running;
    Running.store(0);
    while (Doing!=Cnt.load() || Running.load()!=0) {

        if (Running.load()<(size_t)P[nThread] && Doing<Cnt.load()) {
            Running.fetch_add(1);
            allThread[Doing]=thread(followThisRay,Doing, std::ref(Cnt), std::ref(Estimation), std::ref(Running), std::ref(ReachSurfaces),
                             std::ref(RayHeads), RE, branches, std::cref(specialDepths),
                             std::cref(R),std::cref(Vp),std::cref(Vs),std::cref(Rho),MinInc,
                             std::cref(Regions),std::cref(RegionBounds),std::cref(dVp),std::cref(dVs),std::cref(dRho));
            if (Doing>0 && Doing%10000==0) {
                for (size_t i=Done;i<Doing-P[nThread];++i)
                    allThread[i].join();
                Done=Doing-P[nThread];
            }
            ++Doing;
        }
        else usleep(1000);

    } // End of ray tracing.

    for (size_t i=Done;i<Doing;++i)
        allThread[i].join();

    ofstream fpout(P[ReceiverFileName]);
    fpout << "<Takeoff> <Rayp> <Incident> <Dist> <TravelTime> <DispAmp> <RemainingLegs> <WaveTypeTrain> <RayTrain>" << '\n';
    for (const string &s: ReachSurfaces)
        fpout << s;
    fpout.close();

cout << Cnt.load() << "/" << RayHeads.capacity() << "/" << Estimation.load() << endl;

    return 0;
}
