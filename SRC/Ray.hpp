#ifndef RAY
#define RAY

#include<stdlib.h>
#include<string.h>
#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<vector>
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

#include<Lon2180.hpp>
#include<Lon2360.hpp>

#define _TURNINGANGLE 89.99
#define _RE 6371

// Define the ray node.
class Ray {
    public:
        bool IsP,GoUp,GoLeft;
        std::string Comp,Debug;
        int InRegion,Prev,RemainingLegs,Surfacing,Color;
        double Pt,Pr,TravelTime,TravelDist,RayP,Amp,Inc,Takeoff;

        Ray()=default;
        Ray(bool p, bool g, bool l, std::string cmp,
            int i,int rl, int c, double th, double r, double t, double d, double rp,double to) :
            IsP(p), GoUp(g), GoLeft(l), Comp(cmp), Debug(""),
            InRegion(i), Prev(-1), RemainingLegs(rl), Surfacing(0),Color(c),
            Pt(th),Pr(r),TravelTime(t),TravelDist(d), RayP(rp), Amp(1),Inc(0), Takeoff(to) {}
};

// Declarations.
std::vector<double> MakeRef(const double &depth,const std::vector<std::vector<double>> &dev);
std::size_t findClosetLayer(const std::vector<double> &R, const double &r);
std::size_t findClosetDepth(const std::vector<double> &D, const double &d);
void followThisRay(
    std::size_t i, std::atomic<std::size_t> &Cnt, std::atomic<int> &Estimation, std::atomic<std::size_t> &Running,
    char **ReachSurfaces, int *ReachSurfacesSize, char **RayInfo, int *RayInfoSize,
    double **RaysTheta, int *RaysN, double **RaysRadius,
    std::vector<Ray> &RayHeads, int branches, const std::vector<double> &specialDepths,
    const std::vector<std::vector<double>> &R, const std::vector<std::vector<double>> &Vp,
    const std::vector<std::vector<double>> &Vs,const std::vector<std::vector<double>> &Rho,
    const std::vector<std::vector<std::pair<double,double>>> &Regions, const std::vector<std::vector<double>> &RegionBounds,
    const std::vector<double> &dVp, const std::vector<double> &dVs,const std::vector<double> &dRho,
    const bool &DebugInfo,const bool &TS,const bool &TD,const bool &RS,const bool &RD, const bool &StopAtSurface);
void PreprocessAndRun(
    const std::vector<int> &initRaySteps,const std::vector<int> &initRayComp,const std::vector<int> &initRayColor,
    const std::vector<double> &initRayTheta,const std::vector<double> &initRayDepth,const std::vector<double> &initRayTakeoff,
    const std::vector<double> &gridDepth1,const std::vector<double> &gridDepth2,const std::vector<double> &gridInc,
    const std::vector<double> &specialDepths,const std::vector<std::vector<double>> &Deviation,
    const std::vector<std::vector<double>> &regionProperties,
    const std::vector<std::vector<double>> &regionPolygonsTheta,
    const std::vector<std::vector<double>> &regionPolygonsDepth,
    const double &RectifyLimit, const bool &TS, const bool &TD, const bool &RS, const bool &RD,
    const std::size_t &nThread, const bool &DebugInfo, const bool &StopAtSurface,
    const std::size_t &branches, const std::size_t &potentialSize,
    char **ReachSurfaces, int *ReachSurfacesSize, char **RayInfo, int *RayInfoSize,
    int *RegionN,double **RegionsTheta,double **RegionsRadius,
    double **RaysTheta, int *RaysN, double **RaysRadius,int *Observer);

#endif
