#include<vector>
#include<set>
#include<cmath>
#include<algorithm>
#include<complex>
#include<atomic>
#include<thread>
#include<unistd.h>
#include<string.h>
#ifndef RAY
#define RAY



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
#ifndef ASU_LON2180
#define ASU_LON2180

/*******************************************************
 * This C++ template convert input longitude (in deg) to
 * [-180,180] deg.
 *
 * input(s):
 * const double &lon  ----  Longitude.
 *
 * return(s):
 * double ans  ----  Longitude in -180 ~ 180.
 *
 * Shule Yu
 * Dec 29 2017
 *
 * Key words: coordinates, longitude.
*******************************************************/

inline double Lon2180(const double &lon){

    double ans=lon;
    if (ans>=0) ans-=360*((int)(ans/360));
    else ans+=360*(1+(int)(-ans/360));

    if (ans>180) ans-=360;
    return ans;
}

#endif
#ifndef ASU_LON2360
#define ASU_LON2360

/***********************************************************
 * This C++ template convert input longitude (in deg) to
 * [0,360] deg.
 *
 * input(s):
 * const double &lon  ----  Longitude.
 *
 * return(s):
 * double ans  ----  Longitude in 0 ~ 360.
 *
 * Shule Yu
 * Dec 29 2017
 *
 * Key words: coordinates, longitude.
***********************************************************/

inline double Lon2360(const double &lon){

    double ans=lon;
    if (ans>=0) ans-=360.0*((int)(ans/360));
    else ans+=360.0*(1+(int)(-ans/360));

    if (ans>=360) ans=0;
    return ans;
}

#endif
#ifndef ASU_LL2TP
#define ASU_LL2TP



/*****************************************************
 * This C++ template convert Longitude, Latitude ( in
 * deg. ) to Theta, Phi ( in rad. ).
 *
 * The Cartesian coordinate x,y,z axis is as follows:
 * x axis: lon=0 ,lat=0
 * y axis: lon=90,lat=0.
 * z axis: lat=90.
 *
 * Theta : angle with z, 0 <= theta <= PI.
 * Phi   : angle with x, counter-clockwise.
 *
 * Longitude: -180 < lon <= 180. (Lon=0 at poles.)
 * Latitude : -90 <= lat <= 90.
 *
 * input(s):
 * const double &Lon  ----  Longitude.
 * const double &Lat  ----  Latitude.
 *
 * return(s):
 * pair<double,double> ans  ----  {theta,phi}
 *
 * Shule Yu
 * Dec 29 2017
 *
 * Key words: geography, coordinates, convertion.
*****************************************************/

inline std::pair<double,double> LL2TP(const double &lon, const double &lat){

    // Deal with latitude, theta.
    // Check latitude.
    if ( lat < -90 || lat > 90 ){
        return {};
    }

    if (lat==90) return {0.0,0.0};
    else if (lat==-90) return {M_PI,0.0};

    return {M_PI/2-lat*M_PI/180,Lon2180(lon)*M_PI/180};
}

#endif
#ifndef ASU_SPH2CART
#define ASU_SPH2CART


/***********************************************************
 * This C++ template convert r, theta, phi (in rad) in
 * Spherical Cartesian coordinates to x, y, z in
 * Cartesian coordinates.
 *
 * input(s):
 * const double &r      ----  r
 * const double &theta  ----  theta, angle with positive z direction.  (rad)
 * const double &phi    ----  phi, angle couter-clock wise with positive x direction. (rad)
 *
 * return(s):
 * vector<double> ans  ----  {x,y,z} in Cartesian.
 *
 * Shule Yu
 * Dec 29 2017
 *
 * Key words: coordinates, convertion, spherical, cartesian
 *
 * Reference: https://en.wikipedia.org/wiki/Spherical_coordinate_system
***********************************************************/

inline std::vector<double> Sph2Cart(const double &r, const double &theta, const double &phi){
    return {r*cos(phi)*sin(theta),r*sin(phi)*sin(theta),r*cos(theta)};
}

#endif
#ifndef ASU_DOTDIST
#define ASU_DOTDIST


/*********************************************************
 * This C++ template return the straight line distance
 * between two Cartesian input points.
 *
 * input(s):
 * const double x1  ----  x position of point 1.
 * const double y1  ----  y position of point 1.
 * const double z1  ----  z position of point 1.
 * const double x2  ----  x position of point 2.
 * const double y2  ----  y position of point 2.
 * const double z2  ----  z position of point 2.
 *
 * return(s):
 * double ans  ----  distance between point 1 and point 2.
 *
 * Shule Yu
 * Mar 30 2017
 *
 * Key words: distance.
*********************************************************/

inline double DotDist(const double &x1, const double &y1,const double &z1,
                      const double &x2, const double &y2,const double &z2){
    return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
}

#endif
#ifndef ASU_CROSSPRODUCT
#define ASU_CROSSPRODUCT


/*********************************************************
 * This C++ template take the cross product of two vectors
 * (x1,y1,z1) and (x2,y2,z2).
 *
 * input(s):
 * const double &x1  ----  x value of vector 1.
 * const double &y1  ----  y value of vector 1.
 * const double &z1  ----  z value of vector 1.
 * const double &x2  ----  x value of vector 2.
 * const double &y2  ----  y value of vector 2.
 * const double &z2  ----  z value of vector 2.
 *
 * return(s):
 * vector<double> ans  ----  {x3,y3,z3}
 *
 * Shule Yu
 * Dec 28 2017
 *
 * Key words: cross product.
*********************************************************/

inline std::vector<double> CrossProduct(const double &x1, const double &y1, const double &z1,
                                        const double &x2, const double &y2, const double &z2){
    return {(1.0*y1*z2-y2*z1),(1.0*z1*x2-z2*x1),(1.0*x1*y2-x2*y1)};
}

#endif
#ifndef ASU_CREATEGRID
#define ASU_CREATEGRID


/***********************************************************
 * This C++ template returns the 1D grid meshing results.
 *
 * input(s):
 * const double &lowerbound   ----  Grid lower bound.
 * const double &upperbound   ----  Grid upper bound.
 * const double &Para         ----  meaning-changing parameter.
 *                                  it's meaning depends on which mode:
 *
 * const int    &mode         ----  select mode.
 *                                  Define the meaning of Para and return vector.
 *
 *                                  0: Para means number of gird points.
 *                                     Return a vector (ans) of size Para:
 *
 *                                     ans.size()=Para;
 *                                     ans[0]=lowerbound;
 *                                     ans[Para-1]=upperbound;
 *
 *                                  1: Para means grid increment.
 *                                     Last grid point is less equal to higherbound.
 *                                     Return a vector (ans) of calculated size:
 *
 *                                     ans[0]=lowerbound;
 *                                     ans[1]-ans[0]=Para;
 *                                     upperbound-Para<ans.back();
 *                                     ans.back()<=upperbound;
 *
 *                                 -1: Same as 1. Will only return the grid property:
 *
 *                                     ans.size()=2;
 *                                     ans[0] = calculated grid size.
 *                                     ans[1] = adjusted upperbound.
 *
 *
 *                                  2: Para means an estimation of grid increment.
 *                                     The calculated grid increment (Para*) is (possibly)
 *                                     sightly adjusted such that the higherbound is meet.
 *                                     Return a vector (ans) of calculated size:
 *
 *                                     ans[0]=lowerbound;
 *                                     ans[1]-ans[0]=Para* (<=Para);
 *                                     ans.back()=upperbound;
 *
 *                                 -2: Same as 2. Will only return the grid property:
 *
 *                                     ans.size()=2;
 *                                     ans[0] = calculated grid size.
 *                                     ans[1] = adjusted Para (Para*).
 *
 * return(s):
 * vector<double> ans  ----  Created grid or grid properties, depending on mode.
 *
 * Shule Yu
 * Jan 23 2018
 *
 * Key words: creating grid.
***********************************************************/

std::vector<double> CreateGrid(const double &lowerbound, const double &upperbound,
                               const double &Para, const int mode){

    // check lower and upper bound.
    if (upperbound<lowerbound) {
        return {};
    }

    if (mode==0){

        int N=Para;
        double Inc=1.0*(upperbound-lowerbound)/(N-1);

        std::vector<double> ans(N,lowerbound);
        for (int i=1;i<N;++i) ans[i]=ans[i-1]+Inc;
        return ans;
    }
    if (mode==1 || mode==-1){

        double Inc=Para;

        if (mode==-1) {
            double N=1+floor(1.0*(upperbound-lowerbound)/Inc);
            if (fabs((upperbound-lowerbound-N*Inc))<Inc*1e-6) // float number rounding error?
                N+=1;
            return {N,lowerbound+(N-1)*Inc};
        }

        std::vector<double> ans;

        double Cur=lowerbound;

        while (Cur<=upperbound || Cur-upperbound<Inc*1e-6) {
            ans.push_back(Cur);
            Cur+=Inc;
        }
        return ans;
    }
    if (mode==2 || mode==-2){

        double Inc=Para;
        int N=1+(int)round((upperbound-lowerbound)/Inc);
        Inc=1.0*(upperbound-lowerbound)/(N-1);

        if (mode==-2) {
            return {1.0*N,Inc};
        }

        std::vector<double> ans;
        double Cur=lowerbound;

        for (int i=0;i<N;++i) {
            ans.push_back(Cur);
            Cur+=Inc;
        }

        // Round-off errors eliminator!
        ans.back()=upperbound;
        return ans;
    }

    return {};
}

#endif
#ifndef ASU_LINEJUNCTION
#define ASU_LINEJUNCTION


/**************************************************************
 * This C++ template return the junction point of two straigh
 * lines.
 *
 * input(s):
 * const pair<T1,T2> &p1  ----  first point.
 * const double      &s1  ----  slope through the first point.
 * const pair<T3,T4> &p2  ----  second point.
 * const double      &s2  ----  slope through the second point.
 *
 * return(s):
 * std::pair<bool,pair<double,double>>  ans
 *                       ----  {Result flag, the position of the junction point}
 *                             If lines are parallel(no junction), flag=false.
 *
 * Shule Yu
 * Jan 19 2018
 *
 * Key words: line junction, intercept, cross.
****************************************************************/

template<typename T1,typename T2,typename T3,typename T4>
std::pair<bool,std::pair<double,double>> LineJunction(const std::pair<T1,T2> &p1,const double &s1,
                                                      const std::pair<T3,T4> &p2,const double &s2){

    std::pair<double,double> ans{std::numeric_limits<double>::max(),std::numeric_limits<double>::max()};

    if (s1==s2 || (std::isinf(s1) && std::isinf(s2))) return {false,ans};

    double b1=p1.second-s1*p1.first,b2=p2.second-s2*p2.first;
    if (std::isinf(s1)){
        ans.first=p1.first;
        ans.second=s2*p1.first+b2;
    }
    else if (std::isinf(s2)){
        ans.first=p2.first;
        ans.second=s1*p2.first+b1;
    }
    else{
        ans.second=((s2*b1-s1*b2)/(s2-s1));
        if (s1!=0) ans.first=(ans.second-b1)/s1;
        else ans.first=(ans.second-b2)/s2;
    }

    return {true,ans};
}

#endif
#ifndef ASU_LOCDIST
#define ASU_LOCDIST


/***********************************************************
 * This C++ template caculate the straight line distance
 * between two geographic points.
 *
 * input(s):
 * const double &lon1   ----  p1's longitude.
 * const double &lat1   ----  p1's latitude.
 * const double &r1     ----  p1's distance to Earth's center.
 * const double &lon2   ----  p2's longitude.
 * const double &lat2   ----  p2's latitude.
 * const double &r2     ----  p2's distance to Earth's center.
 *
 * return(s):
 * double ans  ----  Straight line distance between these two points.
 *
 * Shule Yu
 * Dec 31 2017
 *
 * Key words: geography, coordinates, distance.
***********************************************************/

double LocDist(const double &lon1, const double &lat1, const double &r1,
               const double &lon2, const double &lat2, const double &r2){

    auto res1=LL2TP(lon1,lat1);
    auto res2=LL2TP(lon2,lat2);

    auto p1=Sph2Cart(r1,res1.first,res1.second);
    auto p2=Sph2Cart(r2,res2.first,res2.second);

    return DotDist(p1[0],p1[1],p1[2],p2[0],p2[1],p2[2]);
}

#endif
#ifndef ASU_PREM
#define ASU_PREM


// PREM            // original PREM
//
// PREMSmoothed    // Replace first-order discontinuities.
//                     with straight line over +- 30km.
//
// PREMX           // Spline interpolate from surface to 400km.
//                 // Then degree2 from 600 km to the proper depth.

/***********************************************************
 * This C++ template returns PREM model values at certain depth.
 * There's some small template you can use for just one value.
 *
 * input(s):
 * const double &Depth  ----  Given depth (in km)
 * const int    &iso    ----  0/1: (set to 1 if want isotropy result)
 * const int    &ocean  ----  0/1: (set to 1 if want ocean properties near surface)
 * double       &rho    ----  density (g/cm^3)
 * double       &vpv    ----  Vertical polarity P-wave speed (P-wave speed when traveling along the radius direction)
 * double       &vph    ----  Horizontal polarity P-wave speed (P-wave speed when traveling perpendicular to the radius direction)
 * double       &vsv    ----  Vertical polarity S-wave speed (SV-wave speed when traveling perpendicular to the radius direction)
 * double       &vsh    ----  Horizontal polarity S-wave speed (SH-wave speed when traveling perpendicular to the radius direction)
 *                                                             (also, SV- and SH-wave speed when traveling along the radius direction)
 * double       &qu     ----  Isotropic dissipation of shear energy.
 * double       &qk     ----  Isotropic dissipation of compressional energy.
 * double       &yita   ----  Anisotropic parameter.
 *
 * output(s):
 * double &rho  (in-place)
 * double &vpv  (in-place)
 * double &vph  (in-place)
 * double &vsv  (in-place)
 * double &vsh  (in-place)
 * double &qu   (in-place)
 * double &qk   (in-place)
 * double &yita (in-place)
 *
 * Shule Yu
 * Dec 29 2017
 *
 * Key words: PREM, prem
 *
 * Reference :
 *     Physics_of_the_earth_and_planetary_interiors_1981_Dziewonski
***********************************************************/

void PREM(const double &Depth,const int &iso, const int &ocean,
          double &rho, double &vpv, double &vph, double &vsv,
          double &vsh, double &qu, double &qk, double &yita,
          const bool &PrintWarning=true){

    double R=6371.0;
    double r=R-Depth;
    double x=r/R;

    yita=sqrt(-1);
    qk=57823;

    if ( 0<=r && r<1221.5 ) {
    // inner core

        rho= -8.8381*x*x + 13.0885;
        vpv= -6.3640*x*x + 11.2622;
        vph=  vpv;
        vsv= -4.4475*x*x + 3.6678;
        vsh=  vsv;
        qu =  84.6;
        qk =  1327.7;
    }

    else if ( 1221.5<=r && r<3480.0 ) {
    // 2891 ~ 5149.5 km, outter core.

        rho=  -5.5281*x*x*x - 3.6426*x*x - 1.2638*x + 12.5815;
        vpv= -13.5732*x*x*x + 4.8023*x*x - 4.0362*x + 11.0487;
        vph=  vpv;
        vsv=  0;
        vsh=  0;
        qu =  INFINITY;
    }

    else if ( 3480.0<=r && r<3630.0 ) {
    // 2741 ~ 2891 km, D''

        rho= -3.0807*x*x*x + 5.5283*x*x - 6.4761*x + 7.9565;
        vpv= -2.5514*x*x*x + 5.5242*x*x - 5.3181*x + 15.3891;
        vph=  vpv;
        vsv=  0.9783*x*x*x - 2.0834*x*x + 1.4672*x + 6.9254;
        vsh=  vsv;
        qu =  312;
    }

    else if ( 3630.0<=r && r<5600.0 ) {
    // 771 ~ 2741 km

        rho=  -3.0807*x*x*x +  5.5283*x*x -  6.4761*x +  7.9565;
        vpv= -26.6419*x*x*x + 51.4832*x*x - 40.4673*x + 24.9520;
        vph=  vpv;
        vsv=  -9.2777*x*x*x + 17.4575*x*x - 13.7818*x + 11.1671;
        vsh=  vsv;
        qu =  312;
    }

    else if ( 5600.0<=r && r<5701.0 ) {
    // 670 ~ 771 km

        rho= -3.0807*x*x*x + 5.5283*x*x -  6.4761*x +  7.9565;
        vpv= -2.5514*x*x*x + 5.5242*x*x - 23.6027*x + 29.2766;
        vph=  vpv;
        vsv=  0.9783*x*x*x - 2.0834*x*x - 17.2473*x + 22.3459;
        vsh=  vsv;
        qu =  312;
    }

    else if ( 5701.0<=r && r<5771.0 ) {
    // 600 ~ 670 km

        rho= -1.4836*x +  5.3197;
        vpv= -9.8672*x + 19.0957;
        vph=  vpv;
        vsv= -4.9324*x +  9.9839;
        vsh=  vsv;
        qu =  143;
    }

    else if ( 5771.0<=r && r<5971.0 ) {
    // 400 ~ 600 km

        rho=  -8.0298*x + 11.2494;
        vpv= -32.6166*x + 39.7027;
        vph=  vpv;
        vsv= -18.5856*x + 22.3512;
        vsh=  vsv;
        qu =  143;
    }

    else if ( 5971.0<=r && r<6151.0 ) {
    // 220 ~ 400 km

        rho=  -3.8045*x +  7.1089;
        vpv= -12.2569*x + 20.3926;
        vph=  vpv;
        vsv=  -4.4597*x +  8.9496;
        vsh=  vsv;
        qu =  143;
    }

    else if ( 6151.0<=r && r<6291.0 ) {
    // 80 ~ 220 km, LVZ

        rho=  0.6924*x + 2.6910;
        vpv=  7.2180*x + 0.8317;
        vph=  4.6172*x + 3.5908;
        vsv= -1.4678*x + 5.8582;
        vsh=  5.7176*x - 1.0839;
        qu =  80;
        yita= -2.4778*x + 3.3687;
    }

    else if ( 6291.0<=r && r<6346.6 ) {
    // 24.4 ~ 80 km, LID

        rho=  0.6924*x + 2.6910;
        vpv=  7.2180*x + 0.8317;
        vph=  4.6172*x + 3.5908;
        vsv= -1.4678*x + 5.8582;
        vsh=  5.7176*x - 1.0839;
        qu =  600;
        yita= -2.4778*x + 3.3687;
    }

    else if ( 6346.6<=r && r<6356.0 ) {
    // crust 1
        rho=  2.900;
        vpv=  6.800;
        vph=  vpv;
        vsv=  3.900;
        vsh=  vsv;
        qu =  600;
    }

    else if ( 6356.0<=r && r<=R+1e-5 ) {
    // crust 2

        rho=  2.600;
        vpv=  5.800;
        vph=  vpv;
        vsv=  3.200;
        vsh=  vsv;
        qu =  600;
    }

    else { // space

        if (PrintWarning)

        rho=0;
        vpv=0;
        vph=0;
        vsv=0;
        vsh=0;
        qu =0;
        qk =0;
        yita=0;

    }

    if ( ocean==1 && 6368.0<=r && r<=R+1e-5 ) {
    // ocean

        rho=  1.020;
        vpv=  1.450;
        vph=  vpv;
        vsv=  0;
        vsh=  vsv;
        qu =  INFINITY;
    }

    if ( iso==1 && 6151.0<=r && r<6291.0 ) {
    // LVZ iso

        rho=  0.6924*x + 2.6910;
        vpv=  3.9382*x + 4.1875;
        vph=  vpv;
        vsv=  2.3481*x + 2.1519;
        vsh=  vsv;
        qu =  80;
        yita= -2.4778*x + 3.3687;
    }

    if ( iso==1 && 6291.0<=r && r<6346.6 ) {
    // LID iso

        rho=  0.6924*x + 2.6910;
        vpv=  3.9382*x + 4.1875;
        vph=  vpv;
        vsv=  2.3481*x + 2.1519;
        vsh=  vsv;
        qu =  600;
        yita= -2.4778*x + 3.3687;
    }

    return ;
}

double Dvs(const double &Depth, const bool &PrintWarning=true) {
    double rho=0,vpv=0,vph=0,vsv=0,vsh=0,qu=0,qk=0,yita=0;
    PREM(Depth,1,0,rho,vpv,vph,vsv,vsh,qu,qk,yita,PrintWarning);
    return vsv;
}

double Dvp(const double &Depth, const bool &PrintWarning=true) {
    double rho=0,vpv=0,vph=0,vsv=0,vsh=0,qu=0,qk=0,yita=0;
    PREM(Depth,1,0,rho,vpv,vph,vsv,vsh,qu,qk,yita,PrintWarning);
    return vpv;
}

double Drho(const double &Depth, const bool &PrintWarning=true) {
    double rho=0,vpv=0,vph=0,vsv=0,vsh=0,qu=0,qk=0,yita=0;
    PREM(Depth,1,0,rho,vpv,vph,vsv,vsh,qu,qk,yita,PrintWarning);
    return rho;
}

double Rvs(const double &Radius, const bool &PrintWarning=true) {
    return Dvs(6371.0-Radius,PrintWarning);
}

double Rvp(const double &Radius, const bool &PrintWarning=true) {
    return Dvp(6371.0-Radius,PrintWarning);
}

double Rrho(const double &Radius, const bool &PrintWarning=true) {
    return Drho(6371.0-Radius,PrintWarning);
}

// PREMSmoothed

void PREMSmoothed(const double &Depth, double &rho,double &vpv,double &vph,
                  double &vsv,double &vsh,double &qu,double &qk,double &yita,
                  const int &RemoveCrust,const int &Remove220,
                  const int &Remove400, const int &Remove670,
                  const bool &PrintWarning=true){

// Inspired from Mike Thorne.
// By default, no ocean (0), no anisotropy (1).
// Linear interpolate between +- 30km around each discontinuity.

    double r=6371.0-Depth;
    PREM(Depth,1,0,rho,vpv,vph,vsv,vsh,qu,qk,yita,PrintWarning);

    // Remove Crust.
    if ( RemoveCrust==1 && 6346.6<=r ) {
        vsv=4.49100712;
        vpv=8.11061727;
        rho=3.38074821;
    }

    // Remove 220km.
    if ( Remove220==1 && 6121.0<=r && r<6181.0 ){
        vsv=Rvs(6181.0)+(Rvs(6121.0)-Rvs(6181.0))/(6121.0-6181.0)*(r-6181.0);
        vpv=Rvp(6181.0)+(Rvp(6121.0)-Rvp(6181.0))/(6121.0-6181.0)*(r-6181.0);
        rho=Rrho(6181.0)+(Rrho(6121.0)-Rrho(6181.0))/(6121.0-6181.0)*(r-6181.0);
    }

    // Remove 400km.
    if ( Remove400==1 && 5941.0<=r && r<6001.0 ){
        vsv=Rvs(6001.0)+(Rvs(5941.0)-Rvs(6001.0))/(5941.0-6001.0)*(r-6001.0);
        vpv=Rvp(6001.0)+(Rvp(5941.0)-Rvp(6001.0))/(5941.0-6001.0)*(r-6001.0);
        rho=Rrho(6001.0)+(Rrho(5941.0)-Rrho(6001.0))/(5941.0-6001.0)*(r-6001.0);
    }

    // Remove 670km.
    if ( Remove670==1 && 5671.0<=r && r<5731.0 ){
        vsv=Rvs(5731.0)+(Rvs(5671.0)-Rvs(5731.0))/(5671.0-5731.0)*(r-5731.0);
        vpv=Rvp(5731.0)+(Rvp(5671.0)-Rvp(5731.0))/(5671.0-5731.0)*(r-5731.0);
        rho=Rrho(5731.0)+(Rrho(5671.0)-Rrho(5731.0))/(5671.0-5731.0)*(r-5731.0);
    }

    vsh=vsv;
    vph=vpv;

    return;
}

double DvsS(const double &Depth,const int &RemoveCrust, const int &Remove220,
            const int &Remove400, const int &Remove670, const bool &PrintWarning=true) {
    double rho=0,vpv=0,vph=0,vsv=0,vsh=0,qu=0,qk=0,yita=0;
    PREMSmoothed(Depth,rho,vpv,vph,vsv,vsh,qu,qk,yita,RemoveCrust,Remove220,Remove400,Remove670,PrintWarning);
    return vsv;
}

double DvpS(const double &Depth,const int &RemoveCrust, const int &Remove220,
            const int &Remove400, const int &Remove670, const bool &PrintWarning=true) {
    double rho=0,vpv=0,vph=0,vsv=0,vsh=0,qu=0,qk=0,yita=0;
    PREMSmoothed(Depth,rho,vpv,vph,vsv,vsh,qu,qk,yita,RemoveCrust,Remove220,Remove400,Remove670,PrintWarning);
    return vpv;
}

double DrhoS(const double &Depth,const int &RemoveCrust, const int &Remove220,
            const int &Remove400, const int &Remove670, const bool &PrintWarning=true) {
    double rho=0,vpv=0,vph=0,vsv=0,vsh=0,qu=0,qk=0,yita=0;
    PREMSmoothed(Depth,rho,vpv,vph,vsv,vsh,qu,qk,yita,RemoveCrust,Remove220,Remove400,Remove670,PrintWarning);
    return rho;
}

double RvsS(const double &Radius,const int &RemoveCrust, const int &Remove220,
            const int &Remove400, const int &Remove670, const bool &PrintWarning=true) {
    return DvsS(6371.0-Radius,RemoveCrust,Remove220,Remove400,Remove670,PrintWarning);
}

double RvpS(const double &Radius,const int &RemoveCrust, const int &Remove220,
            const int &Remove400, const int &Remove670, const bool &PrintWarning=true) {
    return DvpS(6371.0-Radius,RemoveCrust,Remove220,Remove400,Remove670,PrintWarning);
}

double RrhoS(const double &Radius,const int &RemoveCrust, const int &Remove220,
            const int &Remove400, const int &Remove670, const bool &PrintWarning=true) {
    return DrhoS(6371.0-Radius,RemoveCrust,Remove220,Remove400,Remove670,PrintWarning);
}

// PREMX

void PREMX(const double &Depth, double &rho,double &vpv,double &vph,
           double &vsv,double &vsh,double &qu,double &qk,double &yita,
           const bool &PrintWarning=true) {

    PREM(Depth,1,0,rho,vpv,vph,vsv,vsh,qu,qk,yita,PrintWarning);

    double R=6371.0;
    double r=R-Depth;
    double x=r/R;


    if (Depth<=400) {
        vsh= -1147.226569441673*x*x*x + 3481.648723871561*x*x - 3521.617739418105*x + 1191.686592108216;
        vsv= vsh;
        rho=   734.780041069559*x*x*x - 2071.193615789281*x*x + 1938.047108369886*x -  598.252785440164;
        vph=   -4.9208884249226*x*x*x + 274.0496803031463*x*x - 533.3366953315248*x + 272.3185207233010;
        vpv= vph;
    }

    if (600<=Depth && Depth<=1155.674) {
        vsh= -84.46554942*x*x + 134.4361189*x - 46.95411628;
        vsv= vsh;
    }
    if (600<=Depth && Depth<=1726.323) {
        rho= -13.29902926*x*x + 16.06334012*x + 0.3373366175;
    }
    if (600<=Depth && Depth<=970.426) {
        vph= -189.1928756*x*x + 310.1340153*x - 115.5330402;
        vpv= vph;
    }

    return;
}

double DvsX(const double &Depth, const bool &PrintWarning=true) {
    double rho=0,vpv=0,vph=0,vsv=0,vsh=0,qu=0,qk=0,yita=0;
    PREMX(Depth,rho,vpv,vph,vsv,vsh,qu,qk,yita,PrintWarning);
    return vsv;
}

double DvpX(const double &Depth, const bool &PrintWarning=true) {
    double rho=0,vpv=0,vph=0,vsv=0,vsh=0,qu=0,qk=0,yita=0;
    PREMX(Depth,rho,vpv,vph,vsv,vsh,qu,qk,yita,PrintWarning);
    return vpv;
}

double DrhoX(const double &Depth, const bool &PrintWarning=true) {
    double rho=0,vpv=0,vph=0,vsv=0,vsh=0,qu=0,qk=0,yita=0;
    PREMX(Depth,rho,vpv,vph,vsv,vsh,qu,qk,yita,PrintWarning);
    return rho;
}

double RvsX(const double &Radius, const bool &PrintWarning=true) {
    return DvsX(6371.0-Radius,PrintWarning);
}

double RvpX(const double &Radius, const bool &PrintWarning=true) {
    return DvpX(6371.0-Radius,PrintWarning);
}

double RrhoX(const double &Radius, const bool &PrintWarning=true) {
    return DrhoX(6371.0-Radius,PrintWarning);
}

#endif
#ifndef ASU_PLANEWAVECOEFFICIENTS
#define ASU_PLANEWAVECOEFFICIENTS


/****************************************************************************************************
 * This C++ template returns plane wave reflection/refraction(transmission) coefficients. Reflection
 * related coefficients are denoted as "R", distinguishing from "T" for refraction(trasmission).
 *
 * A. For post-critical incident, reflection coefficients should be an imaginary number with
 *    amplitude 1. For the sake of consistant, will use complex<double> for all return values. Pre-
 *    critical incident will returns complex values with imaginary part 0.
 *
 * B. The Frequency Sign by default is negative to gives an evasive wave at post-critical conditions,
 *    whose amplitude decrease exponentially away from the interface.
 *
 * C. Ray always incident from media 1 (parameters related to index 1).
 *
 * D. In Aki & Richard, the meaning of angles (i,j) are related to wave types(P,S) for PSV. Here
 *    we use i for incident angle; j (with subscripts "_r", "_pp", "_ps", "_sp" and "_ss" for specific
 *    cases) to denote reflection or refraction(transmission) angles.
 *
 * E. Interfaces code names (SS,SA,SL,LS) (Media1 -> Media2).
 *
 *    Return (ans) order, for example, if cmp=="PSV" and mode=="SS":
 *    ans[0]=R_PP, ans[1]=R_PS, ...
 *
 *     Solid -> Solid  (SS):
 *
 *      PSV(8): R_PP, R_PS, R_SP, R_SS, T_PP, T_PS, T_SP, T_SS
 *      SH (2): R_SS, T_SS
 *
 *     Solid -> Liquid (SL):
 *
 *      PSV(6): R_PP, R_PS, R_SP, R_SS, T_PP, T_SP
 *      SH (1): R_SS (R_SS=1)
 *
 *     Solid -> Air    (SA):
 *
 *      PSV(2): R_PP, R_PS, R_SP, R_SS
 *      SH (1): R_SS (R_SS=1) [Notice: the phase shift introduced by internal caustic surface is not considered.]
 *
 *     Liquid -> Solid (LS):
 *
 *      PSV(3): R_PP, T_PP, T_PS
 *      SH (0):
 *
 *     Liquid -> Liquid (LL):
 *
 *      PSV(2): R_PP, T_PP
 *      SH (0):
 *
 *     Liquid -> Air (LA):
 *
 *      PSV(1): R_PP (R_PP=1)
 *      SH (0):
 *
 * F. input(s):
 *
 *     const double &rho1      ----  Density of media 1.
 *     const double &vp1       ----  P-wave velocity of media 1.
 *     const double &vs1       ----  S-wave velocity of media 1.
 *     const double &rho2      ----  Density of media 2.
 *     const double &vp2       ----  P-wave velocity of media 2.
 *     const double &vs2       ----  S-wave velocity of media 2.
 *     (If media2 is liquid or air, rho2, vp2, vs2 sometimes are ignored)
 *
 *     const double &inc       ----  Incident angle(in deg) in media 1.
 *     const string &cmp       ----  Component. "PSV" or "SH".
 *     const string &mode      ----  Structure mode. "SS", "SL", "SA", "LS", "LL", "LA".
 *     const int    &FreqSign  ----  (Optional, physics convention) Default value is -1, sign of frequencies.
 *
 * G. return(s):
 *
 *     vector<complex<double>> ans  ----  displacement coefficients. (see section E)
 *
 * Shule Yu
 * Apr 03 2018
 *
 * Key words: reflection, transmission, coefficient, plane wave, infinite frequency.
 *
 * Notes    : We are using negative freqeuncies so the imaginary part of vertical slowness is negative.
 *
 * Reference: Aki and Richards, Quantitative Seismology, 2nd ed., 2002.
 *            Page 140 (eq 5.34) - boundary conditions.
 *            Page 134,136,139,144,145 - results for "PSV-SA", "SH-SS", "PSV-SS"
****************************************************************************************************/

std::vector<std::complex<double>>
PlaneWaveCoefficients(const double &rho1,const double &vp1 ,const double &vs1,
                      const double &rho2,const double &vp2 ,const double &vs2,
                      const double &inc,const std::string &cmp,
                      const std::string &mode, const int &FreqSign=-1){

    double sini=sin(inc/180*M_PI),cosi=sqrt(1-sini*sini);

    if (cmp=="PSV") {

        if (mode=="SS") {

            // P as incident.
            double p=sini/vp1,sinj_r=p*vs1,sinj_pp=p*vp2,sinj_ps=p*vs2;

            std::complex<double> yp1=cosi/vp1,yp2(0),ys1(0),ys2(0);
            if (sinj_pp<=1) yp2=sqrt(1-sinj_pp*sinj_pp)/vp2;
            else yp2.imag(FreqSign*sqrt(p*p-1.0/vp2/vp2));

            if (sinj_r<=1) ys1=sqrt(1-sinj_r*sinj_r)/vs1;
            else ys1.imag(FreqSign*sqrt(p*p-1.0/vs1/vs1));

            if (sinj_ps<=1) ys2=sqrt(1-sinj_ps*sinj_ps)/vs2;
            else ys2.imag(FreqSign*sqrt(p*p-1.0/vs2/vs2));

            std::complex<double> a,b,c,d,E,F,G,H,D;
            a=rho2*(1-2*vs2*vs2*p*p)-rho1*(1-2*vs1*vs1*p*p);
            b=rho2*(1-2*vs2*vs2*p*p)+2*rho1*vs1*vs1*p*p;
            c=rho1*(1-2*vs1*vs1*p*p)+2*rho2*vs2*vs2*p*p;
            d=2*(rho2*vs2*vs2-rho1*vs1*vs1);

            E=b*yp1+c*yp2;
            F=b*ys1+c*ys2;
            G=a-d*yp1*ys2;
            H=a-d*yp2*ys1;

            D=E*F+G*H*p*p;

            std::complex<double> R_PP,R_PS,T_PP,T_PS;
            R_PP=(b*yp1-c*yp2)*F-(a+d*yp1*ys2)*H*p*p;
            R_PS=-2.0*yp1*(a*b+c*d*yp2*ys2)*p*vp1/vs1;
            T_PP=2.0*rho1*yp1*F*vp1/vp2;
            T_PS=2.0*rho1*yp1*H*p*vp1/vs2;

            R_PP/=D;
            R_PS/=D;
            T_PP/=D;
            T_PS/=D;


            // SV as incident.
            p=sini/vs1;sinj_r=p*vp1;
            double sinj_ss=p*vs2,sinj_sp=p*vp2;

            ys1=cosi/vs1;ys2=0;yp1=0;yp2=0;

            if (sinj_r<=1) yp1=sqrt(1-sinj_r*sinj_r)/vp1;
            else yp1.imag(FreqSign*sqrt(p*p-1.0/vp1/vp1));

            if (sinj_ss<=1) ys2=sqrt(1-sinj_ss*sinj_ss)/vs2;
            else ys2.imag(FreqSign*sqrt(p*p-1.0/vs2/vs2));

            if (sinj_sp<=1) yp2=sqrt(1-sinj_sp*sinj_sp)/vp2;
            else yp2.imag(FreqSign*sqrt(p*p-1.0/vp2/vp2));


            a=rho2*(1-2*vs2*vs2*p*p)-rho1*(1-2*vs1*vs1*p*p);
            b=rho2*(1-2*vs2*vs2*p*p)+2*rho1*vs1*vs1*p*p;
            c=rho1*(1-2*vs1*vs1*p*p)+2*rho2*vs2*vs2*p*p;
            d=2*(rho2*vs2*vs2-rho1*vs1*vs1);

            E=b*yp1+c*yp2;
            F=b*ys1+c*ys2;
            G=a-d*yp1*ys2;
            H=a-d*yp2*ys1;

            D=E*F+G*H*p*p;

            std::complex<double> R_SP,R_SS,T_SP,T_SS;
            R_SP=-2.0*ys1*(a*b+c*d*yp2*ys2)*p*vs1/vp1;
            R_SS=-(b*ys1-c*ys2)*E+(a+d*yp2*ys1)*G*p*p;
            T_SP=-2.0*rho1*ys1*G*p*vs1/vp2;
            T_SS=2.0*rho1*ys1*E*vs1/vs2;

            R_SP/=D;
            R_SS/=D;
            T_SP/=D;
            T_SS/=D;

            return {R_PP,R_PS,R_SP,R_SS,T_PP,T_PS,T_SP,T_SS};
        }
        else if (mode=="SL") {

            // P as incident.
            double p=sini/vp1,sinj_r=p*vs1,sinj_pp=p*vp2;

            std::complex<double> yp1=cosi/vp1,yp2(0),ys1(0);
            if (sinj_pp<=1) yp2=sqrt(1-sinj_pp*sinj_pp)/vp2;
            else yp2.imag(FreqSign*sqrt(p*p-1.0/vp2/vp2));

            if (sinj_r<=1) ys1=sqrt(1-sinj_r*sinj_r)/vs1;
            else ys1.imag(FreqSign*sqrt(p*p-1.0/vs1/vs1));

            std::complex<double> a,b,D;
            a=4*rho1*rho1*pow(vs1,4)*p*p*yp1*ys1;
            b=rho1*rho1*pow((1-2*vs1*vs1*p*p),2);
            D=(a+b)*yp2+rho1*rho2*yp1;

            std::complex<double> R_PP,R_PS,T_PP;
            R_PP=(a-b)*yp2+rho1*rho2*yp1;
            R_PS=2.0*rho1*rho1*vp1*vs1*(1-2*vs1*vs1*p*p)*p*yp1*yp2;
            T_PP=2.0*rho1*rho1*vp1/vp2*(1-2*vs1*vs1*p*p)*yp1;

            R_PP/=D;
            R_PS/=D;
            T_PP/=D;

            // S as incident.
            p=sini/vs1;sinj_r=p*vp1;
            double sinj_sp=p*vp2;

            ys1=cosi/vs1;yp1=0;yp2=0;
            if (sinj_r<=1) yp1=sqrt(1-sinj_r*sinj_r)/vp1;
            else yp1.imag(FreqSign*sqrt(p*p-1.0/vp1/vp1));

            if (sinj_sp<=1) yp2=sqrt(1-sinj_sp*sinj_sp)/vp2;
            else yp2.imag(FreqSign*sqrt(p*p-1.0/vp2/vp2));

            a=4*rho1*rho1*pow(vs1,4)*p*p*yp1*ys1;
            b=rho1*rho1*pow((1-2*vs1*vs1*p*p),2);
            D=(a+b)*yp2+rho1*rho2*yp1;

            std::complex<double> R_SP,R_SS,T_SP;
            R_SP=4*rho1*rho1*pow(vs1,3)/vp1*(1-2*vs1*vs1*p*p)*p*ys1*yp2;
            R_SS=(b-a)*yp2+rho1*rho2*yp1;
            T_SP=-4*rho1*rho1*pow(vs1,3)/vp2*p*ys1*yp1;

            R_SP/=D;
            R_SS/=D;
            T_SP/=D;

            return {R_PP,R_PS,R_SP,R_SS,T_PP,T_SP};
        }
        else if (mode=="SA") {

            // P as incident.
            double p=sini/vp1,sinj=p*vs1;

            std::complex<double> yp1=cosi/vp1,ys1(0);
            if (sinj<=1) ys1=sqrt(1-sinj*sinj)/vs1;
            else ys1.imag(FreqSign*sqrt(p*p-1.0/vs1/vs1));

            std::complex<double> A;
            A=pow(1.0/vs1/vs1-2*p*p,2)+4.0*p*p*yp1*ys1;

            std::complex<double> R_PP,R_PS;
            R_PP=-pow(1/vs1/vs1-2*p*p,2)+4*p*p*yp1*ys1;
            R_PS=4*vp1/vs1*p*yp1*(1.0/vs1/vs1-2*p*p);

            R_PP/=A;
            R_PS/=A;

            // SV as incident.
            p=sini/vs1;sinj=p*vp1;

            ys1=cosi/vs1;yp1=0;
            if (sinj<=1) yp1=sqrt(1-sinj*sinj)/vp1;
            else yp1.imag(FreqSign*sqrt(p*p-1.0/vp1/vp1));

            A=pow(1.0/vs1/vs1-2*p*p,2)+4.0*p*p*ys1*yp1;

            std::complex<double> R_SP,R_SS;
            R_SP=4*vs1/vp1*p*ys1*(1.0/vs1/vs1-2*p*p);
            R_SS=pow(1/vs1/vs1-2*p*p,2)-4*p*p*yp1*ys1;

            R_SP/=A;
            R_SS/=A;

            return {R_PP,R_PS,R_SP,R_SS};
        }
        else if (mode=="LS") {

            // P as incident.
            double p=sini/vp1,sinj_pp=p*vp2,sinj_ps=p*vs2;

            std::complex<double> yp1=cosi/vp1,yp2(0),ys2(0);
            if (sinj_pp<=1) yp2=sqrt(1-sinj_pp*sinj_pp)/vp2;
            else yp2=sqrt(p*p-1.0/vp2/vp2);

            if (sinj_ps<=1) ys2=sqrt(1-sinj_ps*sinj_ps)/vs2;
            else ys2=sqrt(p*p-1.0/vs2/vs2);

            std::complex<double> a,b,D;
            a=4*rho2*rho2*pow(vs2,4)*p*p*yp2*ys2;
            b=rho2*rho2*pow((1-2*vs2*vs2*p*p),2);
            D=(a+b)*yp1+rho1*rho2*yp2;

            std::complex<double> R_PP,T_PP,T_PS;
            R_PP=(a+b)*yp1-rho1*rho2*yp2;
            T_PP=2*rho1*rho2*(1-2*vs2*vs2*p*p)*vp1/vp2*yp1;
            T_PS=-4*rho1*rho2*vp1*vs2*p*yp1*yp2;

            R_PP/=D;
            T_PP/=D;
            T_PS/=D;

            return {R_PP,T_PP,T_PS};
        }
        else if (mode=="LL"){

            double p=sini/vp1,sinj=p*vp2;

            std::complex<double> yp1=cosi/vp1,yp2(0);

            if (sinj<=1) yp2=sqrt(1-sinj*sinj)/vp2;
            else yp2.imag(FreqSign*sqrt(p*p-1.0/vp2/vp2));

            std::complex<double> D(0);
            D=rho2*yp1+rho1*yp2;

            std::complex<double> R_PP,T_PP;
            R_PP=(rho2*yp1-rho1*yp2)/D;
            T_PP=2*rho1*vp1/vp2*yp1/D;

            return {R_PP,T_PP};
        }
        else if (mode=="LA") return {1};
    }
    else if (cmp=="SH"){

        if (mode=="SS") {

            double p=sini/vs1,sinj=p*vs2;

            std::complex<double> ys1=cosi/vs1,ys2(0);
            if (sinj<=1) ys2=sqrt(1.0-sinj*sinj)/vs2;
            else ys2.imag(FreqSign*sqrt(p*p-1.0/vs2/vs2));

            std::complex<double> A,B;
            A=rho1*vs1*vs1*ys1;
            B=rho2*vs2*vs2*ys2;

            std::complex<double> R_SS,T_SS;
            R_SS=(A-B)/(A+B);
            T_SS=2.0*A/(A+B);

            return {R_SS,T_SS};
        }
        else if (mode=="SA" || mode=="SL") return {1};
    }

    return {};
}

#endif
#ifndef ASU_POINTONSEGMENT
#define ASU_POINTONSEGMENT


/**************************************************************
 * This C++ template return whether given point is on a given
 * straight line segment. Due to rounding error of floating
 * number, the criteria is hard-wired that if the angle between
 * the point and two end point is larger than 179.99 deg)
 * [ cos (179.99 deg) = -0.999999984769129 ], the
 * point will be deemed as on the segment.
 *
 * input(s):
 * const pair<T1,T2> &p1   ----  {x1,y1} One end point of the line segment.
 * const pair<T3,T4> &p2   ----  {x2,y2} Another end point of the line segment.
 * const pair<T5,T6> &p    ----  Given point.
 *
 * return(s):
 * bool  ans  ----  true means the point is on the line segment.
 *
 * Shule Yu
 * Feb 02 2018
 *
 * Key words: point on line segment.
****************************************************************/

template<typename T1,typename T2,typename T3,typename T4,typename T5,typename T6>
bool PointOnSegment(const std::pair<T1,T2> &p1,const std::pair<T3,T4> &p2,const std::pair<T5,T6> &p){


    // First check if p lies on the straight line of p1p2.
    // p1p x p1p2
    double res1=(p.first-p1.first)*(p2.second-p1.second)-(p.second-p1.second)*(p2.first-p1.first);
    if (res1!=0) return false;

    double x1=p1.first,x2=p2.first;
    if (x1>x2) std::swap(x1,x2);

    double y1=p1.second,y2=p2.second;
    if (y1>y2) std::swap(y1,y2);

    return (x1<=p.first && p.first<=x2 && y1<=p.second && p.second<=y2);
}

#endif
#ifndef ASU_POINTINPOLYGON
#define ASU_POINTINPOLYGON



/**************************************************************
 * This C++ template decide weather given point on (x-y) plane
 * located inside one given polygon using the Winding Number
 * method.
 *
 * input(s):
 * const vector<pair<double,double>> &Polygon         ----  given polygon.
 * const pair<double,double>         &Point           ----  given point.
 * const int                         &BoundaryMode    ----  (Optional), default is 0
 *                                       0: left and bottom edges are decided to be inside.
 *                                          right and top edges are decided to be outside.
 *                                       1: point on all boundaries are counted as inside.
 *                                      -1: point on all boundaries are counted as outside.
 * const vector<double>              &PolygonBound    ----  (Optional), default is {DOUBLE_MIN,DOUBLE_MAX,DOUBLE_MIN,DOUBLE_MAX}
 *                                                          The given polygon bounds (Xmin, Xmax, Ymin, Ymax).
 * return(s):
 * bool  ans  ----  true means the point is inside the polygon.
 *
 * Notice: We define the inside of the polygon as the left hand
 * side area as one walks on its edge in the order of input.
 *
 * Shule Yu
 * Jan 18 2018
 *
 * Key words: point inside polygon, winding number.
 * Reference: http://geomalgorithms.com/a03-_inclusion.html
****************************************************************/

template<typename T1,typename T2,typename T3, typename T4>
bool PointInPolygon(const std::vector<std::pair<T1,T2>> &Polygon,const std::pair<T3,T4> &Point, const int BoundaryMode=0,
                    const std::vector<double> &PolygonBound={-std::numeric_limits<double>::max(),std::numeric_limits<double>::max(),
                                                             -std::numeric_limits<double>::max(),std::numeric_limits<double>::max()}){

    int n=(int)Polygon.size(),WN=0;
    double px=Point.first,py=Point.second;

    // Check bounds.
    if (px<PolygonBound[0] || px>PolygonBound[1] || py<PolygonBound[2] || py>PolygonBound[3]) return false;

    for (int i=0;i<n;++i){

        double ex1,ey1,ex2,ey2;

        ex1=Polygon[i].first;
        ey1=Polygon[i].second;
        ex2=Polygon[(i+1)%n].first;
        ey2=Polygon[(i+1)%n].second;

        auto flag=PointOnSegment(std::make_pair(ex1,ey1),std::make_pair(ex2,ey2),Point);

        if (BoundaryMode==1 && flag) return true;
        if (BoundaryMode==-1 && flag) return false;

        if (ey1<=py && py<ey2 && CrossProduct(ex1-px,ey1-py,0,ex2-px,ey2-py,0).back()>0) ++WN;
        else if (ey2<=py && py<ey1 && CrossProduct(ex1-px,ey1-py,0,ex2-px,ey2-py,0).back()<0) --WN;
    }

    return (WN!=0);
}

#endif
#ifndef ASU_RAYPATH
#define ASU_RAYPATH



/********************************************************************
 * This C++ template calculate the ray path for the given ray
 * parameters.
 *
 * input(s):
 * const vector<T1> &r             ----  layer radius array. r[0] is the shallowest layer.
 * const vector<T2> &v             ----  velocity array at each radius.
 * const double     &rayp          ----  ray parameter. (in sec/deg, p=Rsin/c/180*PI)
 * const double     &MinDepth      ----  depth to start the ray tracing.
 * const double     &MaxDepth      ----  depth to stop the ray tracing.
 * vector<double>   &degree        ----  Output ray path, angles. If degree[0]<-1e5, will only output degree.back();
 * std::size_t      &radius        ----  Output ray path end index in "r".
 * const double     &TurningAngle  ----  (Optional) the critical angle for turning assessment.
 *                                       default value is 89.9 deg.
 * return(s):
 * pair<pair<double,double>,bool>  ans  ----  {{travel time (s) / pursuit distance (km)}, RayTurns?}
 * vector<double> &degree (in-place)
 * std::size_t    &radius (in-place)
 *
 * Shule Yu
 * Jan 28 2018
 *
 * Key words: ray path, ray tracing, layers.
********************************************************************/

template<typename T1, typename T2>
std::pair<std::pair<double,double>,bool>
RayPath(const std::vector<T1> &r, const std::vector<T2> &v,
        const double &rayp, const double &MinDepth, const double &MaxDepth,
        std::vector<double> &degree, std::size_t &radius, const double &TurningAngle=89.9){

    // check inputs.
    double RE=6371.0;
    if (MaxDepth<=MinDepth || MaxDepth>RE-r.back() || MinDepth<RE-r[0]){
        return {{-1,-1},false};
    }

    auto cmp=[](const T1 &a, const T1 &b){
        return a>=b;
    };
    if (!std::is_sorted(r.begin(),r.end(),cmp)) {
        return {{-1,-1},false};
    }

    // locate our start Layer and end Layer.
    // because the grid could be unevenly spaced, this is a bit lengthy;
    std::size_t P1;
    double CurMin=std::numeric_limits<double>::max();
    for (P1=0;P1<r.size();++P1) {
        double NewMin=fabs(RE-MinDepth-r[P1]);
        if (CurMin<NewMin) {--P1;break;}
        CurMin=NewMin;
    }

    std::size_t P2;
    CurMin=std::numeric_limits<double>::max();
    for (P2=0;P2<r.size();++P2) {
        double NewMin=fabs(RE-MaxDepth-r[P2]);
        if (CurMin<NewMin) {--P2;break;}
        CurMin=NewMin;
    }
    if (P2==r.size()) --P2;

    // prepare output.
    bool OutPutDegree=(degree.empty() || degree[0]>=-1e5);
    degree.clear();

    // start ray tracing.
    //
    //   B,C are angles in the same layer, B=C+D.
    //   B=sin(incident_angle on current layer);
    //   C=sin(takeoff_angle from last layer);
    //   D=sin(trun_angle);

    double deg=0,MaxAngle=sin(TurningAngle*M_PI/180),Rayp=rayp*180/M_PI;
    std::pair<std::pair<double,double>,bool> ans{{0,0},false};
    for (std::size_t i=P1;i<P2;++i){

        double B,C,D;

        B=Rayp*v[i+1]/r[i+1];
        C=Rayp*v[i+1]/r[i];
        D=B*sqrt(1-C*C)-sqrt(1-B*B)*C;

        // Judge turning.
        if (C>=1 || B>1) {
            radius=i;
            degree.push_back(deg);
            ans.second=true;
            return ans;
        }

        double dist=r[i+1]/C*D;
        if (std::isnan(dist)) dist=LocDist(0,0,r[i],asin(D)*180/M_PI,0,r[i+1]);

        // store travel time and distance of this step.
        ans.first.first+=dist/v[i+1];
        ans.first.second+=dist;

        // store the path of this step.
        if (OutPutDegree) degree.push_back(deg);
        deg+=asin(D)*180/M_PI;

        // Judge turning.
        if (B>=MaxAngle) {
            radius=i+1;
            degree.push_back(deg);
            ans.second=true;
            return ans;
        }
    }
    radius=P2;
    degree.push_back(deg);

    return ans;
}

#endif
#ifndef ASU_SEGMENTJUNCTION
#define ASU_SEGMENTJUNCTION


/**************************************************************************************************
 * This C++ template return the junction point of two line segments.
 *
 * input(s):
 * const pair<T1,T2> &p1        ----  First end point of the first segment.
 * const pair<T3,T4> &q1        ----  Second end point of the first segment.
 * const pair<T5,T6> &p2        ----  First end point of the second segment.
 * const pair<T7,T8> &q2        ----  Second end point of the second segment.
 * const bool        &EndPoint  ----  (Optional) End points counted as a junction? Default is true.
 *                                      true : return true for junction happens at any end ponts.
 *                                      false: return false for junction happens at any end ponts.
 *
 * return(s):
 * pair<bool,pair<double,double>>  ans  ----  {flag, position of the junction point}
 *                                            If no junction point, flag is set to false.
 *
 * Shule Yu
 * Jan 19 2018
 *
 * Key words: segment junction, intercept, cross.
**************************************************************************************************/

template<typename T1,typename T2,typename T3,typename T4,
         typename T5,typename T6,typename T7,typename T8>
std::pair<bool,std::pair<double,double>>
SegmentJunction(const std::pair<T1,T2> &p1, const std::pair<T3,T4> &q1,
                const std::pair<T5,T6> &p2, const std::pair<T7,T8> &q2, const bool &EndPoint=true){
    /*
         q2
         
     p1      q1
         
         p2

    */

    // which side is p1, relative to p2q2
    // p2q2 x p2p1
    double res1=(q2.first-p2.first)*(p1.second-p2.second)-(q2.second-p2.second)*(p1.first-p2.first);

    // which side is q1, relative to p2q2
    // p2q2 x p2q1
    double res2=(q2.first-p2.first)*(q1.second-p2.second)-(q2.second-p2.second)*(q1.first-p2.first);

    // which side is p2, relative to p1q1
    // p1q1 x p1p2
    double res3=(q1.first-p1.first)*(p2.second-p1.second)-(q1.second-p1.second)*(p2.first-p1.first);

    // which side is q2, relative to p1q1
    // p1q1 x p1q2
    double res4=(q1.first-p1.first)*(q2.second-p1.second)-(q1.second-p1.second)*(q2.first-p1.first);

    // Check whether p1q1 and p2q2 straddles each other.
    double minX1=p1.first,maxX1=q1.first,minY1=p1.second,maxY1=q1.second;
    if (minX1>maxX1) std::swap(minX1,maxX1);
    if (minY1>maxY1) std::swap(minY1,maxY1);
    double minX2=p2.first,maxX2=q2.first,minY2=p2.second,maxY2=q2.second;
    if (minX2>maxX2) std::swap(minX2,maxX2);
    if (minY2>maxY2) std::swap(minY2,maxY2);

    if (((res1>0 && res2<0) || (res1<0 && res2>0)) && ((res3>0 && res4<0) || (res3<0 && res4>0))) ;
    else if (res1==0 && minX2<=p1.first && p1.first<=maxX2 && minY2<=p1.second && p1.second<=maxY2) return {true,p1};
    else if (res2==0 && minX2<=q1.first && q1.first<=maxX2 && minY2<=q1.second && q1.second<=maxY2) return {true,q1};
    else if (res3==0 && minX1<=p2.first && p2.first<=maxX1 && minY1<=p2.second && p2.second<=maxY1) return {true,p2};
    else if (res4==0 && minX1<=q2.first && q2.first<=maxX1 && minY1<=q2.second && q2.second<=maxY1) return {true,q2};
    else return {false,{}};

    // Find the junction.
    auto res=LineJunction(p1,(p1.second-q1.second)/(p1.first-q1.first),p2,(p2.second-q2.second)/(p2.first-q2.first));
    return {true,{res.second}};
}

#endif


using namespace std;

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
void followThisRay(
    size_t i, atomic<size_t> &Cnt, atomic<int> &Estimation, atomic<size_t> &Running,
    char **ReachSurfaces, int *ReachSurfacesSize, char **RayInfo, int *RayInfoSize,
    double **RaysTheta, int *RaysN, double **RaysRadius,
    vector<Ray> &RayHeads, int branches, const vector<double> &specialDepths,
    const vector<vector<double>> &R, const vector<vector<double>> &Vp,
    const vector<vector<double>> &Vs,const vector<vector<double>> &Rho,
    const vector<vector<pair<double,double>>> &Regions, const vector<vector<double>> &RegionBounds,
    const vector<double> &dVp, const vector<double> &dVs,const vector<double> &dRho,
    const bool &DebugInfo,const bool &TS,const bool &TD,const bool &RS,const bool &RD, const bool &StopAtSurface){

    if (RayHeads[i].RemainingLegs==0) {
        Running.fetch_sub(1);
        return;
    }


    // Locate the begining and ending depths for the next leg.

    /// ... among special depths.

    //// Which special depth is cloest to ray head depth?
    double RayHeadDepth=_RE-RayHeads[i].Pr;
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
    Top=max(Top,_RE-RegionBounds[CurRegion][3]);
    Bot=min(Bot,_RE-RegionBounds[CurRegion][2]);

    // Print some debug info.
    if (DebugInfo) {
        RayHeads[i].Debug+=to_string(1+i)+" --> ";
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

                RayEnd=(int)j;

                // which region is the new leg entering?
                for (size_t k=1;k<Regions.size();++k){
                    if ((int)k==CurRegion) continue;
                    if (PointInPolygon(Regions[k],p,1,RegionBounds[k])) {
                        NextRegion=(int)k;
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
                    RayEnd=(int)j;
                    NextRegion=(int)k;
                    break;
                }
            }
            if (RayEnd!=-1) break; // If ray enters another region.
            else NextRegion=0;
        }
    }


    // Print some debug info.
    if (DebugInfo) {
    }


    // Prepare reflection/refraction flags. (notice "rs" [r]eflection to [s]ame wave type is always possible)
    bool ts=TS,td=(TD && RayHeads[i].Comp!="SH"),rd=(RD && RayHeads[i].Comp!="SH");


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
        size_t L1=0,L2=1,SearchRegion=(NextRegion==0?CurRegion:NextRegion);
        p2={NextPt_R,NextPr_R};
        q2={NextPt_T,NextPr_T};
        pair<bool,pair<double,double>> res;


        for (L1=0;L1<Regions[SearchRegion].size();++L1){ // Search around the polygon.
            L2=(L1+1)%Regions[SearchRegion].size();
            res=SegmentJunction(Regions[SearchRegion][L1],Regions[SearchRegion][L2],p2,q2);
            if (res.first) break;
        }
        if (L1==Regions[SearchRegion].size()) {
        }

        // Find the junction point between ray and polygon boundary.
        JuncPt=res.second.first;
        JuncPr=res.second.second;

        // Print some debug info.


        // Twick travel times and travel distance, compensate for the lost part.
        double dlx=(p2.first-JuncPt)*M_PI*JuncPr/180,dly=p2.second-JuncPr;
        double dl=sqrt(dlx*dlx+dly*dly);
        ans.first.second+=dl;
        ans.first.first+=dl/v[CurRegion][rIndex(RayEnd-1)]; // Use the velocit within current region to avoid possible "inf" travel time.


        // Get the geometry of the boundary.
        const pair<double,double> &p1=Regions[SearchRegion][L1],&q1=Regions[SearchRegion][L2];
        TiltAngle=180/M_PI*atan2(q1.second-p1.second,(q1.first-p1.first)*M_PI/180*JuncPr);

    }
    else { // If ray doesn't end pre-maturelly (stays in the same region and reflect/refract on horizontal intervals)
        // (one end point of the last line segment (index: RayEnd-1) is on the interface)

        RayEnd=(int)degree.size();
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
    if (DebugInfo) {
    }


    // Prepare to calculate reflection/refractoin(transmission) coefficients.
    string Mode,Polarity=(RayHeads[i].Comp=="SH"?"SH":"PSV");
    if (NextPr_R==_RE) Mode="SA"; // At the surface.
    else if (NextPr_R==3480) Mode=(RayHeads[i].GoUp?"LS":"SL"); // At the CMB.
    else if (NextPr_R==1221.5) Mode=(RayHeads[i].GoUp?"SL":"LS"); // At the ICB.
    else if (1221.5<NextPr_R && NextPr_R<3480) Mode="LL"; // Within outer core.
    else Mode="SS";

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
    if (DebugInfo) {

    }

    // Update current RayHead.
    RayHeads[i].TravelTime=ans.first.first;
    RayHeads[i].TravelDist=ans.first.second;
    --RayHeads[i].RemainingLegs;
    RayHeads[i].Inc=Incident;


    // store ray paths.
    stringstream ss;
    ss << RayHeads[i].Color << " "
       << (RayHeads[i].IsP?"P ":"S ") << RayHeads[i].TravelTime << " sec. " << RayHeads[i].Inc << " IncDeg. "
       << RayHeads[i].Amp << " DispAmp. " << RayHeads[i].TravelDist << " km. ";
    string tmpstr=ss.str();
    RayInfoSize[i]=(int)tmpstr.size()+1;
    RayInfo[i]=(char *)malloc((tmpstr.size()+1)*sizeof(char));
    strcpy(RayInfo[i],tmpstr.c_str());

    RaysN[i]=RayEnd;
    RaysTheta[i]=(double *)malloc(RayEnd*sizeof(double));
    RaysRadius[i]=(double *)malloc(RayEnd*sizeof(double));
    for (int j=0;j<RayEnd;++j) {
        RaysTheta[i][j]=RayHeads[i].Pt+M*degree[j];
        RaysRadius[i][j]=R[CurRegion][rIndex(j)];
    }

    // If ray reaches surface, output info at the surface.
    if (NextPr_R==_RE) ++RayHeads[i].Surfacing;
    if (NextPr_R==_RE && (StopAtSurface==0 || RayHeads[i].Surfacing<2)) {

        // Accumulate the travel-time.
        int I=(int)i;
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
            ss << (1+*rit) << ((*rit)==*hh.begin()?"":"->");

        string tmpstr=ss.str();
        if (!tmpstr.empty()) {
            ReachSurfacesSize[i]=(int)tmpstr.size()+1;
            ReachSurfaces[i]=(char *)malloc((tmpstr.size()+1)*sizeof(char));
            strcpy(ReachSurfaces[i],tmpstr.c_str());
        }

        if (StopAtSurface==1) {
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
    if (NextPr_R==_RE) ts=td=false;

    /// if ray goes down to CMB, no transmission S; if ray goes up to ICB, not transmission S.
    if ((!RayHeads[i].GoUp && NextPr_T==3480) || (RayHeads[i].GoUp && NextPr_T==1221.5)) {
        ts&=RayHeads[i].IsP;
        td&=(!RayHeads[i].IsP);
    }

    /// if ray goes down to ICB as P, no reflection S.
    if (!RayHeads[i].GoUp && RayHeads[i].IsP && NextPr_R==1221.5) rd=false;

    /// if ray goes up to CMB as P, no reflection S.
    if (RayHeads[i].GoUp && RayHeads[i].IsP && NextPr_R==3480) rd=false;

    // Add new ray heads to "RayHeads" according to the rules ans reflection/refraction angle calculation results.

    if (ts) {
        Ray newRay=RayHeads[i];
        newRay.Prev=(int)i;
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
        newRay.Prev=(int)i;
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
        newRay.Prev=(int)i;
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

    int y=(TS && !ts)+(TD && !td)+(RD && !rd);
    int z=RayHeads[i].RemainingLegs;
    if (branches>1) z=(1-pow(branches,RayHeads[i].RemainingLegs))/(1-branches);
    Estimation.fetch_sub(y*z);

    // rs is always possible.
    if (RS) {
        Ray newRay=RayHeads[i];
        newRay.Prev=(int)i;
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

void PreprocessAndRun(
        const vector<int> &initRaySteps,const vector<int> &initRayComp,const vector<int> &initRayColor,
        const vector<double> &initRayTheta,const vector<double> &initRayDepth,const vector<double> &initRayTakeoff,
        const vector<double> &gridDepth1,const vector<double> &gridDepth2,const vector<double> &gridInc,
        const vector<double> &specialDepths,const vector<vector<double>> &Deviation,
        const vector<vector<double>> &regionProperties,
        const vector<vector<double>> &regionPolygonsTheta,
        const vector<vector<double>> &regionPolygonsDepth,

        const double &RectifyLimit, const bool &TS, const bool &TD, const bool &RS, const bool &RD,
        const size_t &nThread, const bool &DebugInfo, const bool &StopAtSurface,
        const size_t &branches, const size_t &potentialSize,

        char **ReachSurfaces, int *ReachSurfacesSize, char **RayInfo, int *RayInfoSize,
        int *RegionN,double **RegionsTheta,double **RegionsRadius,
        double **RaysTheta, int *RaysN, double **RaysRadius,int *Observer){

    // Ray nodes.
    vector<Ray> RayHeads;
    if (potentialSize>RayHeads.max_size()) {
    }

    // Create 1D reference layers. (R[0]. 0 means 1D reference model)

    vector<vector<double>> R(1);
    for (size_t i=0;i<gridDepth1.size();++i){
        auto tmpr=CreateGrid(_RE-gridDepth2[i],_RE-gridDepth1[i],gridInc[i],2);
        if (!R[0].empty()) R[0].pop_back();
        R[0].insert(R[0].end(),tmpr.rbegin(),tmpr.rend());
    }


    // Fix round-off-errors:
    // adding the exact double values in A. special depths and B. modefied 1D model to R[0].
    set<double> depthToCorrect(specialDepths.begin(),specialDepths.end());
    for (const auto &item:Deviation) {
        depthToCorrect.insert(item[0]);
        depthToCorrect.insert(item[1]);
    }
    vector<double> tmpArray;
    swap(R[0],tmpArray);
    tmpArray[0]=_RE;tmpArray.back()=0;
    auto it=depthToCorrect.begin();
    for (int i=0;i<(int)tmpArray.size();++i) {
        if (it==depthToCorrect.end())
            R[0].push_back(tmpArray[i]);
        else if (tmpArray[i]==_RE-*it) {
            R[0].push_back(tmpArray[i]);
            ++it;
        }
        else if (tmpArray[i]>_RE-*it) {
            R[0].push_back(tmpArray[i]);
        }
        else {
            R[0].push_back(_RE-*it);
            ++it;
            --i;
        }
    }
    while (it!=depthToCorrect.end()){
        R[0].push_back(_RE-*it);
        ++it;
    }

    // Find the bounds of input polygons.
    vector<vector<double>> RegionBounds{{-numeric_limits<double>::max(),numeric_limits<double>::max(),
        -numeric_limits<double>::max(),numeric_limits<double>::max()}};
    // the 1D reference bounds is as large as possible.

    for (size_t i=0;i<regionPolygonsTheta.size();++i){
        double Xmin=numeric_limits<double>::max(),Xmax=-Xmin,Ymin=Xmin,Ymax=-Ymin;
        for (size_t j=0;j<regionPolygonsTheta[i].size();++j){
            size_t k=(j+1)%regionPolygonsTheta[i].size();
            double theta1=regionPolygonsTheta[i][j],theta2=regionPolygonsTheta[i][k];
            double radius1=_RE-regionPolygonsDepth[i][j],radius2=_RE-regionPolygonsDepth[i][k];

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
            double radius1=_RE-regionPolygonsDepth[i][j],radius2=_RE-regionPolygonsDepth[i][k];
            double theta1=regionPolygonsTheta[i][j],theta2=regionPolygonsTheta[i][k];

            if (radius1==RegionBounds[i+1][2]) radius1=R[0][adjustedYmin[i+1]];
            if (radius1==RegionBounds[i+1][3]) radius1=R[0][adjustedYmax[i+1]];
            if (radius2==RegionBounds[i+1][2]) radius2=R[0][adjustedYmin[i+1]];
            if (radius2==RegionBounds[i+1][3]) radius2=R[0][adjustedYmax[i+1]];

            double Tdist=theta2-theta1,Rdist=radius2-radius1;

            size_t NPTS=2;
            double dL=_RE,dR=Rdist,dT=Tdist;
            while (dL>RectifyLimit){
                NPTS*=2;
                dR=Rdist/(NPTS-1);
                dT=Tdist/(NPTS-1);
                dL=LocDist(theta1,0,radius1,theta1+dT,0,radius1+dR);
            }

            // Add rectified segments to this polygon.
            for (size_t k=0;k+1<NPTS;++k)
                tmpRegion.push_back(make_pair(theta1+k*dT,radius1+k*dR));
        }

        // Add this rectified polygon to region array.
        Regions.push_back(tmpRegion);

        RegionN[i]=(int)tmpRegion.size();
        RegionsTheta[i]=(double *)malloc(tmpRegion.size()*sizeof(double));
        RegionsRadius[i]=(double *)malloc(tmpRegion.size()*sizeof(double));
        for (int j=0;j<RegionN[i];++j){
            RegionsTheta[i][j]=tmpRegion[j].first;
            RegionsRadius[i][j]=tmpRegion[j].second;
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

    // derive properties layers for all regions (include 1D reference).
    vector<vector<double>> Vp(R.size(),vector<double> ()),Vs=Vp,Rho=Vp;
    for (size_t i=0;i<R.size();++i) {
        for (const auto &item:R[i]) {
            auto ans=MakeRef(_RE-item,Deviation);
            Vp[i].push_back(dVp[i]*ans[0]);
            Vs[i].push_back(dVs[i]*ans[1]);
            Rho[i].push_back(dRho[i]*ans[2]);
        }
    }

    // Create initial rays.
    for (size_t i=0;i<initRaySteps.size();++i){

        // Source in any polygons?
        size_t rid=0;
        for (size_t i=1;i<Regions.size();++i)
            if (PointInPolygon(Regions[i],make_pair(initRayTheta[i],_RE-initRayDepth[i]),1,RegionBounds[i])) {rid=i;break;}

        // Calculate ray parameter.
        auto ans=MakeRef(initRayDepth[i],Deviation);
        double v=(initRayComp[i]==0?ans[0]*dVp[rid]:ans[1]*dVs[rid]);
        double rayp=M_PI/180*(_RE-initRayDepth[i])*sin(fabs(initRayTakeoff[i])/180*M_PI)/v;

        // Push this ray into "RayHeads" for future processing.
        RayHeads.push_back(Ray(initRayComp[i]==0,fabs(initRayTakeoff[i])>=90,initRayTakeoff[i]<0,
                    (initRayComp[i]==0?"P":(initRayComp[i]==1?"SV":"SH")),
                    (int)rid,initRaySteps[i],initRayColor[i],
                    initRayTheta[i],_RE-initRayDepth[i],0,0,rayp,initRayTakeoff[i]));
    }

    atomic<size_t> Cnt;
    Cnt.store(RayHeads.size());
    atomic<int> Estimation;
    Estimation.store((int)potentialSize);
    RayHeads.resize(potentialSize);

    // Start ray tracing. (Finally!)
    //
    // Process each "Ray" leg in "RayHeads".
    // For future legs generated by reflction/refraction, create new "Ray" and assign it to the proper position in "RayHeads" vector.
    vector<thread> allThread(RayHeads.capacity());

    size_t Doing=0,Done=0;
    atomic<size_t> Running;
    Running.store(0);
    while (Doing!=Cnt.load() || Running.load()!=0) {
        *Observer=(int)Doing-(int)nThread;

        if (Running.load()<nThread && Doing<Cnt.load()) {
            Running.fetch_add(1);
            allThread[Doing]=thread(
                followThisRay,Doing, std::ref(Cnt), std::ref(Estimation), std::ref(Running),
                ReachSurfaces,ReachSurfacesSize,RayInfo,RayInfoSize, RaysTheta, RaysN, RaysRadius,
                std::ref(RayHeads), branches, std::cref(specialDepths),
                std::cref(R),std::cref(Vp),std::cref(Vs),std::cref(Rho),
                std::cref(Regions),std::cref(RegionBounds),std::cref(dVp),std::cref(dVs),std::cref(dRho),
                std::cref(DebugInfo),std::cref(TS),std::cref(TD),std::cref(RS),std::cref(RD),std::cref(StopAtSurface));
            if (Doing>0 && Doing%10000==0) {
                for (size_t i=Done;i<Doing-nThread;++i)
                    allThread[i].join();
                Done=Doing-nThread;
            }
            ++Doing;
        }
        else usleep(1000);

    } // End of ray tracing.

    for (size_t i=Done;i<Doing;++i)
        allThread[i].join();

}

void rayTracingInSwift(
    int inputRayN, int *inputRaySteps, int *inputRayComp, int *inputRayColor,
    double *inputRayTheta, double *inputRayDepth, double *inputRayTakeoff,
    int inputGridN, double *inputGridDepth1, double *inputGridDepth2, double *inputGridInc,
    int inputDepthN, double *inputSpecialDepths,
    int inputDeviationN, double **inputDeviation,
    int inputRegionN, double **inputRegionProperties,
    int *inputRegionL, double **inputRegionPolygonsTheta, double **inputRegionPolygonsDepth,
    double inputRectifyLimit, bool inputTS, bool inputTD, bool inputRS, bool inputRD, bool inputStopAtSurface,
    int inputNThread,
    char ***ReachSurfaces, int *ReachSurfacesSize, char ***RayInfo, int *RayInfoSize,
    int *RegionN, double **RegionsTheta,double **RegionsRadius,
    double **RaysTheta, int *RaysN, double **RaysRadius,int *Observer){


    // Bridging variables for the C++ code.
    vector<int> initRaySteps,initRayComp,initRayColor;
    vector<double> initRayTheta,initRayDepth,initRayTakeoff,gridDepth1,gridDepth2,gridInc,specialDepths;
    vector<vector<double>> Deviation,regionProperties,regionPolygonsTheta,regionPolygonsDepth;

    for (int i=0;i<inputRayN;++i) {
        initRaySteps.push_back(inputRaySteps[i]);
        initRayComp.push_back(inputRayComp[i]);
        initRayColor.push_back(inputRayColor[i]);
        initRayTheta.push_back(inputRayTheta[i]);
        initRayDepth.push_back(inputRayDepth[i]);
        initRayTakeoff.push_back(inputRayTakeoff[i]);
    }

    for (int i=0;i<inputGridN;++i) {
        gridDepth1.push_back(inputGridDepth1[i]);
        gridDepth2.push_back(inputGridDepth2[i]);
        gridInc.push_back(inputGridInc[i]);
    }

    for (int i=0;i<inputDepthN;++i)
        specialDepths.push_back(inputSpecialDepths[i]);

    for (int i=0;i<inputDeviationN;++i) {
        Deviation.push_back(vector<double> ());
        for (int j=0;j<5;++j)
            Deviation.back().push_back(inputDeviation[i][j]);
    }

    for (int i=0;i<inputRegionN;++i){
        regionProperties.push_back(vector<double> {inputRegionProperties[i][0],inputRegionProperties[i][1],inputRegionProperties[i][2]});
        regionPolygonsTheta.push_back(vector<double> ());
        regionPolygonsDepth.push_back(vector<double> ());
        for (int j=0;j<inputRegionL[i];++j) {
            regionPolygonsTheta.back().push_back(inputRegionPolygonsTheta[i][j]);
            regionPolygonsDepth.back().push_back(inputRegionPolygonsDepth[i][j]);
        }
    }

    double RectifyLimit=inputRectifyLimit;
    bool TS=inputTS,TD=inputTD,RS=inputRS,RD=inputRD,DebugInfo=false,StopAtSurface=inputStopAtSurface;
    size_t nThread=(size_t)inputNThread,potentialSize=0;
    int branches=TS+TD+RS+RD;
    for (size_t i=0;i<initRaySteps.size();++i){
        if (branches<=1) potentialSize+=initRaySteps[i];
        else potentialSize+=(1-pow(branches,initRaySteps[i]))/(1-branches);
    }

    // Spaces for the outputs.
    Observer=(int *)malloc(1*sizeof(int));
    *Observer=-1;

    RegionsTheta=(double **)malloc(regionProperties.size()*sizeof(double *));
    RegionsRadius=(double **)malloc(regionProperties.size()*sizeof(double *));
    RaysTheta=(double **)malloc(potentialSize*sizeof(double *));
    RaysRadius=(double **)malloc(potentialSize*sizeof(double *));
    *ReachSurfaces=(char **)malloc(potentialSize*sizeof(char *));
    *RayInfo=(char **)malloc(potentialSize*sizeof(char *));

    RaysN=(int *)malloc(potentialSize*sizeof(int));
    ReachSurfacesSize=(int *)malloc(potentialSize*sizeof(int));
    RayInfoSize=(int *)malloc(potentialSize*sizeof(int));
    for (int i=0;i<(int)potentialSize;++i) {
        RaysN[i]=0;
        ReachSurfacesSize[i]=0;
        RayInfoSize[i]=0;
    }
    RegionN=(int *)malloc(regionProperties.size()*sizeof(int));
    for (size_t i=0;i<regionProperties.size();++i) RegionN[i]=0;

    // Call the C++ code.
    PreprocessAndRun(
        initRaySteps,initRayComp,initRayColor,
        initRayTheta,initRayDepth,initRayTakeoff,gridDepth1,gridDepth2,gridInc,specialDepths,
        Deviation,regionProperties,regionPolygonsTheta,regionPolygonsDepth,
        RectifyLimit,TS,TD,RS,RD,nThread,DebugInfo,StopAtSurface,(size_t)branches,potentialSize,
        *ReachSurfaces,ReachSurfacesSize,*RayInfo,RayInfoSize,RegionN,RegionsTheta,RegionsRadius,RaysTheta,RaysN,RaysRadius,Observer);
}

