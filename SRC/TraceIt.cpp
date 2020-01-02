#include<Ray.hpp>

#include<ReadParameters.hpp>

using namespace std;

// The main function mostly dealt with I/O.
int main(int argc, char **argv){

    enum PI{DebugInfo,TS,TD,RS,RD,StopAtSurface,nThread,FLAG1};
    enum PS{InputRays,Layers,Depths,Ref,Polygons,ReceiverFileName,PolygonFilePrefix,RayFilePrefix,FLAG2};
    enum PF{RectifyLimit,FLAG3};

    auto P=ReadParameters<PI,PS,PF> (argc,argv,cin,FLAG1,FLAG2,FLAG3);

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
    if (!regionPolygonsTheta.empty() && regionPolygonsTheta.back().empty()) {
        regionPolygonsTheta.pop_back();
        regionPolygonsDepth.pop_back();
    }
    fpin.close();


    // I/O is Done.
    //
    // Currently we have these variables ------ :
    // ReadParameters<PI,PS,PF> P;
    // vector<int> initRaySteps,initRayComp,initRayColor;
    // vector<double> initRayTheta,initRayDepth,initRayTakeoff,gridDepth1,gridDepth2,gridInc,specialDepths;
    // vector<vector<double>> Deviation,regionProperties,regionPolygonsTheta,regionPolygonsDepth;
    //
    // For future I/O modification, you can start from begining and stop here.


    // Estimate the output size and malloc spaces.
    size_t potentialSize=0;
    int branches=(P[TS]+P[TD]+P[RS]+P[RD]);
    for (size_t i=0;i<initRaySteps.size();++i){
        if (branches<=1) potentialSize+=initRaySteps[i];
        else potentialSize+=(1-pow(branches,initRaySteps[i]))/(1-branches);
    }

    int *Observer=(int *)malloc(1*sizeof(int));
    *Observer=-1;

    double **RegionsTheta=(double **)malloc(regionProperties.size()*sizeof(double *));
    double **RegionsRadius=(double **)malloc(regionProperties.size()*sizeof(double *));
    double **RaysTheta=(double **)malloc(potentialSize*sizeof(double *));
    double **RaysRadius=(double **)malloc(potentialSize*sizeof(double *));
    char **ReachSurfaces=(char **)malloc(potentialSize*sizeof(char *));
    char **RayInfo=(char **)malloc(potentialSize*sizeof(char *));

    int *RaysN=(int *)malloc(potentialSize*sizeof(int));
    int *ReachSurfacesSize=(int *)malloc(potentialSize*sizeof(int));
    int *RayInfoSize=(int *)malloc(potentialSize*sizeof(int));
    for (int i=0;i<(int)potentialSize;++i) {
        RaysN[i]=0;
        ReachSurfacesSize[i]=0;
        RayInfoSize[i]=0;
    }
    int *RegionN=(int *)malloc(regionProperties.size()*sizeof(int));
    for (size_t i=0;i<regionProperties.size();++i) RegionN[i]=0;


    PreprocessAndRun(
        initRaySteps,initRayComp,initRayColor,
        initRayTheta,initRayDepth,initRayTakeoff,gridDepth1,gridDepth2,gridInc,specialDepths,
        Deviation,regionProperties,regionPolygonsTheta,regionPolygonsDepth,
        P[RectifyLimit],(P[TS]!=0),(P[TD]!=0),(P[RS]!=0),(P[RD]!=0),(size_t)P[nThread],(P[DebugInfo]!=0),(P[StopAtSurface]!=0),
        branches,potentialSize,
        ReachSurfaces,ReachSurfacesSize,RayInfo,RayInfoSize,RegionN,RegionsTheta,RegionsRadius,RaysTheta,RaysN,RaysRadius,Observer);


    // Outputs.
    // If I/O changes, change this part.

    ofstream fpout(P[ReceiverFileName]);
    fpout << "<Takeoff> <Rayp> <Incident> <Dist> <TravelTime> <DispAmp> <RemainingLegs> <rayTurns> <WaveTypeTrain> <RayTrain>" << '\n';
    for (size_t i=0;i<potentialSize;++i)
        if (ReachSurfacesSize[i]!=0)
            fpout << string(ReachSurfaces[i]) << '\n';
    fpout.close();

    // Output valid part ray paths.
    if (P[RayFilePrefix]!="NONE") {
        for (size_t i=0;i<potentialSize;++i) {
            if (RayInfoSize[i]==0) continue;

            ofstream fpout(P[RayFilePrefix]+to_string(i+1));
            fpout << "> " << string(RayInfo[i]) << '\n';
            for (int j=0;j<RaysN[i];++j)
                fpout << RaysTheta[i][j] << " " << RaysRadius[i][j] << '\n';
            fpout.close();
        }
    }

    // Output.
    // For plotting: 1D reference property deviation depths.
    if (!Deviation.empty() && P[PolygonFilePrefix]!="NONE"){
        ofstream fpout(P[PolygonFilePrefix]+"0");
        for (const auto &item:Deviation) {
            fpout << ">\n";
            for (double t=0;t<360;t=t+0.1) fpout << t << " " << _RE-item[0] << '\n';
            fpout << ">\n";
            for (double t=0;t<360;t=t+0.1) fpout << t << " " << _RE-item[1] << '\n';
        }
        fpout.close();
    }

    // Output rectified regions.
    if (P[PolygonFilePrefix]!="NONE"){
        for (size_t i=0;i<regionProperties.size();++i) {
            FILE *fpout=fopen((P[PolygonFilePrefix]+to_string(i+1)).c_str(),"w");
            for (int j=0;j<RegionN[i];++j)
                fprintf(fpout,"%.7e %.7e\n",RegionsTheta[i][j],RegionsRadius[i][j]);
            fclose(fpout);
        }
    }

    // Free spaces.
    for (size_t i=0;i<potentialSize;++i){
        if (RaysN[i]!=0) {
            free(RaysTheta[i]);
            free(RaysRadius[i]);
        }
        if (ReachSurfacesSize[i]!=0) free(ReachSurfaces[i]);
        if (RayInfoSize[i]!=0) free(RayInfo[i]);
    }
    for (size_t i=0;i<regionProperties.size();++i) {
        if (RegionN[i]!=0) {
            free(RegionsTheta[i]);
            free(RegionsRadius[i]);
        }
    }
    free(RegionN);
    free(RegionsTheta);
    free(RegionsRadius);
    free(RaysN);
    free(RaysTheta);
    free(RaysRadius);
    free(ReachSurfacesSize);
    free(ReachSurfaces);
    free(RayInfoSize);
    free(RayInfo);

    return 0;
}
