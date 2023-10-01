/*
 *  M.o.r.p.h.o.l.o.g.i.e.s.
 *  ------------------------
 *
 *  Monte Carlo Methods used to Move Particles
 *
 ************************************************
 *
 * McMaster Engineering Physics
 * Written by: Matt Bumstead
 *
 * File generated: 2015.10.20
 * Last modified:  2017.02.08 4
 *
 ***********************************************/
//=================================
// include guard
#ifndef morphologies_h
#define morphologies_h

//=================================
// include dependancies
#include <stdio.h>
#include <vector>
#include <numeric>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

//=================================
// include Morphologies Script


#include "initialization.h"
#include "particles.h"
#include "shake.h"
//#include "boundary.h"



//=================================
// include guard
//=================================
// include dependancies
//=================================
// include Morphologies Script


//  Global Definitions
// --------------------
#define PI     3.141592653589793238462643


// #define DEBUG
#define ERR_NONE		0	// Successful return code
#define ERR_BADPARAM	1	// Bad parameter
#define ERR_BADPF		2	// Bad polygon file
#define ERR_BADINITS	3	// Error during creating of output directories
#define ERR_CREATEDIRS	4	// Error during creating of output directories

extern void PrintError(int errorNumber)
{
    cout << "Stopping Simulation: ";
    switch(errorNumber){
        case ERR_BADPARAM : cout    << "Bad Parameter Found! \n";       break;
        case ERR_BADPF : cout       << "Cannot Load Polygon File! \n";  break;
        case ERR_BADINITS : cout    << "Cannot Load Input File! \n";    break;
        case ERR_CREATEDIRS : cout  << "Cannot Create Directory! \n";   break;
    };
};



//extern double currentScaleFactor;
//extern double global_MaxSingleParticleOverlap;
//extern double global_MaxSystemOverlap;
extern enum class confinement_t { hard, soft, periodic, hexagonal  } BoundaryType;


using namespace std;

class GranularSimulation
{
private:
    
    //  Private Boolean
    //-------------------
    bool RunGranSim;
    
    
    
    //  Private String  (I/O)
    //-------------------------
    string outpath;
    
    
    
    // Private Integers
    //-------------------
    unsigned int simStep;                //  Define
    unsigned int runNumber;              //
    unsigned int numberOfShakes;         //
    unsigned int repackingSteps;
    unsigned int totalParticles;         //  Number of TOTAL Particles
    unsigned int totalPolygonSpecies;    //  Define
    unsigned int currentRepackingStep;
    unsigned int maxOptimizeShakeCycles; //


    
    
    
    //  Private Enum
    //----------------
    enum growth_type_t  { growing, nogrow, shrinking } GrowthType;
    confinement_t  gsBoxConfinement = confinement_t::hard;
    
    
    
    // Private Doubles
    //------------------
    double growthRate;
    double maxOverlap;
    double scaleFactor;
    double repackingScale;
    double shakingAmplitude;
    double maxCoveringAreaFraction;
    double orientationShakingAmplitude;


    
    //  Private Multidimensional Arrays
    //-----------------------------------
    // List of Global Polygons  - Not to be updated after Initialization
    //                          - size = # of Polygon Species
    vector<Particle>  globalParticles;
    
    // List of Particle-To-GlobalObject Maps - Can be updated
    //
    vector<int> objToShapeList;
    vector<int> objsPerShape;
    
    // Assign Boundary Conditions
    //-----------------------------
    vector<vector<double> > hardBoxShape;  // Container Box;
    
    
public:
    //  Constructor -------------
    GranularSimulation()
    {
       //  Default Values                   //   Default Interpretations
       //------------------                 //-----------------------------
        RunGranSim = true;                  // Start the Simulation
        
        runNumber = 1;                      // Default Ensemble Number = 1
        
        outpath="/conf/";                   // Set output path for files

        totalParticles = 1;                 // Total Particles in Simulation
        
        totalPolygonSpecies = 1 ;           // Number of Different Particles in Simulation
        
        maxOverlap = 0.00001;               // Percent of Current Inflated Area Allowed to Overlap
        
        shakingAmplitude = 1.0;             // Shaking = 1 Particle Diameter
        
        orientationShakingAmplitude = 1.0;  // Complete Angular Rotation
        
        numberOfShakes = 1000;              // Trial Monte Carlo Moves
        
        maxOptimizeShakeCycles = 10;        // Perform Additional Monte Carlo Moves
        
        repackingSteps=16;                  // Lower the Shaking Amplitude by repackingScale*( 0.8^step)
        
        currentRepackingStep=0;             // Start with No Repacking (full shaking amp)
        
        repackingScale=1.0;                 // Initialize Repacking Amplitude
        
        simStep = 1;                        // Start at the First SimStep
        
        growthRate = 0.00001;               // standardGrowthRate = 10^-5
        
        scaleFactor = simStep * growthRate; // Set the First Scale Factor for Particles
        
        maxCoveringAreaFraction= -1.0;      // Simulation is not stopped at certian area fraction
        
        globalParticles.clear(); objsPerShape.clear(); objToShapeList.clear(); hardBoxShape.clear();
        gsBoxConfinement = confinement_t::hard;  GrowthType = growing;
        
    };
    
    // Read Only Functions
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    //  Integers
    int getN()                                          {return totalParticles; };
    int getMaxOpt()                                     {return maxOptimizeShakeCycles; };
    int getNShakes()                                    {return numberOfShakes; };
    int getSimStep()                                    {return simStep; };
    int getRunNumber()                                  {return runNumber; };
    int getMaxRepacking()                               {return repackingSteps; };
    int getCurrentRepackingStep()                       {return currentRepackingStep; };
    int getNumberOfPolygonSpecies()                     {return totalPolygonSpecies; };
    int getParticlesPerSpecies(int in_whichPoly)        {return objsPerShape[in_whichPoly]; };
    int getParticleSpeciesID(int in_particleNumber)     {return objToShapeList[in_particleNumber]; };
    
    //  Doubles
    double getGrowthRate()      {return growthRate; };
    double getMaxOverlap()      {return maxOverlap; };
    double getShakingAmp()      {return shakingAmplitude; };
    double getRotationAmp()     {return orientationShakingAmplitude; };
    double getRepackScale()     {return repackingScale; };
    double getScaleFactor()     {scaleFactor = simStep*growthRate; return scaleFactor; };

    //  std::vector
    vector<vector<double> > getHardBox()                {return hardBoxShape; };
    vector<int> getParticleSpeciesIDlist()              {return objToShapeList; };
        
    //  ---------------

    
    
    
    //  ---------------
    vector<Particle> getListOfGlobalParticles(){return globalParticles; };
    Particle getGlobalParticle(int in_particleID){return globalParticles[in_particleID];};
    vector<vector<double> > getGlobalParticleVerticies(int in_polySpecifier) {return globalParticles[in_polySpecifier].getInitialShapeVertices();};



    string getOutpath(){return outpath;};
    
    
    
    
    
    
    
    // Set Update Functions
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    void setRunNumber(int in_runNumber)         {runNumber = in_runNumber; };
    void setOutputDirectory(string in_outdir)   {outpath = in_outdir; };

    void setHardBox(vector<vector<double> > in_hardBox) { hardBoxShape = in_hardBox;};
    
    void setGlobalParticles(vector<Particle> in_particles )
    {
        globalParticles.clear();
        globalParticles=in_particles;
        totalPolygonSpecies = in_particles.size();
        cout << "Setting Global Polys:  number " << globalParticles.size() <<"\n";
    };
    
    
    //  ---------------
    void setParticlePolydispersity(vector<int> in_polydisp)
    {
        objsPerShape = in_polydisp;
        totalParticles = accumulate(in_polydisp.begin(), in_polydisp.end(), 0);
        cout << "inside class - total particles " << totalParticles << "\n";
    };

    //  ---------------
    void setObjToShapeList(vector<int> in_objToShapeList) { objToShapeList = in_objToShapeList; };
    
    
    //  ---------------
    void setSimDetails(int in_maxOptimizeShakeCycles,
                       double in_growthRate,
                       double in_maxOverlap,
                       double in_shakingAmplitude,
                       double in_orientationShakingAmplitude,
                       int in_numberOfShakes
                       )
    {
        maxOptimizeShakeCycles = in_maxOptimizeShakeCycles;
        growthRate = in_growthRate;
        maxOverlap = in_maxOverlap;
        shakingAmplitude = in_shakingAmplitude;
        orientationShakingAmplitude = in_orientationShakingAmplitude;
        numberOfShakes = in_numberOfShakes;
    };
    
    
    void setInitialSimStep(int ss){simStep=ss;  scaleFactor= simStep * growthRate;  };
    
    void updateRepacking(){++currentRepackingStep; repackingScale=( pow(0.8, (currentRepackingStep))); };
    
    
    void updateSimStep() {
        switch(GrowthType){
            case growing : ++simStep;  break;
            case nogrow : break;
            case shrinking : --simStep; break;
        };
        scaleFactor = simStep * growthRate ;
    };
    
    void setMaxCoveringArea(double in_maxCoveringArea){maxCoveringAreaFraction = in_maxCoveringArea; };
    
    //---------
   // void setSimulationParticle(Particle in_Bead) {globalPolygonList.pushback(in_Bead);};
    
    

    // Write Simulation Details to Terminal
    // --------------------------------------------------------
    void PrintCurrentSimulationDetails( )
    {
        cout << "------------------------------------------------- \n";
        cout << " Total Number of Particles: " << totalParticles <<"\n";
        cout << " Number of Unique Particles: " << totalPolygonSpecies <<"\n";
        cout << " Sim Step: " << simStep <<"\n";
        cout << " Growth Rate: " << growthRate <<"\n";
        cout << " ScaleFactor: " << scaleFactor <<"\n";
        cout << "------------------------------------------------- \n";
    };
    
    
    // --------------------------------------------------------
    void getParticleSpecies_TotalParticles()
    {cout << "# of Polygon Species: " << totalPolygonSpecies
        << "\n # of Particles being Simulated: " << totalParticles << "\n"; };
};



// ------------------ Writes the passed Matrix to the Terminal
//  http://stackoverflow.com/questions/10368177/printing-contents-of-2d-vector
//
extern void PrintMatrixFunction( vector<vector<double> >  inMatrix )
{
    cout << " Printing Matrix with " << inMatrix.size() << " verticies \n";
    vector<vector<double> >::iterator it=inMatrix.begin(), end=inMatrix.end();
    while (it!=end)
    {
        vector<double>::iterator it1=it->begin(),end1=it->end();
        copy(it1,end1,ostream_iterator<double>(cout, " "));
        cout << endl;
        ++it;
    };
    
};
//==============================================//==============================================









//==============================================//==============================================
//==============================================//==============================================
//   Main Routine for Writing Data Files
//==============================================//==============================================
void WriteConfigurationPoints( CentroidList in_fullParticles, int in_simstep, double in_growthRate, int in_RunNumber, string in_outpath)
{
//  const char *dir_path[] = {in_outpath.c_str()};
//	boost::filesystem::path dir( in_outpath);
//	if(boost::filesystem::create_directory(dir)) {
//		std::cout << "Success" << "\n";
//	}
    
    //cout << "Printing File : ../3_data/" + in_outpath + "/conf_" << to_string(in_RunNumber) << ".dat ..." ;
    ofstream outfile;

    //  Prints File every
    //outfile.open ("environment/build_environment/Rcurrent/Rcurrent_" + to_string(in_RunNumber) + ".dat");
    //outfile.open("../3_data/"+ in_outpath +"/conf/conf_"+ to_string(in_RunNumber) +".dat");
    
    //  Use For Animation Frames
    outfile.open("../3_data/" + in_outpath +"/conf/conf_" + to_string(in_RunNumber) + "_"+ to_string(in_simstep) +".dat");
    

    //cout << "Begin Write Points (a) \n";
    outfile << in_simstep << "\t"
            << in_growthRate << "\t"
            << in_fullParticles.areaFraction << "\t"
            << in_fullParticles.cpuTime << "\n";
    
    for (int i=0; i < in_fullParticles.Centroids.size(); i++)
    {
        outfile <<  std::setprecision(16)
                <<  in_fullParticles.Centroids[i].X()       << "\t"
                <<  in_fullParticles.Centroids[i].Y()       << "\t"
                <<  in_fullParticles.Centroids[i].A()       << "\t"
                <<  in_fullParticles.Centroids[i].P() +1    << "\n";
    };
    
    outfile.close();
    // cout << "Done File \n";
};
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





//==============================================//==============================================
//==============================================//==============================================
extern void WritePolygons( vector<Particle>  in_globalParticles, string in_outpath)
{
    // cout << "Begin Write Points \n";
    int polys=0;
    for (polys=0; polys<in_globalParticles.size();polys++)
    {
        ofstream outfile;
        
        outfile.open ("../3_data/" + in_outpath +"/meta/simulation_obj/particle_" + to_string(polys+1) + ".dat");
        // cout << "Begin Write Points (a) \n";
        vector<double> polyX, polyY;
        polyX.clear(); polyY.clear();
        
        polyX = in_globalParticles[polys].getInitialShapeVerticesX();
        polyY = in_globalParticles[polys].getInitialShapeVerticesY();

        for (int i=0; i < polyX.size(); i++)
        {
            outfile << std::setprecision(16)
            <<  polyX[i] << "\t"
            <<  polyY[i] << "\n";
        };
        
        outfile.close();
    };
    // cout << "Done File \n";
};



//==============================================//==============================================
//==============================================//==============================================
extern void WriteBox( vector<vector<double> > in_boxArray, string in_outpath)
{
    ofstream outfile;
    outfile.open ("../3_data/" + in_outpath +"/meta/simulation_box/container.dat");
    int i=0;
    for (i=0; i < in_boxArray.size(); i++)
    {
        outfile << std::setprecision(16)
        <<  in_boxArray[i][0] << "\t"
        <<  in_boxArray[i][1] << "\n";
    };
};



//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// extra variables
extern SimulationAndShakingDetails_Struct setCurrentSimDetails(GranularSimulation GSIM)
{
    SimulationAndShakingDetails_Struct tempSimDetails;
    
    tempSimDetails.in_MaxShakes=GSIM.getNShakes();
    tempSimDetails.in_OptCycles=GSIM.getMaxOpt();
    tempSimDetails.in_ScaleFactor=GSIM.getScaleFactor();
    tempSimDetails.in_RepackScale=GSIM.getRepackScale();
    
    tempSimDetails.in_ShakeAmp=GSIM.getShakingAmp();
    tempSimDetails.in_RotAmp=GSIM.getRotationAmp();
    tempSimDetails.in_MaxOverlap=GSIM.getMaxOverlap();
    tempSimDetails.in_NumberOfPolygonSpecies=GSIM.getNumberOfPolygonSpecies();
    
    return tempSimDetails;
};

#endif /* gransim2_h */

