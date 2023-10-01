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
 * Last modified:  2017.02.08
 *
 ***********************************************/

//-------//-------//-------//-------//-------//-------//-------//-------//-------
//-------//-------//-------//-------//-------//-------//-------//-------//-------
//  (Main Loop - A.) Run Simulations
//-------//-------//-------//-------//-------//-------//-------//-------//-------
// A.  Prechecking
//    A.1. - Pre-Update the Simulation Step
//    A.2. - Calculate Global Overlap from Inflation
//    A.2.bool   If Overlap ABOVE Tolerance - Shake the Particles
//    A.2.bool   If Overlap BELOW Tolerance - Skip Shaking and Grow Again
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------



//-------//-------//-------//-------//-------//-------//-------//-------//-------
//-------//-------//-------//-------//-------//-------//-------//-------//-------
//  (Main Loop - B.) Shake
//-------//-------//-------//-------//-------//-------//-------//-------//-------
// B.  Pass to Shake.h
//    B.1. - Pre-Update the Simulation Step
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------





//
//  ---   Include OpenMP
//#include <omp.h>
//  ---------------------

#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
//#include <time.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
//#include <random>
//#include <locale>
#include <string.h>

#include <ctime>  // used to time the main loop


// Main Method

#include "morphologies.h"


// Load polygon Files - Simulation Environment
#include "initialization.h"
// particle class - stores information of Particles
#include "particles.h"
#include "shake.h"
//#include "boundary.h"

using namespace std;



int main(int argc, char **argv)
{
    
// ----------------------------------------------
// ----------------------------------------------
//    Starting M.o.r.p.h.o.l.o.g.i.e.s. Main Code
    
//  - Useful Definitions:
//      GranSim.globalParticles   <--- List of Input Polygons
// ----------------------------------------------
    cout << "Starting M.o.r.p.h.o.l.o.g.i.e.s. \n";
    
    // Set Simulation Parameters
    GranularSimulation GranSim;
    
    //  SimulationStep Test - Working
    //cout<< " Testing SimStep: " << GranSim.getSimStep();
    
    
    
    // Read Input File
    read_input input;
    int error = input.read(argc, argv);
    if (error) return error;
    
    //  If Error Loading - Stop Simulation
    bool RunGranSim ;
    if (error == 0) {RunGranSim = true;} else {RunGranSim = false;};
    
    
    
    GranSim.setSimDetails(input.MaxOptimizeShakeCycles,
                          input.GrowthRate,
                          input.MaxOverlap,
                          input.ShakingAmplitude,
                          input.OrientationShakingAmplitude,
                          input.NumberOfShakes
                          );

    
    
    //if (argv[1][0] == '-')
    GranSim.setRunNumber(stoi(argv[1]));
    cout << "Setting to RunNumber: " << GranSim.getRunNumber() <<"\n";
    GranSim.setOutputDirectory(input.outputPath);
    cout << "Writing Output Conf file to: " << GranSim.getOutpath() <<"\n";
    
    // Initialize
    //+++++++++++++++Read Polygon+++++++++++++++++

    cout << "~~~~~~~~~~~~~~~ Testing Global Particles ~~~~~~~~~~~~~~~~~~ \n";
    // Get initial Polygons from File
    input_polygons_struct initial_polydispersity;
    initial_polydispersity = loadPolygonVertexFromPolyLink(input.ParticlesFile);
    
    GranSim.setParticlePolydispersity(initial_polydispersity.read_nPolySpecies);
    GranSim.setGlobalParticles(initial_polydispersity.read_inputPolygons);
    GranSim.setObjToShapeList(initial_polydispersity.read_polySpecifList);
    if (initial_polydispersity.ErrorLoadingQ)
    {
        PrintError(ERR_BADPF);
        return 0;
    };
    
    
    

    

    
    
    
    
    
    //   Initialize Simulation Particle List
     std::cout << "~~~~~~~~~~~~~~~ Assigning Particle Class ~~~~~~~~~~~~~~~~~~ \n";

    List_of_Particles_struct allParticles;
    //vector<Particle>  globalPolygonList;
    //globalPolygonList = GranSim.getListOfGlobalParticles();
    // allParticles.simParticle.clear();
    // allParticles ;


    for (int i=0; i<GranSim.getNumberOfPolygonSpecies(); ++i)
    {
        allParticles.simParticle.push_back(Particle());
        allParticles.simParticle[i].Initialize( i, GranSim.getGlobalParticleVerticies(i));
    };
    //cout << "~~~~~~~~~~~~~~~ Testing Simulation Polygons ~~~~~~~~~~~~~~~~~~ \n";
    allParticles.simParticle[0].PrintParticle();
    
    

    // cout << "~~~~~~~~~~~~~~~ Testing Container ~~~~~~~~~~~~~~~~~~ \n";
    cout << "~~~~Box Confinement - Class Container \n";
    cout << "~~~~Box Confinement type: " << input.BoundaryFile <<  " \n";
    
    // Load box into Memory
    vector<vector<double> > boxPolygon;
    boxPolygon = loadPolygonVertex( input.BoundaryFile );
    GranSim.setHardBox( boxPolygon );
    
    // Assign Proper Edges to Box
    Container Box;
    Box.setBoundaries( boxPolygon );
    Box.getBoxParticle().PrintParticle();
    if ( input.BoxConfinement == input.confinement_t::hard )    { Box.setAs_Hard(); };
    if ( input.BoxConfinement == input.confinement_t::soft )    { Box.setAs_Soft(); };
    if ( input.BoxConfinement == input.confinement_t::periodic ){
        cout << "Boxtype: Periodic. \n" ;
        Box.setAs_Periodic();
        GranSim.setHardBox( Box.getContainerVertices() );
    };
    if ( input.BoxConfinement == input.confinement_t::hexagonal ){
        cout << "Boxtype: Periodic - Hexagonal. \n" ;
        Box.setAs_PeriodicHexagonal();
        GranSim.setHardBox( Box.getContainerVertices() );
    };

    
    
    
    bool a = false;
    a = Box.satifiedBoundaries(MakeBoostPolygon(allParticles.simParticle[0].getInitialShapeVertices()));
    cout << "Testing if box satisfied boundaries  Boundary.h : " ;
    if(a){cout <<"true";}else {cout <<"false";}
    cout <<" ... complete... \n" ;
    

    
    
    CentroidList Points;
    CentroidList PreviousPoints;

    Points.Centroids.clear();
    PreviousPoints.Centroids.clear();
    
    //   Initialize Random Position list
    cout << "~~~~~~~~~~~~~~~ Initialize Random Position list ~~~~~~~~~~~~~~~~~~ \n";
    vector<vector<double> > rand_init_centroids;
    rand_init_centroids = random_Initial_Positions( GranSim.getN(), GranSim.getHardBox() );
    //cout << "Returned Centroid: "<< rand_init_centroids[0] << "\n";
    
    for (int i=0; i<GranSim.getN(); i++)
    {
      //  allParticles.simParticle[i].updateCentroidV(rand_init_centroids[i]);
        Points.Centroids.push_back(Centroid(rand_init_centroids[i]));
        Points.Centroids[i].setP( GranSim.getParticleSpeciesID(i) );
    //    Points.Centroids[i].PrintCentroid();
    };
    

    Points.cpuTime      = 0.0 ;
    Points.areaFraction = 0.0 ;
    
    Points.acceptedQ    = false;
    Points.quickRejectQ = false;

    GranSim.setInitialSimStep(1);
    Points.simstep      = GranSim.getSimStep();
    Points.runNumber    = GranSim.getRunNumber();
    Points.growthRate   = GranSim.getGrowthRate();

    
    //cout << "Points.PolyIndex.size = "<< Points.size() << "  " << Points.PolyIndex[0]  << "\n" ;
    cout << "Points.Centroids.size = "<< Points.Centroids.size() << "  "<< Points.Centroids[0].X() <<"\n" ;
    
    
    

    //cout << "~~~~~~~~~~~~~~~ WRITE SIMULATION ENVIRONMENT ~~~~~~~~~~~~~~~~~~ \n";
    if(GranSim.getRunNumber() == 1){
    WritePolygons(GranSim.getListOfGlobalParticles(), GranSim.getOutpath() );
    WriteBox( GranSim.getHardBox() , GranSim.getOutpath());

    };
    //cout << "~~~~~~~~~~~~~~~ WRITE SIMULATION ENVIRONMENT ~~~~~~~~~~~~~~~~~~ \n";
    
    
    
    



    
    double currentOverlappingArea = 0.0;

    
    
    // Variable that Breaks while Loop
    int finishSimulation = 0;
    double twiceMaxGrowthSteps = (2./GranSim.getGrowthRate());
    

    
    SimulationAndShakingDetails_Struct currentSimDetails;
    currentSimDetails   = setCurrentSimDetails(GranSim);
    


    
    

    
    // Writing first positional configurations
    WriteConfigurationPoints(Points, GranSim.getSimStep(), GranSim.getGrowthRate(), GranSim.getRunNumber(), GranSim.getOutpath() );

    
    

    // Start the Clock
    clock_t start = clock();

    
    while (RunGranSim != false)
    {
        
       // Points.cpuTime= ((double)clock() - start) / CLOCKS_PER_SEC ;
       // cout << "simstep: " << GranSim.getSimStep() << " time: " << Points.cpuTime << " ";
        PreviousPoints = Points;
        
        
        //   Main Routine:
        // --------------------------------------------------------------------------------------
        //                                      Shake();
        //
            Points = Shake(allParticles.simParticle, PreviousPoints, currentSimDetails, Box );
        // --------------------------------------------------------------------------------------

        
        
        
        
        // Check Stopping Conditions
        //---------------------------------------------------------------------------
        //---------------------------------------------------------------------------
        // (1. Number of Monte Carlo Trials - All Moves Rejected )
        // (2. Repacking Steps Completed )
        // (3. No Inflation - Phase Space Explored )
        // (4. Breaking While Loop -  )
        
        
        
        // (Stopping Condition #1. - All MC Moves Rejected)
        //-------------------------------------------------
        if (Points.acceptedQ)
        {
            // Write Files - then Update SimStep Last
            // WriteConfigurationPoints(Points, GranSim.getSimStep(), GranSim.getGrowthRate(), GranSim.getRunNumber(), GranSim.getOutpath() );
            
            
            //   Main Routine:
            // --------------------------------------------------------
            //      Inflate();
            //
                    GranSim.updateSimStep();
                    Points.simstep      =   GranSim.getSimStep();
                    currentSimDetails   =   setCurrentSimDetails(GranSim);
            
            if (Points.quickRejectQ == false)
            {
                Points.cpuTime= ((double)clock() - start) / CLOCKS_PER_SEC ;
                cout << "step: "<< Points.simstep << " area: "<< Points.areaFraction << " time: " << Points.cpuTime << "\n";
                WriteConfigurationPoints(Points, GranSim.getSimStep(), GranSim.getGrowthRate(), GranSim.getRunNumber(), GranSim.getOutpath() );

            };
            // --------------------------------------------------------
        }
        else
        {
            Points = PreviousPoints ;
            if (GranSim.getCurrentRepackingStep() < GranSim.getMaxRepacking() )
            //  If NOT last Repacking Step - Lower Shaking Amplitude // Do not update simstep
            {
                GranSim.updateRepacking();
                cout << "  Tg: " << GranSim.getCurrentRepackingStep() << " scale: " << GranSim.getRepackScale() << " reshake. \n";
                Points.cpuTime= ((double)clock() - start) / CLOCKS_PER_SEC ;
                WriteConfigurationPoints(Points, GranSim.getSimStep(), GranSim.getGrowthRate(), GranSim.getRunNumber(), GranSim.getOutpath() );
                currentSimDetails=setCurrentSimDetails(GranSim);
            }
            else
            //  Else
            {
                cout << "\n____Simulation Complete! ";
                // Stop Simulation
                RunGranSim = false;  break;
            };

        };

        
        // (Stopping Condition #4. - Too Many Steps - Breaks the "While" Loop)
        //-------------------------------------------------
        //
        finishSimulation++;
        if (finishSimulation >= twiceMaxGrowthSteps) { RunGranSim = false; break;};
    }

    Points.cpuTime= ((double)clock() - start) / CLOCKS_PER_SEC ;
    WriteConfigurationPoints(Points, GranSim.getSimStep(), GranSim.getGrowthRate(), GranSim.getRunNumber(), GranSim.getOutpath());
    
    cout << " Time Taken: " << Points.cpuTime << "\n";
    return 0;
}
///======================================================================================
///======================================================================================
///======================================================================================
