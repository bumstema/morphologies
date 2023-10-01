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

//  ---   Include OpenMP
#include <omp.h>

#include <cmath>
#include <math.h>
#include <time.h>
#include <random>     // Required: "random_shuffle()"
#include <vector>
#include <string>
#include <iostream>
#include <iterator>
#include <stdlib.h>
#include <iostream>
#include <algorithm>


//  Morphologies Scripts
#include "shake.h"
#include "particles.h"


//  Boost Random Numbers
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"

BOOST_GEOMETRY_REGISTER_C_ARRAY_CS(cs::cartesian)
typedef boost::mt19937 RandomNumberGenerator_Type;

using namespace std;
using namespace bgeom;




//===========================================================//===========================================================
polygon_2d MakeBoostPolygon( vector<vector <double> >  in_Data2d_To_Polygon ){
    polygon_2d poly;
    vector<point_2d> points;
    int npoints = in_Data2d_To_Polygon.size();
    points.resize(npoints);
    
    for ( int loop_i=0; loop_i < npoints;  loop_i++ )
    {
        points[ loop_i ] = point_2d(in_Data2d_To_Polygon[loop_i][0],in_Data2d_To_Polygon[loop_i][1]  );
    };
    
    assign_points(poly, points);
    correct(poly);
    return poly;
};
//extern double                    CalculateParticleOverlap( );
//===========================================================//===========================================================
double CalculateParticleOverlap( polygon_2d in_StationaryParticle, polygon_2d in_MovedParticle){
    double localOverlap{0.0};
    
    polygon_list intersectionPolygon;
    
    bgeom::intersection( in_StationaryParticle, in_MovedParticle, intersectionPolygon);
    
    for (polygon_list::const_iterator it = intersectionPolygon.begin(); it != intersectionPolygon.end(); ++it)
    {
        localOverlap += bgeom::area(*it);
    };
    return localOverlap;
};
//===========================================================//===========================================================


//===========================================================//===========================================================
//===========================================================//===========================================================
//===========================================================//===========================================================
//   Method for "Shake All Particles"
//===========================================================//===========================================================
CentroidList Shake(vector<Particle> in_Particles, CentroidList in_Points, SimulationAndShakingDetails_Struct Sim , Container in_Box )
{
    //  Steps to Calculate Acceptance Criteria
    //  0. Pre-Check - Inflate w/ No-Shake
    //  1. Check Hard Boundaries
    //  2. Check Soft Boundaries
    //  3. Check Particle-Particle Overlap
    //  4. Check Periodic Boundaries
    //      4.1.  Does Particle Satisfy Soft Boundary on AT LEAST ONE edge?
    //      4.2.  If YES -> Check Particle-Particle Overlap with Translated Particle
    //      4.3.  Save Centroid with (x,y,r) + (dx,dy,dr) + (ptX,ptY,ptZ)
    
    

    // Particles - Centroids and Acceptance
    // ------------------------------------------------
    CentroidList Points_now;            // Centroid Lists
    Points_now = in_Points;
    Points_now.acceptedQ = false;
    int nTotalParticles = in_Points.Centroids.size();
    vector<Particle> StericPotential = in_Particles;

    
    //  Initialize Particles as BoostPolygons
    //  -----------------------------------------------
    //  Vector of Polygons.
    vector<polygon_2d>  BoostPolygonParticles;
    BoostPolygonParticles.resize(nTotalParticles);

    
    
    vector<vector<double> > tempVertices;
    int polyID{0} ;
    double currentArea{0.0};
    double maxArea{0.0};
    double inflationConstant = Sim.in_ScaleFactor;
    
    
    //  Inflate Molecules to the Current Inflation Size
    // -------------------------------------------------
    for ( int i = 0; i < Sim.in_NumberOfPolygonSpecies; ++i )
    {
        StericPotential[ i ].setCurrentResizedVertices(inflationConstant);
    };
    
    
    //  Assign Vector of Polygons for Overlap Check
    // -------------------------------------------------
    for ( int objI = 0; objI < nTotalParticles; ++objI )
    {
        polyID = in_Points.Centroids[objI].P();
        BoostPolygonParticles[ objI ] = StericPotential[ polyID ].MakeBoostParticleAtCentroid( in_Points.Centroids[objI].getCentroid() );
        
        //  Calculate Area of the Biggest Particle Species - Used for Rejection.
        if(polyID > -1)
        {
            currentArea=bgeom::area(BoostPolygonParticles[objI]);
            if(currentArea > maxArea){maxArea=currentArea;};
        };
    };
    
    //  Allowable Overlap Area - Main Rejection Criteria
    //  ---------------------------------------------------------
    double acceptableOverlapArea{0.0};
    acceptableOverlapArea = ( maxArea * Sim.in_MaxOverlap );


    
    //  0. If Overlap < Rejection, Then Skip Shaking.
    //==========================================================================
    //  Pre-Check:  If overlap from growth is less than tolerence: Skip Shaking.
    // -------------------------------------------------------------------------
/*
    double totalOverlapArea{0.0};
    for (int firstID=1; firstID < nTotalParticles; ++firstID )
    {   for ( int secondID=0; secondID < firstID; ++secondID )
        {totalOverlapArea += CalculateParticleOverlap(BoostPolygonParticles[firstID], BoostPolygonParticles[secondID]);};
    };
*/
    
    //  Optimization for Speed
    double ddX{0.0};
    double ddY{0.0};
    double effectiveRadius{0.0};
    double totalOverlapArea{0.0};
    double bindingSphere{4.0 * pow(inflationConstant, 2.0)};

    
    //  Boundary Conditions -  Special Declaration for Boundary Condition Functions
    bool boundariesSatisfiedQ{false};
    bool acceptedMoveFound{true};
    
    for (int firstID=1; firstID < nTotalParticles; ++firstID )
    {
    for ( int secondID=0; secondID < firstID; ++secondID )
    {
        if(acceptedMoveFound == true){
        if(in_Box.PeriodicQ() == false){
        boundariesSatisfiedQ = in_Box.satifiedBoundaries(BoostPolygonParticles[firstID]) ;
        if(boundariesSatisfiedQ == false){ acceptedMoveFound = false; };
        };
        
        
        if( totalOverlapArea < acceptableOverlapArea)
           {
              // if( nTotalParticles > 100)
              // {
                   ddX = Points_now.Centroids[firstID].X() - Points_now.Centroids[secondID].X()  ;
                   ddY = Points_now.Centroids[firstID].Y() - Points_now.Centroids[secondID].Y()  ;
                   effectiveRadius = (pow(ddX, 2.0) + pow(ddY, 2.0)) ;
            //    cout << "effective Radius: " << effectiveRadius << " " <<  bindingSphere << "\n";
                   if( effectiveRadius < bindingSphere )
                   {
                       totalOverlapArea += CalculateParticleOverlap(BoostPolygonParticles[firstID], BoostPolygonParticles[secondID]);
                   };
               // };
           }else{ acceptedMoveFound = false; };
        };
    };
        
    };
    
    

    if ( acceptedMoveFound == true )
    {   Points_now.acceptedQ = true;
        Points_now.quickRejectQ = true;
        return Points_now;
    };
    //~~~~~~~~~~~~~~~~ Continue...
    //==========================================================================
    
    
    
    
    //  1. If Overlap > Rejection, Then Shake.
    //==========================================================================
    //  SHAKE :  If overlap from growth is greater than tolerence: Shake.
    //-------------------------------------------------------------------------
    
    
    
    //  Initialize Global Shaking Variables
    // -------------------------------------
    
    //  Current Simulation Details - SimStep and GrowthRate
    double rotateConstant  =  Sim.in_RotAmp   * Sim.in_RepackScale ;
    double shakingConstant =  Sim.in_ShakeAmp * Sim.in_RepackScale * inflationConstant;
    

    

    //  Initialize:  Loop Variables for Shaking
    // ------------------------------------------
    // Random Displacement Variables - {x,y,a}
    vector<double> dX(nTotalParticles, 0.0), dY(nTotalParticles, 0.0), dA(nTotalParticles, 0.0);
    
    
    
    // Initialize loop variables only once
    int acceptedMoveID{0};

    double  X_now, Y_now, A_now;
    vector<double> disp_Centroid(3,0.0);

    
    //  Randomize Particle Index for Randomly Permuted Shaking
    vector<int> RandomizedList;
    RandomizedList.clear();
    for (int i=0; i<nTotalParticles; ++i) RandomizedList.push_back(i);
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////
    // <----><----><----><----><----><----><----><----><----><----><----><----><----><----><---->
    //   :: Main Loop for Shaking ::
    //
    
    
//#pragma omp parallel for if(cycles>0) schedule(dynamic) \
    private(shakeN, disp_Obj, disp_Centroid,  totalOverlapArea) \
    shared(npart, X_now,Y_now,A_now,  acceptedMoveFound,  BoostPolygonParticles)

//#pragma omp parallel //schedule(dynamic)
//    {
//#pragma omp critical
//    {
//        cout << "Thread Number: " << omp_get_thread_num() <<"\n";
//        };
//    };
    
//#pragma omp parallel default(none) firstprivate(RandomizedList,nTotalParticles,Points_now,Sim) private(cycles, npart_,npart,dX,dY,dA,totalOverlapArea,acceptedMoveFound) shared(d_min,perm_min)
  //  {
    /*
     OpenMP Variables
     ------------
     Firstprivate:
     - Makes a copy of the distance matrix in each thread. Allows no overwritten data to matrix
     - Needs to be initialized and private within thread for correct initial starting point
     
     Private:
     - Loop Variables:  k,l  thread specific
     - needs to be private so no race-conditions when solving in parallel
     - seed: used to create new seed for threads
     
     Shared:
     - d_min:  Provides global minimum found in parallel
     - perm_min:  Global permutation of cities
     
     
     */
    
    
    //cout << //"step: " << GranSim.getSimStep() << " time: " << Points_now.cpuTime << " cycles: ";
   
    //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    for (int cycles = 0; cycles <= Sim.in_OptCycles; ++cycles)
    {
        //  Print Shuffle Cycle to Terminal
       // cout << cycles ;
        
        // Randomly Shuffled - Sample Centroids in Random Order
        random_shuffle ( RandomizedList.begin(), RandomizedList.end() );
        
        
        //  Monte Carlo Displacement Loops.
        // =======================================
        //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        for ( int npart_ = 0; npart_ < nTotalParticles; ++npart_)
        {
            //  Randomly Sampled Particle from Permutated Index List
            int npart{RandomizedList[npart_]};
            
            
            X_now = Points_now.Centroids[npart].X() ;
            Y_now = Points_now.Centroids[npart].Y() ;
            A_now = Points_now.Centroids[npart].A() ;
            polyID = Points_now.Centroids[npart].P();
            
            
            
            // Interactions
            //---------------
            // Populate Displacement Vectors
            dX = PopulateRandomDisplacements( Sim.in_MaxShakes, shakingConstant, npart  );
            dY = PopulateRandomDisplacements( Sim.in_MaxShakes, shakingConstant, npart_ );
            dA = PopulateRandomDisplacements( Sim.in_MaxShakes, rotateConstant , npart_ + npart + 2);

            
            // Apply Background Potential to Displacement Vectors
            /*
            double charge_X {0.0};
            if(polyID==0){ charge_X=(-0.3)*shakingConstant;}else{charge_X=(0.3)*shakingConstant;};
            for( int i_ =0; i_< dX.size(); ++i_){dX[i_] = dX[i_] + (charge_X) ; };
            */
            
            
            //  Set the First Trial = Last Accepted Move (i.e. move only particles that overlap)
            dX[Sim.in_MaxShakes-1]=0.0;  dY[Sim.in_MaxShakes-1]=0.0;  dA[Sim.in_MaxShakes-1]=0.0;
            //dX[0]=0.0;  dY[0]=0.0;  dA[0]=0.0;
            
            // Write out the Random Monte Carlo Trials
            /*
            if (npart == 0 ){
                ofstream outfileA;
                outfileA.open ("../3_data/data/jname/conf/RandomRolls.dat");
                for (int ig=0; ig < dX.size(); ++ig){ outfileA <<  dX[ig] << "\t"<< dY[ig]<<  "\t"<< dA[ig] <<"\n"; };
                outfileA.close();
            };
             */

            
            
            
            
            //  Reset loop variables
            //-----------------------
            acceptedMoveFound = false;
            acceptedMoveID = 0;
            totalOverlapArea = 0.0;

        //   SHAKE() = Results in Yes/No for "Do particle indicies need to be permuted?"
        //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        for (int shakeN = 0; shakeN <  Sim.in_MaxShakes; ++shakeN)
        {
            if ( acceptedMoveFound == false)
            {   //  First Break - If the Move is accepted
                //==============================================================

                // Displaced Object
                polygon_2d disp_Obj;
                vector<polygon_2d>  displacedBoostParticle;
                vector<vector<double> > pbcVectors;
                vector<vector<double> > pbcCentroid;
                
                //  Check Boundaries - Yes/No for "Is the move within boundary conditions?"
                // ------------------------------------------------------------------------
                displacedBoostParticle.clear();
                if (in_Box.PeriodicQ()==true)
                {
                    //  Fix centroid for pbc boundary (outside->homeotropic):
                    pbcVectors = in_Box.getPBCvectors();
                    for ( int i = 0; i < in_Box.nPBCvectors(); ++i  )
                    {
                        pbcCentroid.push_back({(X_now +  dX[shakeN] + pbcVectors[i][0]),(Y_now + dY[shakeN] + pbcVectors[i][1]),(A_now + dA[shakeN])});
                        displacedBoostParticle.push_back( StericPotential[ polyID ].MakeBoostParticleAtCentroid( pbcCentroid[i] ) );
                        if( in_Box.satifiedBoundaries(displacedBoostParticle[i]) == true){disp_Centroid = pbcCentroid[i]; boundariesSatisfiedQ = true; };
                    };
                    
                    
                }else
                {
                    disp_Centroid = {(X_now +  dX[shakeN]),(Y_now + dY[shakeN]),(A_now + dA[shakeN])};
                    displacedBoostParticle.push_back(  StericPotential[ polyID ].MakeBoostParticleAtCentroid( disp_Centroid ) );
                    boundariesSatisfiedQ = in_Box.satifiedBoundaries(displacedBoostParticle[0]);
                 //   if(boundariesSatisfiedQ==false){ cout << "Boundaries SatisfiedQ? " << boundariesSatisfiedQ <<  " - " <<  npart <<  " - " << shakeN << "\n";};
                };
                
                
                // Rejects the move if the boundaries are not satisfied
                if (boundariesSatisfiedQ == true)
                {

                       
                        //  Moved particle can be checked for overlap
                        //  Loop over particles and compare overlap <= allowable amount
                        //==============================================================
                        totalOverlapArea=0.0;
                        for ( int overlapObjID = 0; overlapObjID < nTotalParticles; ++overlapObjID)
                        {
                            
                            //  SKIP: If its NOT the same Particle   |OR|   Distance is less than 4*radius Distance away
                            if (   overlapObjID != npart ) //  Do not check particle overlap with itself
                            {
                                
                                //==============================================================
                                ///////////////////////////////////////////////////
                                //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//
                                for( int i = 0; i < displacedBoostParticle.size(); ++i )
                                {
                                    if(in_Box.PeriodicQ() == true)
                                    {
                                        ddX = Points_now.Centroids[overlapObjID].X() - pbcCentroid[i][0] ;
                                        ddY = Points_now.Centroids[overlapObjID].Y() - pbcCentroid[i][1] ;
                                        
                                    }
                                    else
                                    {
                                        ddX = Points_now.Centroids[overlapObjID].X() - disp_Centroid[0] ;
                                        ddY = Points_now.Centroids[overlapObjID].Y() - disp_Centroid[1] ;
                                    };
                                    effectiveRadius = (pow(ddX, 2.0) + pow(ddY, 2.0)) ;
                                //    cout << "effective Radius: " << effectiveRadius << " " <<  bindingSphere << "\n";
                                if( effectiveRadius < bindingSphere )
                                {
                                    if( totalOverlapArea < acceptableOverlapArea)
                                    {
                                        totalOverlapArea += CalculateParticleOverlap(BoostPolygonParticles[overlapObjID], displacedBoostParticle[i]);
                                    };
                                };
                               
                                };
                            };
                                //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//
                                //-----------------------------------------------------------
                            
                        };

                    
                        // If the overlap is bellow the acceptable amount...
                        if ( totalOverlapArea < acceptableOverlapArea  )
                        {
                            acceptedMoveFound = true;
                            acceptedMoveID = shakeN;
                        } ;//                                {Meets overlap criteria}
                        //-----------------------------------------------------------
          
                };  //  End "If" - {Skip particles if Particle is outside box}
                //-----------------------------------------------------------
            };    //cout << "...Done \n";  // End "if" - {Accepted Move == true}
            //-----------------------------------------------------------
        };      // End "For" Loop - {Randomly Displaced Centroid}
       //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        
    
        
        // If an ACCEPTED MOVE is FOUND - Update Current State
        // -----------------------------------------------------
        if (acceptedMoveFound == true)
        {
        // Saves the Centroid To List
            Points_now.Centroids[npart].setCentroid(disp_Centroid);
        // Save New BoostPolygon Back into List for Next Pass
            BoostPolygonParticles[npart] = StericPotential[ polyID ].MakeBoostParticleAtCentroid( disp_Centroid );
        };
            
        };
  
        
    //};//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
     // End of Allowable Randomized-Index-Permutaion Cycles
      
   
    
        
        
    
    
   
        
        
        
        
        
        
        
    //  --Cycles Loop------
        
    
    //   Test for Additional Optimization Cycles
    // -------------------------------------------
    
        // Recalculate Global Overlaps (consistancy through redundancy)
        //--------------------------------------------------------------
        
        
        totalOverlapArea=0.0;
        for (int firstID=1; firstID < nTotalParticles; ++firstID )
        {
            for ( int secondID=0; secondID < firstID; ++secondID )
            {
                        ddX = Points_now.Centroids[firstID].X() - Points_now.Centroids[secondID].X()  ;
                        ddY = Points_now.Centroids[firstID].Y() - Points_now.Centroids[secondID].Y()  ;
                        effectiveRadius = (pow(ddX, 2.0) + pow(ddY, 2.0)) ;
                        if( effectiveRadius < bindingSphere )
                        {
                            totalOverlapArea += CalculateParticleOverlap(BoostPolygonParticles[firstID], BoostPolygonParticles[secondID]);
                        };
            };
        };

        
        
        
        
        
        /*
        totalOverlapArea=0.0;
        for (int firstID=1; firstID < nTotalParticles; ++firstID )
        {   for ( int secondID=0; secondID < firstID; ++secondID )
            { totalOverlapArea += CalculateParticleOverlap(BoostPolygonParticles[firstID], BoostPolygonParticles[secondID]); };
        };
        */
        // cout << "Total overlap after shaking: " << totalOverlapArea << " is true=? " << (totalOverlapArea < Sim.in_MaxOverlap) << "\n";
        
        
        
        

        
        // If: Overlap is within Tolerance - Return New Positions
        //--------------------------------------------------------
        if (totalOverlapArea < (nTotalParticles * acceptableOverlapArea) )
        //if (Points_now.acceptedQ == true)
        {
            Points_now.acceptedQ = true;
            Points_now.quickRejectQ = false;
            Points_now.areaFraction = in_Box.CalculateAreaFraction(BoostPolygonParticles);
            return Points_now;
        }
        else
        // Else: Perform Next Opt Cycle or Exit if done.
        {
           if (cycles >= Sim.in_OptCycles)
           {
               //else check omp flag to quit
               
               
        //      cout << "- Rejected! ";
               in_Points.acceptedQ = false;
               in_Points.quickRejectQ = false;
               return in_Points;
           };
        };
    }; // End Optimization Cycles
    //
    //
    //
    //    End Of Monte Carlo Movement / Overlap Method Complete
    // <----><----><----><----><----><----><----><----><----><----><----><----><----><----><---->
    //////////////////////////////////////////////////////////////////////////////////////////////
    // If it reaches here before ending sooner, it has failed
    
    /* ===========================================================
     Critical Region changes the shared variables
     ============================================================== */
    //#pragma omp critical
    //{
    //    if (acceptedMoveFound==true ) {};
    // };
    
    //};  // End OpenMP Parallel
    
    
    
    
//End OpenMP Section
//    cout << " If you see this message, the simulation has failed.";
//    if (Points_now.acceptedQ == true)   { return Points_now; }
//    else                                { return in_Points ; };
    // ------------
    
    
    
//    cout << " If you see this message, the simulation has failed.";
//    in_Points.acceptedQ = false;
//    return in_Points;
};
//===========================================================//===========================================================

//===========================================================//===========================================================





//extern vector<double>            PopulateRandomDisplacements( );
//===========================================================//===========================================================
vector<double> PopulateRandomDisplacements(int in_numberOfRolls, double in_scale, int objNum)
{
    // typedef boost::variate_generator<base_generator_type&, distribution_type> gen_type
    int extraSeed = (int)(in_scale*100000+objNum);
    vector<double> randomNumberVector;
    RandomNumberGenerator_Type rng( time(0)+in_numberOfRolls*(extraSeed) );
    
    // boost::uniform_real<> uni_dist(0.0,in_scale);
    boost::normal_distribution<> nd(0.0, in_scale);
    boost::variate_generator<boost::mt19937&,  boost::normal_distribution<> >  randomize(rng, nd);
    
    randomNumberVector.resize(in_numberOfRolls);
    
    
    for ( int i=0; i < in_numberOfRolls; ++i)
    { randomNumberVector[i] = randomize(); };
    
    // cout << "Done Random Numbers \n" ;
    return randomNumberVector;
};
//===========================================================//===========================================================
//extern vector<vector<double> >   random_Initial_Positions( );
//===========================================================//===========================================================
vector<vector<double> > random_Initial_Positions( int num_Particles , vector<vector<double> > box)
{
    
    vector<double> random_XYR (3);
    vector<vector<double> > randomizedCentroid;
    //random_XYR.resize(3);
    randomizedCentroid.resize(0);
    
    vector<double> xbox (2,0.0) ;
    vector<double> ybox (2,0.0);
    double gMax,gMin;
    
    xbox={-10000,10000};
    ybox={-10000,10000};
    
    // cout << "xmax: " << xbox[0] << " xmin: " << xbox[1] << "\n";
    // cout << "ymax: " << ybox[0] << " ymin: " << ybox[1] << "\n";
    
    
    int temp;
    for ( temp=0; temp<box.size(); temp++)
    {if ( box[temp][0] > xbox[0]) xbox[0] = box[temp][0] ; // xmax
        if ( box[temp][0] < xbox[1]) xbox[1] = box[temp][0] ; // xmin
        if ( box[temp][1] > ybox[0]) ybox[0] = box[temp][1] ; // ymax
        if ( box[temp][1] < ybox[1]) ybox[1] = box[temp][1] ; // ymin
    };
    
    gMax= xbox[0];  if (  xbox[0] < ybox[0] ) {gMax=ybox[0];};
    gMin= xbox[1];  if (  xbox[1] > ybox[1] ) {gMin=ybox[1];};
    
    // cout << "xmax: " << xbox[0] << " xmin: " << xbox[1] << "\n";
    // cout << "ymax: " << ybox[0] << " ymin: " << ybox[1] << "\n";
    
    // cout << "gMax: " << gMax << " gMin: " << gMin << "\n";
    //  random_XYR = {0.1, 0.04, 02.2};
    
    // typedef boost::variate_generator<base_generator_type&, distribution_type> gen_type
    typedef boost::mt19937 RNGType;
    RNGType rng( time(0) );
    
    boost::uniform_real<> uni_dist(gMin,gMax);
    boost::variate_generator<boost::mt19937&,  boost::uniform_real<> >  randomize(rng, uni_dist);
    
    // random_XYR = {randomize(), randomize(), randomize()};
    
    // cout << "RandomizeTest: " << random_XYR[0] << " : " << random_XYR[1] << " : " << random_XYR[2]  << "\n";
    
    // BoostPolygon containerPolygon;
    polygon_2d containerPolygon;
    containerPolygon = MakeBoostPolygon(box);
    
    
    ///--------------------------------///
    //cout << "bpl Area a is: " << bpl::area(containerPolygon) << endl;
    cout << "bpl Area a is: " << area(containerPolygon) << endl;
    
    // BoostPoint testPoint(random_XYR[0] , random_XYR[1]);
    point_2d testPoint(random_XYR[0] , random_XYR[1]);
    
    temp=0;
    while( randomizedCentroid.size() < num_Particles) {
        
        random_XYR = {randomize(), randomize(), randomize()};
        // BoostPoint testPoint(random_XYR[0] , random_XYR[1]);
        point_2d testPoint(random_XYR[0] , random_XYR[1]);
        
        if ( within(testPoint, containerPolygon) )
        {
            randomizedCentroid.push_back(random_XYR);
            //  cout << "Centroid accepted : " << temp << " \t " << random_XYR[0] << " \t " << random_XYR[1] << endl;
            temp++;
        }
        
        
    } ;// i goes out of scope
    //testPoint = BoostPoint(randomize(), randomize());
    // std::cout << "within: " << (  boost::geometry::within(testPoint, containerPolygon) ? "yes" : "no") << std::endl;
    
    
    return randomizedCentroid;
};
//===========================================================//===========================================================

















