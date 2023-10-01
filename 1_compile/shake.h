/*
 *  M.o.r.p.h.o.l.o.g.i.e.s.
 *  ------------------------
 *  - Random Processes
 *  - Vector/Matrix Operactions
 *  - Polygon Overlap
 ************************************************
 *
 * McMaster Engineering Physics
 * Written by: Matt Bumstead
 *
 * File generated: 2015.10.20
 * Last modified:  2017.02.08
 *
 ***********************************************/
//=================================
// include guard
#ifndef SHAKE_h
#define SHAKE_h

//=================================
// include dependancies
#include <vector>
#include <math.h>
#include <iostream>
#include <time.h>
#include <random>
#include <fstream>

//=================================
//  Boost Library Dependancies
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/adapted/c_array.hpp>
#include <boost/geometry/multi/geometries/multi_polygon.hpp>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>

//  Boost Random Numbers
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include <list>
#include <cassert>
#include <deque>

#include "particles.h"



using namespace boost::geometry;
namespace bgeom = boost::geometry;

typedef model::d2::point_xy<double> point_2d;
typedef model::polygon<point_2d> polygon_2d;
typedef model::box<point_2d> box_2d;
typedef vector<polygon_2d> polygon_list; // Used For "Intersection"

using namespace std;




//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
// repeated variables
struct SimulationAndShakingDetails_Struct
{
    int in_MaxShakes;
    int in_OptCycles;
    double in_ScaleFactor;
    double in_RepackScale;
    double in_ShakeAmp;
    double in_RotAmp;
    double in_MaxOverlap;
    int in_NumberOfPolygonSpecies;
};





//extern polygon_2d                MakeBoostPolygon( );
//===========================================================//===========================================================
extern polygon_2d MakeBoostPolygon( vector<vector <double> >  in_Data2d_To_Polygon );
//extern double                    CalculateParticleOverlap( );
//===========================================================//===========================================================
extern double CalculateParticleOverlap( polygon_2d in_StationaryParticle, polygon_2d in_MovedParticle);
//===========================================================//===========================================================



//===============================================================
//  Container Class
//===============================================================
class Container
{
private:
    // Deals with the geometry of the container
    //-----------------------------------------
    vector<vector<double> > containerVertices_;
    Particle boundaryAsParticle_;
    polygon_2d boundaryAsBoost_;
    double containerArea_;
    
    // Types of Confinement
    //----------------
    bool soft_                  {false};
    bool hard_                  {false};
    bool pbcHex_                {false};
    bool pbcSqr_                {false};
    bool periodic_              {false};
    bool boundarySatisfiedQ_    {false};
    
    
    // Built in PBC Vectors
    //-------------------------------------------------------------------
    vector<vector<double> > boxSquare_ {{{0,0},{0,1},{1,1},{1,0},{0,0}}};
    vector<vector<double> > pbcTranslationalBasisSquare_ {{{0,0},{0,1},{1,0},{0,-1.},{-1.,0},{1,1},{-1.,-1},{1.,-1.},{-1.,1.}}};
    
    vector<vector<double> > boxHexagon_ {{{1., 0.5}, {0.75, 0.933013}, {0.25, 0.933013}, {0., 0.5}, {0.25, 0.0669873}, {0.75, 0.0669873}, {1., 0.5}}};
    vector<vector<double> > pbcTranslationalBasisHexagon_ {{{0.0, 0.0}, {0.0, 0.866025}, {-0.75, 0.433013}, {-0.75, -0.433013}, {0.0, -0.866025}, {0.75, -0.433013}, {0.75, 0.433013}}};
    
    vector<vector<double> > periodicTranslation_;
    
public:
    // constructor
    Container(){};
    
    bool HardQ()    { return hard_; };
    bool SoftQ()    { return soft_; };
    bool PeriodicQ(){ return periodic_ ; };
    bool PeriodicSqrQ() { return pbcSqr_ ; };
    bool PeriodicHexQ() { return pbcHex_ ; };

    void setAs_Hard()               { hard_ = true; };
    void setAs_Soft()               { soft_ = true; };
    void setAs_Periodic()           { periodic_ = true; pbcSqr_ = true; periodicTranslation_ = pbcTranslationalBasisSquare_ ; setBoundaries(boxSquare_ ); };
    void setAs_PeriodicHexagonal()  { periodic_ = true; pbcHex_ = true; periodicTranslation_ = pbcTranslationalBasisHexagon_; setBoundaries(boxHexagon_); };


    
    void setBoundaries(vector<vector<double> > in_BoxVertices)
    {
        containerVertices_ = in_BoxVertices;
        boundaryAsParticle_.Initialize( 0, in_BoxVertices);
        boundaryAsBoost_ = MakeBoostPolygon(in_BoxVertices);
        containerArea_ = bgeom::area(boundaryAsBoost_);
    };

    
    
    
    
    //  Container: Return Functions
    //------------------------------
    double getBoxArea()         { return containerArea_; };
    Particle getBoxParticle()   { return boundaryAsParticle_; };
    vector<vector<double> > getContainerVertices(){ return containerVertices_ ; };
    
    

    
    
    
    bool satifiedBoundaries(polygon_2d in_BoostParticle){
        bool SatisfiedQ {false};
        if( HardQ() == true )       { SatisfiedQ = HardBoundary(in_BoostParticle);          };
        if( SoftQ() == true )       { SatisfiedQ = HomeotropicBoundary(in_BoostParticle);   };
        if( PeriodicQ() == true )   { SatisfiedQ = HomeotropicBoundary(in_BoostParticle);   };
        return SatisfiedQ;
    };
    
    
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //          ||      ||      ||      ||      ||      ||
    //          ||      ||      ||      ||      ||      ||
    //===============================================================
    //   (hard)  Hard Boundary
    //===============================================================
    //
    //  Q.  Is the Particle Completely Contained Inside?
    //      .{true}  Return True  (Possible Accepted Move)
    //      .{false} Return False (Rejected Move)
    //
    bool HardBoundary(polygon_2d in_BoostParticle){
        
        boundarySatisfiedQ_ = bgeom::within(in_BoostParticle, boundaryAsBoost_);

    return boundarySatisfiedQ_;
    };
    //
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //          ||      ||      ||      ||      ||      ||
    //          ||      ||      ||      ||      ||      ||
    //===============================================================
    //   (soft)  Homeotropic Boundary
    //==============================================================
    //
    //  Q.  Is the Centroid Still Inside Container?
    //      .{true}  Return True  (Possible Accepted Move)
    //      .{false} Return False (Rejected Move)
    //
    bool HomeotropicBoundary(polygon_2d in_BoostParticle) {

        point_2d p;
        bgeom::centroid(in_BoostParticle, p);
        boundarySatisfiedQ_ = bgeom::within( p, boundaryAsBoost_);
            
    return boundarySatisfiedQ_;
    };
    //
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //          ||      ||      ||      ||      ||      ||
    //          ||      ||      ||      ||      ||      ||
    //===============================================================
    //   (pbc)  Periodic Boundary
    //==============================================================
    //
    int nPBCvectors()                       { return periodicTranslation_.size(); };
    vector<vector<double> > getPBCvectors() { return periodicTranslation_; };
    //
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 
    //          ||      ||      ||      ||      ||      ||
    //===============================================================
    double CalculateAreaFraction( polygon_list in_BoostParticles )
    {
        double insideArea{0.0};
        double areaFraction{0.0};
        
        
        if (soft_ != true )//  Soft Boundaries need to be clipped to be inside the box.
            {
                for (polygon_list::const_iterator it = in_BoostParticles.begin(); it != in_BoostParticles.end(); ++it)
                {   insideArea += bgeom::area(*it);};
            } else             //  All other Boundaries do not require clipping.
            {
                for (polygon_list::const_iterator it = in_BoostParticles.begin(); it != in_BoostParticles.end(); ++it)
                {   insideArea += CalculateParticleOverlap(*it, boundaryAsBoost_  ); };
                };
        
        areaFraction  = insideArea / containerArea_ ;
        return areaFraction;
    };
    //===============================================================
    //          ||      ||      ||      ||      ||      ||

};

extern CentroidList Shake(vector<Particle> in_Particles, CentroidList in_Points, SimulationAndShakingDetails_Struct Sim , Container in_Box );


//extern vector<double>            PopulateRandomDisplacements( );
//===========================================================//===========================================================
extern vector<double> PopulateRandomDisplacements(int in_numberOfRolls, double in_scale, int objNum);

//===========================================================//===========================================================
//extern vector<vector<double> >   random_Initial_Positions( );
//===========================================================//===========================================================
extern vector<vector<double> > random_Initial_Positions( int num_Particles , vector<vector<double> > box);

//===========================================================//===========================================================

#endif


