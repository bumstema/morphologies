/*
 *  M.o.r.p.h.o.l.o.g.i.e.s.
 *  ------------------------
 *
 *  Definitions on Particles for Simulations
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

//=================================
// include guard
#ifndef PARTICLES_h
#define PARTICLES_h

//=================================
// include dependancies
#include <stdio.h>
#include <vector>
#include <iostream>
#include <algorithm>
//#include <mkl.h>


//=================================
//  Boost Library Dependancies
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/adapted/c_array.hpp>
#include <boost/geometry/multi/geometries/multi_polygon.hpp>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
using namespace boost::geometry;
namespace bgeom = boost::geometry;

typedef model::d2::point_xy<double> point_2d;
typedef model::polygon<point_2d> polygon_2d;
typedef model::box<point_2d> box_2d;

using namespace std;
typedef vector<polygon_2d> polygon_list; // Used For "Intersection"



#define TwoPi  6.283185307179586


//-------------------------------------------------------------
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class Centroid
{
private:
    double  x_  {0.0};
    double  y_  {0.0};
    double  a_  {0.0};
    int     p_  {0};
    
    double  q_ {0.0};
    
    
public:
    //  Constructor:  Initialize to Zeros
    Centroid() = default;
    //  Constructor:  Initialize with Centroids
    Centroid(double ix, double iy, double ia): x_{ix}, y_{iy}, a_{ia} {} ;
    
    //  Constructor:  Initialize with Centroids and Polygon Specifier
    Centroid(double ix, double iy, double ia, int ip): x_{ix}, y_{iy}, a_{ia}, p_{ip} {} ;
    
    //  Constructor:  Initialize with Centroids
    Centroid(vector<double> c): x_{c[0]}, y_{c[1]}, a_{c[2]} {} ;
    
    //  Constructor:  Initialize with Centroids and Polygon Specifier
    Centroid(vector<double> c, int poly): x_{c[0]}, y_{c[1]}, a_{c[2]}, p_{poly} {} ;
    
    
    //------------------------------------------------//
    //  Read Functions
    
    double  X() {return x_;};
    double  Y() {return y_;};
    double  A() {return a_;};
    int     P() {return p_;};
    //----
    double  Q() {return q_;};

    void setX(double in_X)  {x_ = in_X;};
    void setY(double in_Y)  {y_ = in_Y;};
    void setA(double in_A)
    {
        //Reduce the angle to the interval from 0 to 2 pi
        while (in_A > TwoPi) in_A -= TwoPi;
        while (in_A < 0.0)  in_A += TwoPi;
        a_ = in_A;
    };
    void setP(int in_P)     {p_ = in_P;};
    void setQ(int in_Q)     {q_ = in_Q;};
    
    
    vector<double>  getCentroid()
    {
        vector<double> v;
        v={x_,y_,a_};
        return v;
    };
    
    void setCentroid(vector<double> v)
    {
        setX(v[0]);
        setY(v[1]);
        setA(v[2]);
    };
    
    void PrintCentroid()
    {
        cout << x_ << " " << y_ << " " << a_ << " " << p_ << "\n";
    };
};
//-----------------------------------------------------


//-------------------------------------------------------------
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
struct CentroidList
{
    vector<Centroid> Centroids;
 
    int simstep;
    int runNumber;
    
    bool acceptedQ;
    bool quickRejectQ;

    double cpuTime;
    double growthRate;
    double areaFraction;
};
//-----------------------------------------------------








//-------------------------------------------------------------
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
class Particle
{
private:
    //------------------------------------------------//
    //  Description of the Polygon Current SimStep
    //------------------------------------------------//
    int polyType_;
    int nVertex_;
    

    //------------------------------------------------//
    //  Each Object Can Have its own Simulation Details
    //------------------------------------------------//
    double particle_shakeAmp;
    double particle_rotAmp;
    double particle_overlap;
    //------------------------------------------------//
    
    
    //------------------------------------------------//
    //  Description of the Shape of Particles
    //------------------------------------------------//
    //  Notation: "Shape" defines the initial input polygon.
    vector<vector<double> > initialShapeVertices_;
    
    // Single X-Y : vectors used for Linear Algebra Routines
    vector<double>  initialShapeVerticesX_ ;
    vector<double>  initialShapeVerticesY_ ;
    
    // Single X-Y : vectors used for Linear Algebra Routines
    vector<double>  currentResizedVerticesX_ ;
    vector<double>  currentResizedVerticesY_ ;
    
    
public:
    //  Particle Constructor
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    Particle() {     };
    
    //  Read Functions
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    int getPolyType()       {return polyType_;};
    int getNumVertices()    {return nVertex_;} ;
    
    vector<double>          getInitialShapeVerticesX()      {return initialShapeVerticesX_;     };
    vector<double>          getInitialShapeVerticesY()      {return initialShapeVerticesY_;     };
    vector<vector<double> > getInitialShapeVertices()       {return initialShapeVertices_;      };
    
    vector<double>          getCurrentResizedVerticesX()    {return currentResizedVerticesX_;   };
    vector<double>          getCurrentResizedVerticesY()    {return currentResizedVerticesY_;   };
    
    
    //  Update Functions
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    //=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    void setPolyType(int in_polyType)   {polyType_ = in_polyType;};
    
    void Initialize(int in_polyType, vector<vector<double> > in_VertexList)
    {
        setInitialShapeVertices(in_VertexList);
        polyType_   =   in_polyType;
    };
 
    

    //  Routines to Set Initial Input Polygons and Current Simulation Polygon
    void setInitialShapeVertices( vector<vector<double> > in_VertexList )
    {
        nVertex_ = in_VertexList.size();
        initialShapeVertices_ = in_VertexList;
        
        initialShapeVerticesX_.resize(nVertex_);
        initialShapeVerticesY_.resize(nVertex_);
        
        currentResizedVerticesX_.resize(nVertex_);
        currentResizedVerticesY_.resize(nVertex_);
        
        for (int i=0;i<nVertex_; ++i)
        {
            initialShapeVerticesX_[i]   =   in_VertexList[i][0];
            initialShapeVerticesY_[i]   =   in_VertexList[i][1];
        };
    };
    
    //-----------------------------------------------------------------   
    void setCurrentResizedVertices( double in_rescale )
    {
        vector<double>  tempX, tempY;
        tempX = getInitialShapeVerticesX();
        tempY = getInitialShapeVerticesY();
        
        for (int i=0; i<tempX.size(); ++i)
        {
            currentResizedVerticesX_[i]=  tempX[i] * in_rescale ;
            currentResizedVerticesY_[i]=  tempY[i] * in_rescale ;
        };
    };

    
    ////////////////////////////////////////////////////////////////////////////////
    //  *    Translate + (Rotate the Polygon[ [Scale the Polygon[polygon]]])  *  //
    //////////////////////////////////////////////////////////////////////////////
    //-----------------------------------------------------------------
    polygon_2d MakeBoostParticleAtCentroid( vector<double> in_Centroid )
    {
       // vector<double> shape_verticesX ;
      //  vector<double> shape_verticesY ;
      //  shape_verticesX = getCurrentResizedVerticesX();
       // shape_verticesY = getCurrentResizedVerticesY();
        
        // vector<vector<double> > return_Vertices ;
        // return_Vertices = getInitialShapeVertices();
        
        double m_dAngle = in_Centroid[2];
        double rotate_x = cos(m_dAngle) ;
        double rotate_y =  -sin(m_dAngle);
        
        
        double polyX, polyY;
        polygon_2d poly;
        vector<point_2d> points;
        points.resize(nVertex_);

        for (int i=0; i<nVertex_; ++i)
        {
            polyX = in_Centroid[0] + ((rotate_x * currentResizedVerticesX_[i]) + (rotate_y * currentResizedVerticesY_[i]));
            polyY = in_Centroid[1] + ((rotate_x * currentResizedVerticesY_[i]) - (rotate_y * currentResizedVerticesX_[i]));
        
            points[ i ] = point_2d( polyX, polyY );
        };
        
        assign_points(poly, points);
        correct(poly);
        
        return poly;
    };

    
   //  ---------------
    
    
    //  Misc Functions
    // ------------------ Writes the passed Matrix to the Terminal
    //  http://stackoverflow.com/questions/10368177/printing-contents-of-2d-vector
    //  ------------------------------------------------------------------------
    void PrintMatrix( vector<vector <double> > inMatrix)
    {
        cout << " Printing Matrix : \n";
        vector<vector<double> >::iterator it=inMatrix.begin(), end=inMatrix.end();
        while (it!=end)
        {
            vector<double>::iterator it1=it->begin(),end1=it->end();
            copy(it1,end1,ostream_iterator<double>(cout, " "));
            cout << endl;
            ++it;
        };
    };
    
    
    void PrintParticle()
    {
        cout << " Printing Particle : \n";
        cout << " polyType :" << polyType_ << " nVertex: " << nVertex_ << " \n";
        cout << " Initial Polygon Vertices : \n";
        PrintMatrix(initialShapeVertices_);
    };
    
};



//===================================
struct List_of_Particles_struct
{
     vector< Particle > simParticle;
} ;



#endif /* particles_h */

