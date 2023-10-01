/************************************************
 *  M.o.r.p.h.o.l.o.g.i.e.s.
 *  ------------------------
 *
 *  Definitions on Particles for Simulations
 ***********************************************/
//=================================
// include guard
#ifndef INITIALIZATION_H
#define INITIALIZATION_H

//=================================
// include dependancies
#include <vector>
#include <math.h>

//=================================
// include Morphologies Script
#include "particles.h"


using namespace std;

//===================================
struct input_polygons_struct {
    vector<int> read_nPolySpecies;
    vector<Particle> read_inputPolygons;
    vector<int> read_polySpecifList;
    bool ErrorLoadingQ;
} ;





//=================
class read_input {
    
public:
    
    int read(int argc, char* argv[]);
    
    //-------------------------------------------------
    //   GranSim 1 Variables
    //-------------------------------------------------
        
    int countlines(char* filename);
    int LoadPolygonFile(char* filename);
    int loadPolygons(char* VerticesFile);
    int loadPolygonVertexFromPolyLink(char* filename );

    
    bool Exchange;
    bool Optimize;
    bool Shake;
    bool RelativeShaking;
    bool VF;
    bool PrintInitialFinalOnly;

    
    char boxFile[512];
    char OutputDir[512];
    char ResumeFile[512];
    char outputPath[512];
    char VerticesFile[512];
    char BoundaryFile[512];
    char ParticlesFile[512];
    char SimulationName[512];

    

    enum optimization_t { none, boundingsphere, pbc } SpeedOptimization;
    enum confinement_t { hard, periodic, soft, hexagonal } BoxConfinement;

    
    
    unsigned int SimulationStep;
    unsigned int NoIterations;
    unsigned int NumberOfShakes;
    unsigned int MaxOptimizeShakeCycles;
    unsigned int MaxPolydisperity;
    unsigned int NumberOfBeads;
    unsigned int PolygonVertices;
    
    
    double GrowthRate;
    double MaxOverlap;
    double ShakingAmplitude;
    double OrientationShakingAmplitude;
    double NormPolygonArea;
    
    double maxCoveringAreaFraction;
    
    double* PolygonXCoordinates;
    double* PolygonYCoordinates;
    
    
};














extern int countlines(char *filename);
//extern int loadPolygonVertex(char *filename);
extern std::vector <std::vector<double> > loadPolygonVertex(char *filename);


//extern int loadPolygonVertexFromPolyLink(char *filename );
input_polygons_struct loadPolygonVertexFromPolyLink(char *filename );




//extern std::vector<Particle> ext_inputPolygons;
//extern int ext_NumberOfPolygonSpecies;
//extern std::vector<int> ext_NumberOfParticlesPerSpecies;

#endif