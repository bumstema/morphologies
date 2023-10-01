/*
 *  M.o.r.p.h.o.l.o.g.i.e.s.
 *  ------------------------
 *
 *  Initialize Particles with User Inputs
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

#include <vector>
#include <cmath>
#include <stdlib.h>
#include <iostream>

#include <math.h>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <cstring>

#include "particles.h"
#include "initialization.h"



//using namespace std;
// =======================================================================
// =======================================================================


// Source File for input
//      SPHERES
//================================================================
int read_input::read(int argc, char * argv[])
{
    int error = 0;
    
    //------------------------------------------------------
    //  Load GranSim 1 Data File
    //------------------------------------------------------
    
    
    FILE* pFile = NULL;
    char* szParameter = new char[64];
    double dParameter;
    unsigned int uiParameter;
    
    pFile = fopen("./defaults.dat", "r");
    
    // Error Check - No Parameters
    if (pFile == NULL)
    {
        fprintf(stderr, "No default parameters found!\n");
        return error;
    }
    
    //  Setting Variables to Defaults (Can be absent from input file)
    //  maxCoveringAreaFraction
   // maxCoveringAreaFraction = 0.0;
    
    
    
    // Error Check - Initial Conditions Found
    printf("File containing default parameters found.\n");
    
    // Open The File and Load Parameters
    while (fscanf(pFile, "%s", szParameter) > 0)
    {
        // Parameter:  Type of Boundary
        // ------------------------------------
        // ------------------------------------
        if (strcmp(szParameter, "box") == 0)
        {
            if (fscanf(pFile, " := %s", szParameter) > 0)
            {
                if (strcmp(szParameter, "hard") == 0)
                {
                    BoxConfinement = confinement_t::hard;
                }
                else if (strcmp(szParameter, "periodic") == 0)
                {
                    BoxConfinement = confinement_t::periodic;
                }
                else if (strcmp(szParameter, "soft") == 0)
                {
                    BoxConfinement = confinement_t::soft;
                }
                else if (strcmp(szParameter, "hexagonal") == 0)
                {
                    BoxConfinement = confinement_t::hexagonal;
                }
                else
                {
                    fprintf(stderr, "Parameter \"%s\" for option -box is not valid!\n",
                            szParameter);
                }
            }
        }
        //  Parameter:  Simulation Growth Rate
        // ------------------------------------
        // ------------------------------------
        else if (strcmp(szParameter, "gr") == 0)
        {
            if (fscanf(pFile, " := %lf", &dParameter) > 0)
            {
                GrowthRate = dParameter;
                printf("Set the growth rate to %le\n", GrowthRate);
            }
        }
        //  Parameter:  Acceptable Area of Overlap - Global
        // ------------------------------------
        // ------------------------------------
        else if (strcmp(szParameter, "mo") == 0)
        {
            if (fscanf(pFile, " := %lf", &dParameter) > 0)
            {
                MaxOverlap = dParameter;
                printf("Set the maximum overlap to %le\n", MaxOverlap);
            }
        }
        // ------------------------------------
        // ------------------------------------
        else if (strcmp(szParameter, "ns") == 0)
        {
            if (fscanf(pFile, " := %u", &uiParameter) > 0)
            {
                NumberOfShakes = uiParameter;
                printf("Set the number of shakes to %d\n", NumberOfShakes);
            }
        }
        //  Parameter:  Define Rotational Shaking - Amount of Particle Spin
        // ------------------------------------
        // ------------------------------------
        else if (strcmp(szParameter, "osa") == 0)
        {
            if (fscanf(pFile, " := %lf", &dParameter) > 0)
            {
                OrientationShakingAmplitude = dParameter;
                printf("Set the orientational shaking amplitude to %le\n",
                       OrientationShakingAmplitude);
            }
        }
        // ------------------------------------
        // ------------------------------------
        else if (strcmp(szParameter, "osc") == 0)
        {
            if (fscanf(pFile, " := %u", &uiParameter) > 0)
            {
                MaxOptimizeShakeCycles = uiParameter;
                printf("Set the maximum number of optimization and shake cycles to "
                       "%d\n", MaxOptimizeShakeCycles);
            }
        }
        // ------------------------------------
        // ------------------------------------
        else if (strcmp(szParameter, "pifo") == 0)
        {
            PrintInitialFinalOnly = true;
            printf("Only first and last configuration will be printed.\n");
        }
        // ------------------------------------
        // ------------------------------------
        else if (strcmp(szParameter, "sa") == 0)
        {
            if (fscanf(pFile, " := %lf", &dParameter) > 0)
            {
                ShakingAmplitude = dParameter;
                printf("Set the shaking amplitude to %le\n", ShakingAmplitude);
            }
        }
        //   Initial Conditions Passed by User
        // ------------------------------------
        // ------------------------------------
        else if (strcmp(szParameter, "boxFile") == 0)
        {
            if (fscanf(pFile, " := %s", BoundaryFile) > 0)
            {
                //BoundaryFile = szParameter;
                printf("Using Boundary Polygon from \"%s\"\n", BoundaryFile);
            }
        }
        
        
        // ------------------------------------
        // ------------------------------------
        else if (strcmp(szParameter, "polyLink") == 0)
        {
            if (fscanf(pFile, " := %s", ParticlesFile) > 0)
            {
                
                //ParticlesFile = szParameter;
                printf("Using Polygon Link Files from \"%s\"\n", ParticlesFile);
            }
        }
        // ------------------------------------
        // ------------------------------------
        else if (strcmp(szParameter, "maxPF") == 0)
        {
            if (fscanf(pFile, " := %lf", &maxCoveringAreaFraction) > 0)
            {
                 printf("Expanding until Covering Area Fraction of  %le \n", maxCoveringAreaFraction);
            }
        }
        // ------------------------------------
        // ------------------------------------
        else if (strcmp(szParameter, "outdir") == 0)
        {
            if (fscanf(pFile, " := %s", outputPath) > 0)
            {
               // outputPath = szParameter;
                printf("Setting Output Path as: \"%s\"\n", outputPath);
            }
        }
        // ------------------------------------
        else
        {
            printf("Unknown option %s\n", szParameter);
        }
        
    }
    
    fclose(pFile);
    delete[] szParameter;
    
    printf("Initialization Routine Complete.\n");
    return error;
};






//===========================================================
//===========================================================
int countlines(char *filename){

    unsigned int number_of_lines = 0;

    
    std::ifstream f(filename);
    std::string line;
    for (int i = 0; std::getline(f, line); ++i) {number_of_lines++ ;}
    
    return number_of_lines;
};


//===========================================================
//===========================================================
//   Load a single polygon file
//
std::vector <std::vector<double> > loadPolygonVertex(char *filename)
{
    // Temp Variable for {X,Y} positions
    double tempX, tempY;
    unsigned int i;
    unsigned int NumOfVert=0;
    unsigned int X=0;
    unsigned int Y=1;
    // Set Array Size
    std::vector<double> PolygonVertX;
    std::vector<double> PolygonVertY;
    
    std::vector <std::vector<double> > PolygonVertexList;

    
    PolygonVertX.clear();  PolygonVertY.clear(); PolygonVertexList.clear();

    FILE *pPolygonFile = fopen(filename, "r");
    
    
     std::cout << "Called loadPolygonVertex "<< filename <<" Successfully. ";
    
        //  Test if File Exists (==Null means No File Opened)
        if (pPolygonFile == NULL)
        {
             std::cout << "Failure to open:   \n";
        }
        else
        {

            while (fscanf(pPolygonFile, "%le %le\n", &tempX, &tempY) == 2)
            {
                PolygonVertX.push_back(tempX);
                PolygonVertY.push_back(tempY);
                NumOfVert++;
            };
                fclose(pPolygonFile);
            
            
            if (  fabs(PolygonVertX[0] - PolygonVertX[NumOfVert-1])>0.00001 ||
                fabs(PolygonVertY[0] - PolygonVertY[NumOfVert-1])>0.00001 )
            {
                PolygonVertX.push_back(PolygonVertX[0]);
                PolygonVertY.push_back(PolygonVertY[0]);
                NumOfVert++;
                std::cout << " Closing the polygon! \n";
            };
        };
            
         //   std::cout << "Polygon Loaded! \n";
            PolygonVertexList.resize( NumOfVert , std::vector<double>( 2 , 0.0 ) );
 
    
             for (int vertID =0; vertID <NumOfVert; vertID++ )
             {
             PolygonVertexList[vertID][X]=PolygonVertX[vertID];
             PolygonVertexList[vertID][Y]=PolygonVertY[vertID];

             };
    
     std::cout << "Returning to Main! \n";
     return PolygonVertexList;
};




//===========================================================
//===========================================================
input_polygons_struct loadPolygonVertexFromPolyLink(char *polyLink_Filename)
{
    int seperatePolygonFiles = countlines(polyLink_Filename);
    int nLoadedPolys = -1;  // Actual Nonzero Number of Particles found.
    
    // Define input beads
    std::vector <std::vector<double> > loadingPoly;
    loadingPoly.clear();
    
    // Define input beads
    std::vector<Particle> read_globalPolygonList;
    
    // Structure to Return Input Polygons, #of_Objects, and ParticleToPolygon List
    input_polygons_struct return_parameter;
    return_parameter.read_nPolySpecies.clear();
    return_parameter.read_inputPolygons.clear();
    return_parameter.read_polySpecifList.clear();
    return_parameter.ErrorLoadingQ = false;
    
    std::cout << "PolygonLinkFile: correct countlines?  "<< seperatePolygonFiles << std::endl;
    
    char polygonLinkName[512];
    int number_of_typePolygons[seperatePolygonFiles];

    
    
    FILE *PolygonListFile = fopen(polyLink_Filename, "r");
    for (int i = 0; i<seperatePolygonFiles; i++)
    {
        
        
        // Read each Line of the File : #ofPolys  \t  Path_to_Poly
        fscanf(PolygonListFile, "%d \t %s", &number_of_typePolygons[i], polygonLinkName) ;
        
        if (number_of_typePolygons[i] >=1){
            nLoadedPolys++;
        std::cout << "PolygonLinkFile: " <<  polygonLinkName << " with "
        << number_of_typePolygons[i]<<" nBeads "<<std::endl;
            
        loadingPoly =  loadPolygonVertex(polygonLinkName);
            if (loadingPoly.size() < 3)
            {
                return_parameter.ErrorLoadingQ = true;
            };
            

        return_parameter.read_inputPolygons.push_back(Particle());
        return_parameter.read_inputPolygons[i].Initialize( nLoadedPolys, loadingPoly);
        return_parameter.read_nPolySpecies.push_back(number_of_typePolygons[i]);
        
        
        
        for (int ii =0; ii<number_of_typePolygons[i]; ii++)
        {return_parameter.read_polySpecifList.push_back(nLoadedPolys);};

        };
        
    };
    
    
    fclose(PolygonListFile);
  
    
    std::cout << "Polydisperse Particles Input Correctly - Returning to Main!"<<std::endl;
    return return_parameter;

};
