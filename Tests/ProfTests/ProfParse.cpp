// --------------------------------------------------------------------------
// ProfParse.cpp
// --------------------------------------------------------------------------
//  this file tests the flex/bison parsing of profiling data.
// --------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <AMReX_BoxArray.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_BLProfUtilities.H>
using std::cout;
using std::endl;

// --------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    amrex::Initialize(argc, argv);

    bool RPSuccess;
    BLProfStats blProfStats_H;
    CommProfStats commProfStats_H;
    RegionsProfStats regProfStats_H;

    std::string filename = "bl_prof_D";
    RPSuccess = amrex::BcastAndParseFile(filename, &blProfStats_H);
    std::cout << "RPSuccess = " << RPSuccess << std::endl;
    std::cout << "NProcs = " << blProfStats_H.GetNProcs() << std::endl;

    filename = "bl_prof_H";
    RPSuccess = amrex::BcastAndParseFile(filename, &blProfStats_H);
    std::cout << "RPSuccess = " << RPSuccess << std::endl;
    std::cout << "NProcs = " << blProfStats_H.GetNProcs() << std::endl;

    filename = "bl_prof_H";
    RPSuccess = amrex::ParseFile(filename, &blProfStats_H);
    std::cout << "RPSuccess = " << RPSuccess << std::endl;
    std::cout << "NProcs = " << blProfStats_H.GetNProcs() << std::endl;


    amrex::Finalize();
    return 0;
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
