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
    std::string filename;
    BLProfStats blProfStats_H;
    CommProfStats commProfStats_H;
    RegionsProfStats regProfStats_H;

    filename = "bl_prof_H";
    std::cout << "RPSuccess = " << amrex::BcastAndParseFile(filename, &blProfStats_H) << std::endl;
    std::cout << " ------------ " << std::endl;

    filename = "bl_call_stats_H";
    std::cout << "RPSuccess = " << amrex::BcastAndParseFile(filename, &blProfStats_H) << std::endl;
    std::cout << " ------------ " << std::endl;

    filename = "bl_comm_prof_H";
    std::cout << "RPSuccess = " << amrex::ParseFile(filename, &blProfStats_H) << std::endl;
    std::cout << " ------------ " << std::endl;

    filename = "bl_prof_H";
    std::cout << "RPSuccess = " << amrex::ParseFile(filename, &blProfStats_H) << std::endl;
    std::cout << " ------------ " << std::endl;

/*
    filename = "bl_prof_D";
    RPSuccess = amrex::BcastAndParseFile(filename, &blProfStats_H);
    std::cout << "RPSuccess = " << RPSuccess << std::endl;
    std::cout << "NProcs = " << blProfStats_H.GetNProcs() << std::endl;
    std::cout << " ------------ " << std::endl;
*/
    filename = "bl_prof_H";
    RPSuccess = amrex::BcastAndParseFile(filename, &blProfStats_H);
    std::cout << "RPSuccess = " << RPSuccess << std::endl;
//    std::cout << "NProcs = " << blProfStats_H.GetNProcs() << std::endl;
    std::cout << " ------------ " << std::endl;
/*
    filename = "bl_prof_H";
    RPSuccess = amrex::BcastAndParseFile(filename, &blProfStats_H);
    std::cout << "RPSuccess = " << RPSuccess << std::endl;
    std::cout << "NProcs = " << blProfStats_H.GetNProcs() << std::endl;
    std::cout << " ------------ " << std::endl;

    filename = "bl_prof_D";
    RPSuccess = amrex::ParseFile(filename, &blProfStats_H);
    std::cout << "RPSuccess = " << RPSuccess << std::endl;
    std::cout << "NProcs = " << blProfStats_H.GetNProcs() << std::endl;
    std::cout << " ------------ " << std::endl;

    filename = "bl_prof_H";
    RPSuccess = amrex::ParseFile(filename, &blProfStats_H);
    std::cout << "RPSuccess = " << RPSuccess << std::endl;
    std::cout << "NProcs = " << blProfStats_H.GetNProcs() << std::endl;
    std::cout << " ------------ " << std::endl;
*/
    filename = "bl_prof_H";
    RPSuccess = amrex::ParseFile(filename, &blProfStats_H);
    std::cout << "RPSuccess = " << RPSuccess << std::endl;
//    std::cout << "NProcs = " << blProfStats_H.GetNProcs() << std::endl;
    std::cout << " ------------ " << std::endl;
    amrex::Finalize();
    return 0;
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
