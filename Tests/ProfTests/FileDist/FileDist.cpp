// --------------------------------------------------------------------------
// FileDist.cpp
// --------------------------------------------------------------------------
//  this file tests the file distribution method for ProfVis.
// --------------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <AMReX_Vector.H>
//#include <AMReX_BoxArray.H>
//#include <AMReX_MultiFab.H>
//#include <AMReX_ParallelDescriptor.H>
//#include <AMReX_Utility.H>
//#include <AMReX_BLProfUtilities.H>
using std::cout;
using std::endl;

// --------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    amrex::Initialize(argc, argv);

    // Goal: evenly distribute files across ranks, but ensure the i/o proc
    //       has the minimal amount, if applicable.
    //       Also, isolate the I/O file (isolatedFile) on the i/o proc.

    // Methodology: Distribute the files forward on ranks starting with I/O proc + 1. 
    //              Swap the isolated file with the first file to be placed on the I/O proc.
    //              Unless isolated file would be on the I/O proc, then do nothing special. 

    int ioRank = 3;
    int numFiles = 34;
    int isolatedFile = 21;
    int numProcs = 11;
    amrex::Vector<int> ranks_for_file(numFiles, -1);
    amrex::Vector<amrex::Vector<int>> files_on_rank;

    files_on_rank.resize(numProcs);
    for (int i=0; i<files_on_rank.size(); ++i)
    {
       files_on_rank[i].resize(0);
    }

    // Method to calculate based file number based on local rank number.
    //   (Distributed, direct method).
    // ----------------------------------------------------------------
    if ((numProcs == 1) || (numFiles == 1))
    {
       for (int j=0; j<numProcs; ++j)     // Mock each rank's work.
       {
         if (j == ioRank)
         {
           // Serial version. No calc, just designate them all to IOrank.
           // Should be optimal unless there is I/O work to do besides the database lookup.
           // Then, possible to put "work" on a different rank than "I/O".
           for (int i=0; i<numFiles; ++i)
           {
              files_on_rank[ioRank].push_back(i);
           }
         }
       }
    }

    for (int j=0; j<numProcs; ++j)    // Mock each rank's work. 
    {
       int replacementFile = ( (numProcs-1) % numProcs );   // First file on ioRank.
 
       int file = ( (j+(numProcs-(ioRank+1))) % numProcs );    // First file on this rank. (maybe the same).
       int isolatedRank = ( (isolatedFile + ioRank + 1) % numProcs );

       if (j == ioRank)  
       {
         if (isolatedRank != ioRank)
         {
           files_on_rank[j].push_back(isolatedFile);
           file += numProcs;
         }
         else
         {
           replacementFile = isolatedFile;
         }
       }

       while(file < numFiles)
       {
         files_on_rank[j].push_back(file);
         file += numProcs;
         if (file == isolatedFile)
         {
            files_on_rank[j].push_back(replacementFile);
            file += numProcs;
         }
       } 
    }

    // Output data & check results
    // ---------------------------------------------------

    std::cout << " =========================" << std::endl;

    std::cout << "ioRank = " << ioRank << std::endl;
    std::cout << "numFiles = " << numFiles << std::endl;
    std::cout << "isolatedFile = " << isolatedFile << std::endl;
    std::cout << "numProcs = " << numProcs << std::endl;

    std::cout << " =========================" << std::endl;

    for (int i=0; i<files_on_rank.size(); ++i)
    {
      for (int j=0; j<files_on_rank[i].size(); ++j)
       {
          std::cout << "[" << i << " ][ " << j << " ] = " << files_on_rank[i][j] << std::endl;
          if (ranks_for_file[files_on_rank[i][j]] != -1)
          {
             amrex::Abort("***** Error: File assigned to multiple ranks.");
          }
          ranks_for_file[files_on_rank[i][j]] = i;
       }
    }

    std::cout << " =========================" << std::endl;


    for (int i=0; i<files_on_rank.size(); ++i)
    {
       std::cout << "files on rank " << i << " = " << files_on_rank[i].size() << std::endl;
       if (files_on_rank[i].size() < files_on_rank[ioRank].size())
       {
          amrex::Abort("***** Error: rank has less files than ioRank.");
       }
    }

    std::cout << " =========================" << std::endl;

    for (int i=0; i<ranks_for_file.size(); ++i)
    {
       std::cout << "file # " << i << " to rank " << ranks_for_file[i] << std::endl;
       if (ranks_for_file[i] == -1) {
           amrex::Abort("****** Error: File not assigned a rank.");
       }
       if ((i == isolatedFile) && (ranks_for_file[i] != ioRank)){
           amrex::Abort("****** Error: Isolated file not placed on the ioRank.");
       }
    }

    amrex::Finalize();
    return 0;
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
