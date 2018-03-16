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

    // Assumption:  ioRank is fixed at rank = 0.

    // Goal: evenly distribute files across ranks, but ensure the i/o proc
    //       has the minimal amount, if applicable.
    //       Also, isolate the I/O file (isolatedFile) on the i/o proc.

    // Methodology: If isolated file is at the end of the list, distribute all files backwards
    //              across the ranks and leave the last 'fair share' of files for the io proc.

    //              Otherwise, Distribute files backwards across all available ranks except the
    //              ioRank. Once at the isolated file, place the next 'fair share' of files on
    //              the ioRank to equalize the balance. Then, continue from where you left off,
    //              including all processors until the files are distributed.

    long ioRank = 0;           // Fixed for now. Add possible accomodation for adjustment later.
    long numFiles = 65000;
    long isolatedFile = 1357;
    long numProcs = 128;
    amrex::Vector<long> ranks_for_file(numFiles, -1);
    amrex::Vector<amrex::Vector<long>> files_on_rank;
    amrex::Vector<amrex::Vector<long>> files_on_rank_B;

    files_on_rank.resize(numProcs);
    files_on_rank_B.resize(numProcs);
    for (int i=0; i<files_on_rank.size(); ++i)
    {
       files_on_rank[i].resize(0);
       files_on_rank_B[i].resize(0);
    }
   
    // Method to calculate rank for all files.
    // ---------------------------------------------------
    long fpr = numFiles/euu(numProcs-1);
    if (numFiles < (isolatedFile + fpr))            // Isolated file is near the end of the list.
    {
       for (long i(0); i<(numFiles-fpr); ++i)       // Files before the isolated group.
       {
          ranks_for_file[i] = (numProcs-1) - (i%(numProcs-1));
       }
       for (long i=numFiles-fpr; i<numFiles; ++i)   // The isolated group.
       {
          ranks_for_file[i] = ioRank;
       }
    }
    else
    {
       for (long i=0; i< isolatedFile; ++i)         // All files before the isolated file.
       {
          ranks_for_file[i] = (numProcs-1) - (i%(numProcs-1));
       }
       for (long i=isolatedFile; i< isolatedFile + isolatedFile/(numProcs-1); ++i)   // Catch up the IO proc
       {
          ranks_for_file[i] = ioRank;
       }
       for (long i=(isolatedFile+(isolatedFile/numProcs)); i<numFiles; ++i)     // All files once all ranks are involved.
       {
          ranks_for_file[i] = (numProcs - 1) - (i%numProcs);
       }
    }

    // Output data & check results
    // ---------------------------------------------------
    for (long i=0; i<ranks_for_file.size(); ++i)
    {
       std::cout << "file # " << i << " to rank " << ranks_for_file[i] << std::endl;
       if (ranks_for_file[i] == -1) {
           amrex::Abort("Error: File not assigned a rank.");
       }
       if ((i == isolatedFile) && (ranks_for_file[i] != ioRank)){
           amrex::Abort("Error: Isolated file not placed on the ioRank.");
       }
       files_on_rank[ranks_for_file[i]].push_back(i);
    }

    std::cout << " =========================" << std::endl;

    for (long i=0; i<files_on_rank.size(); ++i)
    {
       std::cout << "files on rank " << i << " = " << files_on_rank[i].size() << std::endl;
    }

    // Method to calculate based file number based on local rank number.
    // ----------------------------------------------------------------
    long numLocalFiles = -1; 
    long numLocalFilesBeforeIso = -1;
    long endingRank = 0;

    if (rank >= ranks - (files%ranks))
    {
      numLocalFiles = 
    } 


    if (numFiles < (isolatedFile + fpr)     // Isolated file is at the end.
    {
      if (rank == ioRank)                   // IO proc gets the last bunch of files.
      {
         for (int i = numFiles-(numFiles/(numProcs-1)); i < files; ++i)
         {
           files_on_rank_B[i].push_back(i);
         }
      }
      else
      {
         for (int i=0; i<numFiles/
      }
    }
    else
    {
      if (rank == ioRank)
      {

      }
    }

    // Comparison of the two methods.
    // -------------------------------
    for (long i=0; i<files_on_rank.size(); ++i)
    {
       for (long j=0; j<files_on_rank[i].size(); ++j)
          if (files_on_rank[i][j] != files_on_rank_B[i][j]){
            amrex::Abort("Error: Methods don't agree");
          }
    }
    std::cout << "Methods match." << std::endl;

    amrex::Finalize();
    return 0;
}
// --------------------------------------------------------------------------
// --------------------------------------------------------------------------
