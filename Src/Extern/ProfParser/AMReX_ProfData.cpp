// ---------------------------------------------------------------
// ProfData.cpp
// ---------------------------------------------------------------
#ifndef _PROFDATA_CPP_
#define _PROFDATA_CPP_

//#include <AMReX_AmrvisConstants.H>
#include <AMReX_ProfData.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_BLProfUtilities.H>
//#include <AMReX_XYPlotDataList.H>

#ifdef BL_USE_PROFPARSER
//#include <AMReX_BLWritePlotFile.H>
#include <AMReX_BLProfStats.H>
//#include <AMReX_CommProfStats.H>
//#include <AMReX_RegionsProfStats.H>
//#include <iomanip>
//#include <cstdarg>
#endif

#include <iostream>
#include <fstream>
#include <cstdio>

using std::ios;
using std::ofstream;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;

namespace amrex {

#ifdef BL_USE_PROFPARSER
namespace {
  const int XDIR(0);
  const int YDIR(1);
  const int ZDIR(2);
  const int NTIMESLOTS(25600);
}


#define SHOWVAL(val) { cout << #val << " = " << val << endl; }
#endif

// ---------------------------------------------------------------
ProfData::ProfData()
{

  dirName = "";
  bProfDataOk =  false;
  bProfDataAvailable = false;
  bRegionDataAvailable = false;
  bTraceDataAvailable = false;
  bCommDataAvailable = false;

}
// ---------------------------------------------------------------
ProfData::ProfData(const string &dirname)
             : dirName(dirname), bProfDataOk(false)
{

  bProfDataAvailable = false;
  bRegionDataAvailable = false;
  bTraceDataAvailable = false;
  bCommDataAvailable = false;

  Init();

}


// ---------------------------------------------------------------
~ProfData::ProfData() {
}

// ---------------------------------------------------------------
void ProfData::Init() {

  BL_PROFILE("ProfData::Init()");

  bProfDataOk = false;

#if (BL_SPACEDIM == 2)
    int  displayProc(0);  // TEMPORARY. Needs to be passed into init in the future. Data rank displayed to screen.

 
    bool bIOP(ParallelDescriptor::IOProcessor());
//    int  myProc(ParallelDescriptor::MyProc());
//    int  nProcs(ParallelDescriptor::NProcs());

    BLProfStats::SetDirName(dirName);

    // -------- parse the main blprof header file.  I/O reads and distributes for local parsing. 
    if(bIOP) { cout << "Parsing main blprof header file." << endl; }
    string blpFileName_H(BLProfStats::GetProfFilePrefix()+"_H");   // bl_prof_H
    string blpFullFileName_H(dirName + '/' + blpFileName_H);
    if( amrex::BcastAndParseFile(blpFullFileName_H, &blProfStats_H))
    {
      bProfDataAvailable = true;
      fileNumbers = CalcFileNumbers(displayProc, blProfStats_H.GetNOutFiles(), BLProfStats::GetProfFilePrefix()+"_D_");
      blProfStats_H.SetLocalData();
    }
    else
    {
      cerr << "ProfData::Init: 0: Cannot open file:    " << blpFullFileName_H << endl;
    }

    // -------- 




    DistributeWork

    fileNumbers = CalcFileNumbers(displayProc, blProfStats_H.GetNOutFiles(), BLProfStats::GetProfFilePrefix()+"_D_");
    


    // REDO STREAM INDEXES!!!
    // Create and store file names of all databases for BlockInit calls.
    // Store proc numbers in databases.

    // Build work box (cut region into "proc" boxes and combine as needed).
    // Read correct remaining headers.


    // -------- parse the main call stats header file.  everyone does this for now
/*
    if(bIOP) { cout << "Parsing main call stats header file." << endl; }
    string regPrefix_H(RegionsProfStats::GetRegFilePrefix()+"_H");
    std::string regFileName_H(fileName + '/' + regPrefix_H);
    if( ! (yyin = fopen(regFileName_H.c_str(), "r"))) {
      if(bIOP) {
        cerr << "DataServices::Init:  1:  Cannot open file:  " << regFileName_H << endl;
      }
      bRegionDataAvailable = false;
      bTraceDataAvailable  = false;
    } else {
      bRegionDataAvailable = true;
      bTraceDataAvailable  = true;
    }
    if(bTraceDataAvailable) {
      yyparse(&regOutputStats_H);
      fclose(yyin);
    } else {
    }

    // ---- make a box for distributing work
    int dataNProcs(BLProfStats::GetNProcs());
    Box procBox(IntVect(0, 0), IntVect(0, dataNProcs - 1));
    IntVect procMaxGrid(1, (dataNProcs / nProcs) + ((dataNProcs % nProcs) > 0 ? 1 : 0));
    BoxArray procBoxArrayTemp(procBox);
    procBoxArrayTemp.maxSize(procMaxGrid);
    // ---- now ensure the boxarray is nprocs long
    Vector<Box> procBoxes;
    int needMoreBoxes(nProcs - procBoxArrayTemp.size());
    for(int ipb(0); ipb < procBoxArrayTemp.size(); ++ipb) {
      Box b(procBoxArrayTemp[ipb]);
      if(needMoreBoxes) {
        Box chopBox(b.chop(YDIR, (b.smallEnd(YDIR) + b.bigEnd(YDIR)) / 2));
        procBoxes.push_back(chopBox);
        --needMoreBoxes;
      }
      procBoxes.push_back(b);
    }
    procBoxArray.resize(procBoxes.size());
    for(int i(0); i < procBoxes.size(); ++i) {
      procBoxArray.set(i, procBoxes[i]);
    }

    if(procBoxArray.size() != nProcs) {
      SHOWVAL(nProcs);
      SHOWVAL(dataNProcs);
      SHOWVAL(procBoxArray.size());
      if(bIOP) cout << "---- procBoxArray = " << procBoxArray << endl;
      amrex::Abort("procBoxArray::Error 0");
    }

    if(bTraceDataAvailable) {
      // -------- parse the data headers.  everyone does this for now
      if(bIOP) { cout << "Parsing data headers." << endl; }
      const Vector<string> &regHeaderFileNames = regOutputStats_H.GetHeaderFileNames();
      for(int i(0); i < regHeaderFileNames.size(); ++i) {
        std::string regFileName_H_nnnn(fileName + '/' + regHeaderFileNames[i]);
        if( ! (yyin = fopen(regFileName_H_nnnn.c_str(), "r"))) {
          if(bIOP) {
            cerr << "DataServices::Init:  2:  Cannot open file:  " << regFileName_H_nnnn
                 << " ... continuing." << endl;
          }
          continue;
        }
        BL_PROFILE_VAR("DataServices::Init(), parsing data headers.", yydheaders);
        yyparse(&regOutputStats_H);
        BL_PROFILE_VAR_STOP(yydheaders);

        fclose(yyin);
      }

      if(regOutputStats_H.TraceDataValid()) {
//        if(bIOP) {
//	  cout << "Calling InitRegionTimeRanges." << endl;
//	}
        RegionsProfStats::OpenAllStreams(fileName);
        Box myBox(procBoxArray[myProc]);
//        bRegionDataAvailable = regOutputStats_H.InitRegionTimeRanges(myBox);
        regOutputStats_H.SyncFNamesAndNumbers();
        RegionsProfStats::CloseAllStreams();
        regOutputStats_H.SetFNames(blProfStats_H.BLPFNames());
//        if(bIOP) {
//	  cout << "Finished InitRegionTimeRanges." << endl;
//	}
      } else {
        bTraceDataAvailable = false;
      }

    }
*/
/*
    bool bParseFilterFile(true);
    if(bParseFilterFile) {
      if(bIOP) cout << "Parsing filter file." << endl;
      ParseFilterFile();
    } 

    if( ! regOutputStats_H.TimeRangeInitialized()) {
      regOutputStats_H.InitFilterTimeRanges();
      if(bIOP) {
        cout << ">>>> timerangelist =" ;
        PrintTimeRangeList(regOutputStats_H.GetFilterTimeRanges()[0]);
      }

    }
*/
    // ----------------------------------------------- comm headers

    // -------- parse the main header file.  everyone does this for now
/*
    if(bIOP) { cout << "Parsing main comm header file." << endl; }
    std::string commPrefix_H(CommProfStats::GetCommFilePrefix()+"_H");
    std::string commFileName_H(fileName + '/' + commPrefix_H);
    if( ! (yyin = fopen(commFileName_H.c_str(), "r"))) {
      if(bIOP) {
        cerr << "DataServices::Init:  3:  Cannot open file:  " << commFileName_H << endl;
      }
      bCommDataAvailable = false;
    } else {
      bCommDataAvailable = true;
    }

    if(bCommDataAvailable) {
      yyparse(&commOutputStats_H);
      fclose(yyin);
    } else {
    }
*/
/*
    if(bRegionDataAvailable) {
      commOutputStats_H.SetRegionTimeRanges(regOutputStats_H.GetRegionTimeRanges());
      commOutputStats_H.SetFilterTimeRanges(regOutputStats_H.GetFilterTimeRanges());
    }
*/
    bProfDataOk = true;
#endif
}

// ---------------------------------------------------------------
ProfData::~ProfData() {
}

// ---------------------------------------------------------------
amrex::Vector<int> ProfData::CalcFileNumbers(int displayProc, int numFiles, std::string fileNamePrefix)
{
    // Goal: Evenly distribute files across ranks, but ensure the i/o proc
    //           has the minimal amount, if applicable.
    //       Also, isolate the I/O file (isolatedFile) on the i/o proc.

    // Methodology: Distribute the files forward on ranks starting with I/O proc + 1. 
    //              Swap the display file with the first file to be placed on the I/O proc.
    //              Unless display file would be on the I/O proc, then do nothing special. 
 
    amrex::Vector<int> local_fileNumbers(0);  // Explicitly define with size = 0 
    int ioRank(ParallelDescriptor::IOProcessorNumber());
    int nProcs(ParallelDescriptor::NProcs());
    int myRank(ParallelDescriptor::MyProc());

    // Everything is local. Don't bother with the calc.
    if ((nProcs == 1) || (numFiles == 1))
    {
       if (myRank == ioRank)
       {
          for (int i=0; i<numFiles; ++i)
          {
             local_fileNumbers.push_back(i);
          }
       }
       return local_fileNumbers;
    }

    // Parse the displayFile for the number.
    // Needed in the distribution calculation.
    std::string displayFileString( blProfStats_H.GetFileNameFromProc(displayProc) );
    std::stringstream displayFileNumStr ( displayFileString.substr(fileNamePrefix.size(), displayFileString.size()) );
    int displayFile;
    displayFileNumStr >> displayFile;

    // replacementFile = first file on ioRank.
    // file = first file on this rank. Will be incremented.
    int replacementFile = ( (nProcs-1) % nProcs ); 
    int file = ( (myRank+(nProcs-(ioRank+1))) % nProcs ); 
    int displayRank = ( (displayFile + ioRank + 1) % nProcs );

    if (myRank == ioRank)
    {
       if (displayRank != ioRank)
       {
          local_fileNumbers.push_back(displayFile);
          file += nProcs;
       }
       else
       {
          replacementFile = displayFile;
       }
    }

    while(file < numFiles)
    {
       local_fileNumbers.push_back(file);
       file += nProcs;
       if (file == displayFile)
       {
          local_fileNumbers.push_back(replacementFile);
          file += nProcs;
       }
    }

    return local_fileNumbers;
}
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
}  // namespace amrex 
#endif        // _PROFDATA_CPP_
