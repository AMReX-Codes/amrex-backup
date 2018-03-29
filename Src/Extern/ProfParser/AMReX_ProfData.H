// ---------------------------------------------------------------
// ProfData.H
// ---------------------------------------------------------------
#ifndef _PROFDATA_H_
#define _PROFDATA_H_

//#include <AMReX_Array.H>
//#include <AMReX_Vector.H>

//#include <cstdarg>
//#include <string>

#ifdef BL_USE_PROFPARSER
//#include <map>
#include <AMReX_BLProfUtilities.H>
#include <AMReX_BLProfStats.H>
//#include <BLProfParser.tab.H>
#include <AMReX_CommProfStats.H>
#include <AMReX_RegionsProfStats.H>
#endif


using std::string;

namespace amrex {

class ProfData { 

  public:

    ProfData();
    ProfData(const string &dirname);
//    ProfData();
    ~ProfData();

    void Init();
//    void InitRegionTimeRanges();

//    static void SetBatchMode();
//    static void Dispatch(DSRequestType requestType, DataServices *ds, ...);

//    static void SetFabOutSize(int iSize);

  private:

    // Dispatchable functions

  public:

    // information functions
    const string &GetDirName() const { return dirName; }
    bool ProfDataOk() { return bProfDataOk; }  

    // profiling functions
//#ifdef BL_USE_PROFPARSER
    bool ProfDataAvailable()   const { return bProfDataAvailable;   }
    bool RegionDataAvailable() const { return bRegionDataAvailable; }
    bool TraceDataAvailable()  const { return bTraceDataAvailable;  }
    bool CommDataAvailable()   const { return bCommDataAvailable;   }

    BLProfStats &GetBLProfStats()               { return blProfStats_H; }
    RegionsProfStats &GetRegionsProfStats()     { return regOutputStats_H; }
    CommProfStats &GetCommOutputStats()         { return commOutputStats_H; }
    const amrex::BoxArray &ProcBoxArray() const { return procBoxArray; }
//#endif

// Store files & calculation procs that will be locally processed. 

  private:
    string dirName;
    bool bProfDataOk;
/*  int dsArrayIndex;
    int numberOfUsers;
    int iWriteToLevel;

    static Vector<DataServices *> dsArray;
    static int dsArrayIndexCounter;
    static int dsFabOutSize;
    static bool dsBatchMode;
*/
//#ifdef BL_USE_PROFPARSER
    bool bProfDataAvailable;
    bool bRegionDataAvailable;
    bool bTraceDataAvailable;
    bool bCommDataAvailable;

    BLProfStats blProfStats_H;
    RegionsProfStats regOutputStats_H;
    CommProfStats commOutputStats_H;

    amrex::BoxArray procBoxArray;
    amrex::Vector<int> fileNumbers;     // List of local file numbers.
                                        // Each filename is stored in the database object.
                                        //   since there are separate names for each type.
//#endif
};

}

#endif
