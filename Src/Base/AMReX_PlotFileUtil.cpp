
#include <fstream>
#include <iomanip>

#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>

#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#endif

#ifdef AMREX_USE_HDF5
#include "AMReX_HDF5.H"
#include "AMReX_MultiFabUtil.H"
#include "hdf5.h"
#endif

namespace amrex {

std::string LevelPath (int level, const std::string &levelPrefix)
{
    return Concatenate(levelPrefix, level, 1);  // e.g., Level_5
}

std::string MultiFabHeaderPath (int level,
                                const std::string &levelPrefix,
                                const std::string &mfPrefix)
{
    return LevelPath(level, levelPrefix) + '/' + mfPrefix;  // e.g., Level_4/Cell
}

std::string LevelFullPath (int level,
                           const std::string &plotfilename,
                           const std::string &levelPrefix)
{
    std::string r(plotfilename);
    if ( ! r.empty() && r.back() != '/') {
	r += '/';
    }
    r += LevelPath(level, levelPrefix);  // e.g., plt00005/Level_5
    return r;
}

std::string MultiFabFileFullPrefix (int level,
                                    const std::string& plotfilename,
                                    const std::string &levelPrefix,
                                    const std::string &mfPrefix)
{
    std::string r(plotfilename);
    if ( ! r.empty() && r.back() != '/') {
	r += '/';
    }
    r += MultiFabHeaderPath(level, levelPrefix, mfPrefix);
    return r;
}


void
PreBuildDirectorHierarchy (const std::string &dirName,
                           const std::string &subDirPrefix,
                           int nSubDirs, bool callBarrier)
{
  UtilCreateCleanDirectory(dirName, false);  // ---- dont call barrier
  for(int i(0); i < nSubDirs; ++i) {
    const std::string &fullpath = LevelFullPath(i, dirName);
    UtilCreateCleanDirectory(fullpath, false);  // ---- dont call barrier
  }

  if(callBarrier) {
    ParallelDescriptor::Barrier();
  }
}


void
WriteGenericPlotfileHeader (std::ostream &HeaderFile,
                            int nlevels,
                            const Vector<BoxArray> &bArray,
                            const Vector<std::string> &varnames,
                            const Vector<Geometry> &geom,
                            Real time,
                            const Vector<int> &level_steps,
                            const Vector<IntVect> &ref_ratio,
                            const std::string &versionName,
                            const std::string &levelPrefix,
                            const std::string &mfPrefix)
{
        BL_PROFILE("WriteGenericPlotfileHeader()");

        BL_ASSERT(nlevels <= bArray.size());
        BL_ASSERT(nlevels <= geom.size());
        BL_ASSERT(nlevels <= ref_ratio.size()+1);
        BL_ASSERT(nlevels <= level_steps.size());

        int finest_level(nlevels - 1);

	HeaderFile.precision(17);

	// ---- this is the generic plot file type name
        HeaderFile << versionName << '\n';

        HeaderFile << varnames.size() << '\n';

        for (int ivar = 0; ivar < varnames.size(); ++ivar) {
	    HeaderFile << varnames[ivar] << "\n";
        }
        HeaderFile << AMREX_SPACEDIM << '\n';
        HeaderFile << time << '\n';
        HeaderFile << finest_level << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            HeaderFile << geom[0].ProbLo(i) << ' ';
	}
        HeaderFile << '\n';
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            HeaderFile << geom[0].ProbHi(i) << ' ';
	}
        HeaderFile << '\n';
        for (int i = 0; i < finest_level; ++i) {
            HeaderFile << ref_ratio[i][0] << ' ';
	}
        HeaderFile << '\n';
	for (int i = 0; i <= finest_level; ++i) {
	    HeaderFile << geom[i].Domain() << ' ';
	}
        HeaderFile << '\n';
        for (int i = 0; i <= finest_level; ++i) {
            HeaderFile << level_steps[i] << ' ';
	}
        HeaderFile << '\n';
        for (int i = 0; i <= finest_level; ++i) {
            for (int k = 0; k < AMREX_SPACEDIM; ++k) {
                HeaderFile << geom[i].CellSize()[k] << ' ';
	    }
            HeaderFile << '\n';
        }
        HeaderFile << (int) geom[0].Coord() << '\n';
        HeaderFile << "0\n";

	for (int level = 0; level <= finest_level; ++level) {
	    HeaderFile << level << ' ' << bArray[level].size() << ' ' << time << '\n';
	    HeaderFile << level_steps[level] << '\n';

	    for (int i = 0; i < bArray[level].size(); ++i)
	    {
		const Box &b(bArray[level][i]);
		RealBox loc = RealBox(b, geom[level].CellSize(), geom[level].ProbLo());
		for (int n = 0; n < AMREX_SPACEDIM; ++n) {
		    HeaderFile << loc.lo(n) << ' ' << loc.hi(n) << '\n';
		}
	    }

	    HeaderFile << MultiFabHeaderPath(level, levelPrefix, mfPrefix) << '\n';
	}
}


void
WriteMultiLevelPlotfile (const std::string& plotfilename, int nlevels,
                         const Vector<const MultiFab*>& mf,
                         const Vector<std::string>& varnames,
                         const Vector<Geometry>& geom, Real time, const Vector<int>& level_steps,
                         const Vector<IntVect>& ref_ratio,
                         const std::string &versionName,
                         const std::string &levelPrefix,
                         const std::string &mfPrefix,
                         const Vector<std::string>& extra_dirs)
{
    BL_PROFILE("WriteMultiLevelPlotfile()");

    BL_ASSERT(nlevels <= mf.size());
    BL_ASSERT(nlevels <= geom.size());
    BL_ASSERT(nlevels <= ref_ratio.size()+1);
    BL_ASSERT(nlevels <= level_steps.size());
    BL_ASSERT(mf[0]->nComp() == varnames.size());

    int finest_level = nlevels-1;

//    int saveNFiles(VisMF::GetNOutFiles());
//    VisMF::SetNOutFiles(std::max(1024,saveNFiles));

    bool callBarrier(false);
    PreBuildDirectorHierarchy(plotfilename, levelPrefix, nlevels, callBarrier);
    if (!extra_dirs.empty()) {
        for (const auto& d : extra_dirs) {
            const std::string ed = plotfilename+"/"+d;
            amrex::PreBuildDirectorHierarchy(ed, levelPrefix, nlevels, callBarrier);
        }
    }
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor()) {
      VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
      std::string HeaderFileName(plotfilename + "/Header");
      std::ofstream HeaderFile;
      HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
      HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
	                                      std::ofstream::trunc |
                                              std::ofstream::binary);
      if( ! HeaderFile.good()) {
        FileOpenFailed(HeaderFileName);
      }

      Vector<BoxArray> boxArrays(nlevels);
      for(int level(0); level < boxArrays.size(); ++level) {
	boxArrays[level] = mf[level]->boxArray();
      }

      WriteGenericPlotfileHeader(HeaderFile, nlevels, boxArrays, varnames,
                                 geom, time, level_steps, ref_ratio, versionName, levelPrefix, mfPrefix);
    }


    for (int level = 0; level <= finest_level; ++level)
    {
        const MultiFab* data;
        std::unique_ptr<MultiFab> mf_tmp;
        if (mf[level]->nGrow() > 0) {
            mf_tmp.reset(new MultiFab(mf[level]->boxArray(),
                                      mf[level]->DistributionMap(),
                                      mf[level]->nComp(), 0, MFInfo().SetArena(The_Pinned_Arena()),
                                      mf[level]->Factory()));
            MultiFab::Copy(*mf_tmp, *mf[level], 0, 0, mf[level]->nComp(), 0);
            data = mf_tmp.get();
        } else {
            data = mf[level];
        }
	VisMF::Write(*data, MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix));
    }

//    VisMF::SetNOutFiles(saveNFiles);
}

// write a plotfile to disk given:
// -plotfile name
// -vector of MultiFabs
// -vector of Geometrys
// variable names are written as "Var0", "Var1", etc.    
// refinement ratio is computed from the Geometry vector
// "time" and "level_steps" are set to zero
void WriteMLMF (const std::string &plotfilename,
                const Vector<const MultiFab*>& mf,
                const Vector<Geometry> &geom)
{
    int nlevs = mf.size();
    int ncomp = mf[0]->nComp();

    // variables names are "Var0", "Var1", etc.
    Vector<std::string> varnames(ncomp);
    for (int i=0; i<ncomp; ++i) {
        varnames[i] = "Var" + std::to_string(i);
    }

    // compute refinement ratio by looking at hi coordinates of domain at each level from
    // the geometry object
    Vector<IntVect> ref_ratio(nlevs-1);
    for (int i = 0; i < nlevs-1; ++i) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            int rr = (geom[i+1].Domain()).bigEnd(d)/(geom[i].Domain()).bigEnd(d);
            ref_ratio[i][d] = rr;
        }
    }

    // set step_array to zero
    Vector<int> step_array(nlevs,0);

    // set time to zero
    Real time = 0.;
    
    WriteMultiLevelPlotfile(plotfilename, nlevs, mf, varnames,
                            geom, time, step_array, ref_ratio);   
    
}    


void
WriteMultiLevelPlotfileHeaders (const std::string & plotfilename, int nlevels,
                                const Vector<const MultiFab*> & mf,
                                const Vector<std::string>     & varnames,
                                const Vector<Geometry>        & geom,
                                Real time, const Vector<int>  & level_steps,
                                const Vector<IntVect>     & ref_ratio,
                                const std::string         & versionName,
                                const std::string         & levelPrefix,
                                const std::string         & mfPrefix,
                                const Vector<std::string> & extra_dirs)
{
    BL_PROFILE("WriteMultiLevelPlotfile()");

    int finest_level = nlevels-1;

    bool callBarrier(false);
    PreBuildDirectorHierarchy(plotfilename, levelPrefix, nlevels, callBarrier);
    if (!extra_dirs.empty()) {
        for (const auto& d : extra_dirs) {
            const std::string ed = plotfilename+"/"+d;
            amrex::PreBuildDirectorHierarchy(ed, levelPrefix, nlevels, callBarrier);
        }
    }
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor()) {
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::string HeaderFileName(plotfilename + "/Header");
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                        std::ofstream::trunc |
                        std::ofstream::binary);
        if( ! HeaderFile.good()) {
            FileOpenFailed(HeaderFileName);
        }

        Vector<BoxArray> boxArrays(nlevels);
        for(int level(0); level < boxArrays.size(); ++level) {
            boxArrays[level] = mf[level]->boxArray();
        }

        WriteGenericPlotfileHeader(HeaderFile, nlevels, boxArrays, varnames,
                                   geom, time, level_steps, ref_ratio, versionName, levelPrefix, mfPrefix);
    }

    for (int level = 0; level <= finest_level; ++level) {
        const MultiFab * data;
        data = mf[level];
        VisMF::WriteOnlyHeader(*data, MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix));
    }

}



void
WriteSingleLevelPlotfile (const std::string& plotfilename,
                          const MultiFab& mf, const Vector<std::string>& varnames,
                          const Geometry& geom, Real time, int level_step,
                          const std::string &versionName,
                          const std::string &levelPrefix,
                          const std::string &mfPrefix,
                          const Vector<std::string>& extra_dirs)
{
    Vector<const MultiFab*> mfarr(1,&mf);
    Vector<Geometry> geomarr(1,geom);
    Vector<int> level_steps(1,level_step);
    Vector<IntVect> ref_ratio;

    WriteMultiLevelPlotfile(plotfilename, 1, mfarr, varnames, geomarr, time,
                            level_steps, ref_ratio, versionName, levelPrefix, mfPrefix, extra_dirs);
}

#ifdef AMREX_USE_EB
void
EB_WriteSingleLevelPlotfile (const std::string& plotfilename,
                             const MultiFab& mf, const Vector<std::string>& varnames,
                             const Geometry& geom, Real time, int level_step,
                             const std::string &versionName,
                             const std::string &levelPrefix,
                             const std::string &mfPrefix,
                             const Vector<std::string>& extra_dirs)
{
    Vector<const MultiFab*> mfarr(1,&mf);
    Vector<Geometry> geomarr(1,geom);
    Vector<int> level_steps(1,level_step);
    Vector<IntVect> ref_ratio;

    EB_WriteMultiLevelPlotfile(plotfilename, 1, mfarr, varnames, geomarr, time,
                               level_steps, ref_ratio, versionName, levelPrefix, mfPrefix, extra_dirs);
}

void
EB_WriteMultiLevelPlotfile (const std::string& plotfilename, int nlevels,
                            const Vector<const MultiFab*>& mf,
                            const Vector<std::string>& varnames,
                            const Vector<Geometry>& geom, Real time, const Vector<int>& level_steps,
                            const Vector<IntVect>& ref_ratio,
                            const std::string &versionName,
                            const std::string &levelPrefix,
                            const std::string &mfPrefix,
                            const Vector<std::string>& extra_dirs)
{
    BL_PROFILE("WriteMultiLevelPlotfile()");

    BL_ASSERT(nlevels <= mf.size());
    BL_ASSERT(nlevels <= geom.size());
    BL_ASSERT(nlevels <= ref_ratio.size()+1);
    BL_ASSERT(nlevels <= level_steps.size());
    BL_ASSERT(mf[0]->nComp() == varnames.size());
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(mf[0]->hasEBFabFactory(),
                                     "EB_WriteMultiLevelPlotfile: does not have EB Factory");

    int finest_level = nlevels-1;

//    int saveNFiles(VisMF::GetNOutFiles());
//    VisMF::SetNOutFiles(std::max(1024,saveNFiles));

    bool callBarrier(false);
    PreBuildDirectorHierarchy(plotfilename, levelPrefix, nlevels, callBarrier);
    if (!extra_dirs.empty()) {
        for (const auto& d : extra_dirs) {
            const std::string ed = plotfilename+"/"+d;
            amrex::PreBuildDirectorHierarchy(ed, levelPrefix, nlevels, callBarrier);
        }
    }
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor()) {
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::string HeaderFileName(plotfilename + "/Header");
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
	                                        std::ofstream::trunc |
                                                std::ofstream::binary);
        if( ! HeaderFile.good()) {
            FileOpenFailed(HeaderFileName);
        }

        Vector<BoxArray> boxArrays(nlevels);
        for(int level(0); level < boxArrays.size(); ++level) {
            boxArrays[level] = mf[level]->boxArray();
        }

        Vector<std::string> vn = varnames;
        vn.push_back("vfrac");
        WriteGenericPlotfileHeader(HeaderFile, nlevels, boxArrays, vn,
                                   geom, time, level_steps, ref_ratio, versionName,
                                   levelPrefix, mfPrefix);

        for (int lev = 0; lev < nlevels; ++lev) {
            HeaderFile << "1.0e-6\n";
        }
    }


    for (int level = 0; level <= finest_level; ++level)
    {
        const int nc = mf[level]->nComp();
        MultiFab mf_tmp(mf[level]->boxArray(),
                        mf[level]->DistributionMap(),
                        nc+1, 0);
        MultiFab::Copy(mf_tmp, *mf[level], 0, 0, nc, 0);
        auto const& factory = dynamic_cast<EBFArrayBoxFactory const&>(mf[level]->Factory());
        MultiFab::Copy(mf_tmp, factory.getVolFrac(), 0, nc, 1, 0);
	VisMF::Write(mf_tmp, MultiFabFileFullPrefix(level, plotfilename, levelPrefix, mfPrefix));
    }

//    VisMF::SetNOutFiles(saveNFiles);
}

#endif


#ifdef AMREX_USE_HDF5

void writeLevelAttrHDF5 (H5& h5, int level, const Geometry& geom, const IntVect& ref_ratio,
                         Real level_time, Real level_dt, Real level_min_dt, int n_cycle,
                         int level_steps, int level_count,
                         const BoxArray& grids, const DistributionMapping& dmap)
{
    // level
    h5.writeAttribute("level", level, H5T_NATIVE_INT);
    
    // cell spacing
    const Real* dx = geom.CellSize();
    hid_t realvect_id = makeH5RealVec();
    real_h5_t vec_dx = writeH5RealVec(dx);
    h5.writeAttribute("vec_dx", vec_dx, realvect_id);
    h5.writeAttribute("dx", vec_dx.x, H5T_NATIVE_DOUBLE);
    H5Tclose(realvect_id);
    
    // refinement ratio
    hid_t intvect_id = makeH5IntVec();
    int rr[AMREX_SPACEDIM];
    const int* vrr = ref_ratio.getVect();
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        rr[i] = vrr[i];
    }
    h5.writeAttribute("vec_ref_ratio", rr, intvect_id);
    h5.writeAttribute("ref_ratio", rr[0], H5T_NATIVE_INT);

    // data time stamp
    double time = level_time;
    h5.writeAttribute("time", time, H5T_NATIVE_DOUBLE);

    // time step size
    double dt = level_dt;
    h5.writeAttribute("dt", dt, H5T_NATIVE_DOUBLE);

    double min_dt = level_min_dt;
    h5.writeAttribute("min_dt", min_dt, H5T_NATIVE_DOUBLE);

    // cycles
    int nc = n_cycle;
    h5.writeAttribute("n_cycle", nc, H5T_NATIVE_INT);

    // steps
    int steps = level_steps;
    h5.writeAttribute("level_steps", steps, H5T_NATIVE_INT);

    // count
    int count = level_count;
    h5.writeAttribute("level_count", count, H5T_NATIVE_INT);

    // logical problem domain
    writeBoxOnHDF5(geom.Domain(), h5, "prob_domain");

    // real problem domain
    writeRealBoxOnHDF5(geom.ProbDomain(), h5, "real_domain");

    // box layout etc
    SortedGrids sg(grids, dmap);
    sg.sortedGrids.writeOnHDF5(h5, "boxes");

    // number of ghost cells (always zero in plotfiles)
    H5 attr = h5.createGroup("data_attributes");
    IntVect nGrow(0);
    int_h5_t gint = writeH5IntVec(nGrow.getVect());
    attr.writeAttribute("outputGhost", gint, intvect_id);
    H5Tclose(intvect_id);
    attr.closeGroup();
}

void WriteMultiLevelPlotfileHDF5 (const std::string& plotfilename, 
				  int nlevels,
                         	  const Vector<const MultiFab*>& mf,
                         	  const Vector<std::string>& varnames,
                         	  const Vector<Geometry>& geom, 
				  Real time, 
				  const Vector<int>& level_steps,
                         	  const Vector<IntVect>& ref_ratio,
                         	  const std::string &versionName,
                         	  const std::string &levelPrefix,
                         	  const std::string &mfPrefix,
                         	  const Vector<std::string>& extra_dirs)
{
    BL_PROFILE("WriteMultiLevelPlotfileHDF5");

    BL_ASSERT(nlevels <= mf.size());
    BL_ASSERT(nlevels <= geom.size());
    BL_ASSERT(nlevels <= ref_ratio.size()+1);
    BL_ASSERT(nlevels <= level_steps.size());
    BL_ASSERT(mf[0]->nComp() == varnames.size());

    std::string output_name = plotfilename + ".hdf5";

    H5 output_h5;
    output_h5.createFile(output_name, ParallelDescriptor::Communicator());

    std::map<std::string, int> vMInt;
    std::map<std::string, Real> vMReal;
    std::map<std::string, std::string> vMString;
        
    // generate and populate Chombo_global (required for visIt)
    vMInt["SpaceDim"] = amrex::SpaceDim;
    vMReal["testReal"] = 0.0;
    vMString["testString"] = "test string";
    H5 global = output_h5.createGroup("Chombo_global");
    global.writeAttribute(vMInt, vMReal, vMString);
    global.closeGroup();
        
    vMInt.clear();
    vMReal.clear();
    vMString.clear();
        
    // populate the root folder
    vMReal["time"] = time;
    vMInt["iteration"] = level_steps[0];
    vMInt["max_level"] = nlevels-1;
    vMInt["num_levels"] = nlevels-1;
        
    vMString["filetype"] = "VanillaAMRFileType";
    vMInt["num_levels"] = nlevels;
    vMInt["num_components"] = mf[0]->nComp();
    for (int ivar(0); ivar < mf[0]->nComp(); ++ivar) {
        char labelChSt[100];
        sprintf(labelChSt, "component_%d", ivar);
        std::string label(labelChSt);
        vMString[label] = varnames[ivar];
    }
    output_h5.writeAttribute(vMInt, vMReal, vMString);

    for (int level = 0; level <= nlevels; ++level)
    {
        // create a group for all of the levels data
        H5 level_grp = output_h5.createGroup("/level_" + num2str(level));
    
        writeLevelAttrHDF5(level_grp, level, geom[level], ref_ratio[level], time, 0.0,
                           0.0, 0, level_steps[level], 0, mf[level]->boxArray(),
                           mf[level]->DistributionMap());
        writeMultiFab(level_grp, mf[level], time);
        level_grp.closeGroup();
    }

    output_h5.closeFile();
}

void
WriteSingleLevelPlotfileHDF5 (const std::string& plotfilename,
                              const MultiFab& mf, const Vector<std::string>& varnames,
                              const Geometry& geom, Real time, int level_step,
                              const std::string &versionName,
                              const std::string &levelPrefix,
                              const std::string &mfPrefix,
                              const Vector<std::string>& extra_dirs)
{
    Vector<const MultiFab*> mfarr(1,&mf);
    Vector<Geometry> geomarr(1,geom);
    Vector<int> level_steps(1,level_step);
    Vector<IntVect> ref_ratio;
    
    WriteMultiLevelPlotfileHDF5(plotfilename, 1, mfarr, varnames, geomarr, time,
                                level_steps, ref_ratio, versionName, levelPrefix, mfPrefix, extra_dirs);
}
    


#endif
}
