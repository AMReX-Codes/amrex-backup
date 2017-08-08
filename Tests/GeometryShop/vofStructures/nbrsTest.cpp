#include <cmath>
#include <type_traits>

#include "AMReX_EBIndexSpace.H"
#include "AMReX_EBLevelGrid.H"
#include "AMReX_GeometryShop.H"
#include "AMReX_BoxIterator.H"
#include "AMReX_ParmParse.H"
#include "AMReX_RealVect.H"
#include "AMReX_SphereIF.H"
#include "AMReX_FlatPlateGeom.H"
#include "AMReX_EBGraph.H"
#include "AMReX_EBISBox.H"
#include "AMReX_MultiFab.H"
#include "AMReX_iMultiFab.H"
#include "AMReX_BLProfiler.H"
#include "AMReX_BLFort.H"
#include "AMReX_VisMF.H"
#include "AMReX_PlotFileUtil.H"
#include "AMReX_MeshRefine.H"

#include "AMReX_ArrayLim.H"
#include "AMReX_EBRedist.H"

using namespace amrex;

static void
get_EBLG(EBLevelGrid& eblg)
{
    GridParameters params;
    getGridParameters(params, 0, true);

    ParmParse pp;

    EBIndexSpace* ebisPtr = AMReX_EBIS::instance();
    int which_geom = params.whichGeom;
    bool verbosity=true;
    pp.query("verbosity",verbosity);
    if (which_geom == 0)
    {
        RealVect origin = RealVect::Zero;
        Real radius = 0.5;
        pp.get("sphere_radius", radius);
        std::vector<Real> centervec(SpaceDim);
        pp.getarr("sphere_center", centervec, 0, SpaceDim);
        RealVect center;
        for(int idir = 0; idir < SpaceDim; idir++)
        {
            center[idir] = centervec[idir];
        }

        //bool insideRegular = false;
        bool insideRegular = true;
        SphereIF sphere(radius, center, insideRegular);

        GeometryShop gshop(sphere, verbosity);
        ebisPtr->define(params.coarsestDomain, origin, params.coarsestDx, gshop);
    }
    else if (which_geom == 1) 
    {
        RealVect origin = RealVect::Zero;
        std::vector<Real>  platelovec(SpaceDim);
        std::vector<Real>  platehivec(SpaceDim);

        int normalDir;
        Real plateLoc;
    
        pp.getarr("plate_lo", platelovec, 0, SpaceDim);
        pp.getarr("plate_hi", platehivec, 0, SpaceDim);
        pp.get("plate_location", plateLoc);
        pp.get("plate_normal", normalDir);

        RealVect plateLo, plateHi;
        for(int idir = 0; idir < SpaceDim; idir++)
        {
            plateLo[idir] = platelovec[idir];
            plateHi[idir] = platehivec[idir];
        }
        FlatPlateGeom flat_plate(normalDir, plateLoc, plateLo, plateHi);
        ebisPtr->define(params.coarsestDomain, origin, params.coarsestDx, flat_plate);
    }
    else
    {
        amrex::Abort("Unknown geom_type");
    }

    std::vector<EBLevelGrid> eblg_tmp;
    getAllIrregEBLG(eblg_tmp, params);
    eblg = eblg_tmp[0];

}

static void
Copy(intDIM& out, const IntVect& in)
{
    
    for (int d=0; d<3; ++d) {
	out[d] = 0;
    }
    for (int d=0; d<BL_SPACEDIM; ++d) {
        out[d] = in[d];
    }    
}


static void
Copy(IntVect& out, const intDIM& in)
{
    for (int d=0; d<BL_SPACEDIM; ++d) {
        out[d] = in[d];
    }    
}

extern "C"
{
    void fill_redist_stencil(const int*  lo, const int*  hi,
                             const int* slo, const int* shi,
			     NbrSten* redistSten, const int* Nsten,
			     BL_FORT_FAB_ARG_3D(vfrac));

    void fill_flux_interp_stencil(const int* lo,  const int*  hi,
                                  const int* slo, const int* shi,
                                  FaceSten* faceSten, const int* Nsten, const int* idir,
                                  const BL_FORT_FAB_ARG_3D(fcent));

    void apply_redist_stencil(const int* lo,  const int*  hi,
                              const int* slo, const int* shi,
                              NbrSten* redistSten, const int* Nsten,
			      const BL_FORT_FAB_ARG_3D(vin),
			      BL_FORT_FAB_ARG_3D(vout));

    void apply_flux_interp_stencil(const int* lo,  const int*  hi,
                                   const int* slo, const int* shi,
                                   FaceSten* redistSten, const int* Nsten, const int* idir,
                                   const BL_FORT_FAB_ARG_3D(vin),
                                   BL_FORT_FAB_ARG_3D(vout));
}

int myTest()
{
    EBLevelGrid eblg;
    get_EBLG(eblg);
    const BoxArray& ba = eblg.getBoxArray();
    const DistributionMapping& dm = eblg.getDM();
    const Box& domain = eblg.getDomain();
    Geometry geom(domain); 

    MultiFab vfrac(ba,dm,1,1);
    for(MFIter  mfi(vfrac); mfi.isValid(); ++mfi)
    {
        const EBISBox& ebis = eblg.getEBISL()[mfi];
        const Box& vbox = mfi.validbox();
	FArrayBox& vfab = vfrac[mfi];
	vfab.setVal(0);
	
        if (ebis.isAllRegular())
	{
	    vfab.setVal(1);
	}
	else
        {
	    for (IntVect iv(vbox.smallEnd()); iv<=vbox.bigEnd(); vbox.next(iv))
	    {
		const std::vector<VolIndex> vofs = ebis.getVoFs(iv);
		Real vtot = 0;
                BL_ASSERT(vofs.size()<=1);
		for (int j=0; j<vofs.size(); ++j)
		{
                    const VolIndex& vof = vofs[j];
		    vtot += ebis.volFrac(vof);
		}
		vfab(iv,0) = vtot;
	    }
	}
    }
    vfrac.FillBoundary();
    //Array<std::string> name(1,"vfrac");
    //WriteSingleLevelPlotfile("pltfile",vfrac,name,geom,0,0,"CartGrid-V2.0");

    MultiFab afrac[3];
    MultiFab fcent[3];
    for (int idir=0; idir<BL_SPACEDIM; ++idir)
    {
        afrac[idir].define(BoxArray(ba).surroundingNodes(idir),dm,1,1);
        fcent[idir].define(BoxArray(ba).surroundingNodes(idir),dm,BL_SPACEDIM,1);

        for(MFIter  mfi(afrac[idir]); mfi.isValid(); ++mfi)
        {
            const EBISBox& ebis = eblg.getEBISL()[mfi];
            const Box& vbox = ba[mfi.index()];
            FArrayBox& afab = afrac[idir][mfi];
            FArrayBox& cfab = fcent[idir][mfi];
            afab.setVal(0);
            cfab.setVal(0.5);
	
            if (ebis.isAllRegular())
            {
                afab.setVal(1);
            }
            else
            {
                for (IntVect iv(vbox.smallEnd()); iv<=vbox.bigEnd(); vbox.next(iv))
                {
                    std::vector<VolIndex> vofs = ebis.getVoFs(iv);
                    BL_ASSERT(vofs.size()<=1);
                    for (int j=0; j<vofs.size(); ++j)
                    {
                        bool do_hi_side = iv[idir] == vbox.bigEnd()[idir];
                        const VolIndex& vof = vofs[j];
                        for (SideIterator sit; sit.ok(); ++sit)
                        {
                            if (sit() == Side::Lo || do_hi_side)
                            {
                                Real atot = 0;
                                const std::vector<FaceIndex> faces = ebis.getFaces(vof,idir,sit());
                                BL_ASSERT(faces.size() <= 1);
                                for (int k=0; k<faces.size(); ++k)
                                {
                                    atot += ebis.areaFrac(faces[k]);
                                }

                                const int iside = sit() == Side::Lo ? 0 : 1;
                                IntVect iv_face = iv + iside*BASISV(idir);
                                afab(iv_face,0) = atot;
                                
                                if (faces.size()>0)
                                {
                                    const RealVect cent = ebis.centroid(faces[0]);
                                    for (int idir1=0; idir1<BL_SPACEDIM; ++idir1)
                                    {
                                        cfab(iv_face,idir1) = cent[idir1];
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        afrac[idir].FillBoundary();
    }


    amrex::Print() << "Building stencil\n"; 

    const Box stenBox(IntVect(D_DECL(-1,-1,-1)),
                      IntVect(D_DECL(+1,+1,+1)));
    intDIM idx3;

    FArrayBox stenMask(stenBox,1);
    std::map<int,std::vector<NbrSten>> redistSten;
    for(MFIter  mfi(vfrac); mfi.isValid(); ++mfi)
    {
        const EBISBox& ebis = eblg.getEBISL()[mfi];
        const Box& vbox = mfi.validbox();
	int i = mfi.index();
	FArrayBox& vfab = vfrac[mfi];
	
        if ( !(ebis.isAllRegular())) 
        {
            IntVectSet isIrreg = ebis.getIrregIVS(vbox);
            for (IVSIterator it(isIrreg); it.ok(); ++it)
            {
                const IntVect& iv=it();
                std::vector<VolIndex> vofs = ebis.getVoFs(iv);
                for (int ivof=0; ivof<vofs.size(); ++ivof)
                {

                    // Note here that if cells are multivalued, we can probably 
                    // stack these stencils here... an exercise for the reader...
		    redistSten[i].push_back(NbrSten());
		    NbrSten& nbr = redistSten[i].back();
		    Copy(nbr.iv,iv);
 
                    // set stenMask values to +1 if connected, -1 if not
                    stenMask.setVal(-1);
                    const VolIndex& vof = vofs[ivof];
                    stenMask(IntVect::TheZeroVector(),0) = +1;

                    for (int idir=0; idir<BL_SPACEDIM; ++idir)
                    {
                        for (SideIterator sit; sit.ok(); ++sit)
                        {
                            std::vector<VolIndex> nbrs = ebis.getVoFs(vof,idir,sit(),1);
                            for (int inbr=0; inbr<nbrs.size(); ++inbr)
                            {
                                const VolIndex& vof1 = nbrs[inbr];
                                stenMask(vof1.gridIndex() - iv,0) = +1;

                                for (int idir1=0; idir1<BL_SPACEDIM; ++idir1)
                                {
                                    if (idir != idir1)
                                    {
                                        for (SideIterator sit1; sit1.ok(); ++sit1)
                                        {
                                            std::vector<VolIndex> nnbrs = ebis.getVoFs(vof1,idir1,sit1(),1);
                                            for (int innbr=0; innbr<nnbrs.size(); ++innbr)
                                            {
                                                const VolIndex& vof2 = nnbrs[innbr];
                                                IntVect ivLoc = vof2.gridIndex() - iv;
                                                BL_ASSERT(stenBox.contains(ivLoc));
                                                stenMask(ivLoc,0) = +1;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    for (IntVect ivs=stenBox.smallEnd(); ivs<=stenBox.bigEnd(); stenBox.next(ivs))
                    {
                        Copy(idx3,ivs);
                        nbr.val[idx3[2]+1][idx3[1]+1][idx3[0]+1] = stenMask(ivs,0);
                    }
                }
            }

	    int Nsten = redistSten[i].size();
	    fill_redist_stencil(BL_TO_FORTRAN_BOX(vbox),
                                BL_TO_FORTRAN_BOX(stenBox),
                                redistSten[i].data(), &Nsten,
                                BL_TO_FORTRAN_3D(vfab));
        }
    }


    Box fbox[3];
    std::array<std::map<int,std::vector<FaceSten>>,3> faceSten;

    for (int idir=0; idir<BL_SPACEDIM; ++idir)
    {
        fbox[idir] = amrex::bdryLo(Box(IntVect(D_DECL(0,0,0)),
                                       IntVect(D_DECL(0,0,0))),idir,1);
        for (int idir1=0; idir1<BL_SPACEDIM; ++idir1)
        {
            if (idir1 != idir) fbox[idir].grow(idir1,1);
        }

        FArrayBox fsten(fbox[idir],1);

        for(MFIter  mfi(fcent[idir]); mfi.isValid(); ++mfi)
        {
            const EBISBox& ebis = eblg.getEBISL()[mfi];
            int i = mfi.index();
            const Box& vbox = ba[i];
            const Box& ebox = mfi.validbox();
            const FArrayBox& centroid = fcent[idir][mfi];

            std::set<IntVect> cut_faces;
            IntVectSet isIrreg = ebis.getIrregIVS(vbox);
            for (IVSIterator it(isIrreg); it.ok(); ++it)
            {
                const IntVect& iv=it();
                std::vector<VolIndex> vofs = ebis.getVoFs(iv);
                BL_ASSERT(vofs.size()<=1);
                for (int j=0; j<vofs.size(); ++j)
                {
                    const VolIndex& vof = vofs[j];
                    for (SideIterator sit; sit.ok(); ++sit)
                    {
                        const std::vector<FaceIndex> faces = ebis.getFaces(vof,idir,sit());
                        BL_ASSERT(faces.size() <= 1);
                        for (int k=0; k<faces.size(); ++k)
                        {
                            const FaceIndex& face = faces[k];
                            if (ebis.areaFrac(face) < 1) {
                                // IntVect of face associated with that of hi-side cell
                                cut_faces.insert(faces[k].gridIndex(Side::Hi));
                            }
                        }
                    }
                }
            }
            for (std::set<IntVect>::const_iterator it=cut_faces.begin(); it!=cut_faces.end(); ++it)
            {
                faceSten[idir][i].push_back(FaceSten());
                FaceSten& sten = faceSten[idir][i].back();
                Copy(sten.iv,*it);
            }

            int Nsten = faceSten[idir][i].size();
            fill_flux_interp_stencil(BL_TO_FORTRAN_BOX(ebox),
                                     BL_TO_FORTRAN_BOX(fbox[idir]),
                                     faceSten[idir][i].data(), &Nsten, &idir,
                                     BL_TO_FORTRAN_3D(centroid));
        }
    }
    

    amrex::Print() << "Applying stencil\n";
    MultiFab tin(ba,dm,1,1);
    tin.setVal(2);
    MultiFab tout(ba,dm,1,1);
    tout.setVal(0);

    for(MFIter  mfi(tin); mfi.isValid(); ++mfi)
    {
        const Box& vbox = mfi.validbox();
	int i = mfi.index();
	const FArrayBox& infab = tin[mfi];
	FArrayBox& outfab = tout[mfi];
	int Nsten = redistSten[i].size();
	
        if (Nsten > 0)
        {
	    apply_redist_stencil(BL_TO_FORTRAN_BOX(vbox),
                                 BL_TO_FORTRAN_BOX(stenBox),
                                 redistSten[i].data(), &Nsten,
                                 BL_TO_FORTRAN_3D(infab),
                                 BL_TO_FORTRAN_3D(outfab));
        }
    }

    for (int idir=0; idir<BL_SPACEDIM; ++idir)
    {
        amrex::Print() << "Applying face stencil " << idir << "\n";

        MultiFab fin(BoxArray(ba).surroundingNodes(idir),dm,1,1);
        fin.setVal(0);
        for(MFIter  mfi(fin); mfi.isValid(); ++mfi)
        {
            int i = mfi.index();
            FArrayBox& infab = fin[mfi];
            int Nsten = faceSten[idir][i].size();
            for (int n=0; n<Nsten; ++n)
            {
                IntVect iv;
                Copy(iv,faceSten[idir][i][n].iv);
                infab(iv,0) = 1;
            }
        }


        MultiFab fout(BoxArray(ba).surroundingNodes(idir),dm,1,1);
        fout.setVal(0);

        for(MFIter  mfi(tin); mfi.isValid(); ++mfi)
        {
            int i = mfi.index();
            const Box& ebox = mfi.validbox();
            const FArrayBox& infab = fin[mfi];
            FArrayBox& outfab = fout[mfi];
            int Nsten = faceSten[idir][i].size();
            
            if (Nsten > 0)
            {
                apply_flux_interp_stencil(BL_TO_FORTRAN_BOX(ebox),
                                          BL_TO_FORTRAN_BOX(fbox[idir]),
                                          faceSten[idir][i].data(), &Nsten, &idir,
                                          BL_TO_FORTRAN_3D(infab),
                                          BL_TO_FORTRAN_3D(outfab));
            }
        }
    }


    return 0;
}

int
main(int argc,char **argv)
{
    amrex::Initialize(argc,argv);
    {
        BL_PROFILE_VAR("main()", pmain);

        int retFlag = myTest();

        if (retFlag != 0)
        {
            amrex::Print() << "non zero return detected = " << retFlag << '\n';
            amrex::Print() << "sphere test failed" << '\n';
        }
        else
        {
            amrex::Print() << "test passed" << '\n';
        }

        BL_PROFILE_VAR_STOP(pmain);
    }
    amrex::Finalize();
    return 0;
}
