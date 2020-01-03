#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <iostream>


#include "myfunc.H"
#include "myfunc_F.H"  // includes advance.cpp
#include "AMReX_SDCstruct.H"
int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    
    main_main();
    
    amrex::Finalize();
    return 0;
}

void main_main ()
{
    Real a;  // advection coef.
    Real d;  // diffusion coef.
    Real r;  // reaction coef. 
    Real test_norm;
    // AMREX_SPACEDIM: number of dimensions
    int n_cell, max_grid_size, Nsteps, plot_int, Nprob;
    Vector<int> bc_lo(AMREX_SPACEDIM,0);
    Vector<int> bc_hi(AMREX_SPACEDIM,0);

    // What time is it now?  We'll use this to compute total run time.
    Real strt_time = ParallelDescriptor::second();

    // ParmParse is way of reading inputs from the inputs file
    ParmParse pp;
    
    // We need to get n_cell from the inputs file - this is the number of cells on each side of 
    //   a square (or cubic) domain.
    pp.get("n_cell",n_cell);
    
    // The domain is broken into boxes of size max_grid_size
    pp.get("max_grid_size",max_grid_size);
    
    // Default plot_int to -1, allow us to set it to something else in the inputs file
    //  If plot_int < 0 then no plot files will be writtenq
    plot_int = -1;
    pp.query("plot_int",plot_int);

    // Set  plot_err = 1 to  output the error to plot files instead of the solution
    int plot_err = 1;
    
    // Read in number of steps and final time
    pp.query("Nsteps",Nsteps);
    Real Tfin=0.0;
    pp.query("Tfin",Tfin);    
    Real dt = Tfin/Nsteps;  // Set the time step

    // read in BC; see Src/Base/AMReX_BC_TYPES.H for supported types
    pp.queryarr("bc_lo", bc_lo);
    pp.queryarr("bc_hi", bc_hi);
    
    //  Read in the coefficients for A-D-R
    pp.query("a",a);
    pp.query("d",d);
    pp.query("r",r);
    pp.query("Nprob",Nprob);
    int Lord=222;
    pp.query("Lord",Lord);
    
    // Manufactured solution parameters
    //    Real k_freq =3.14159265358979323846;
    Real k_freq =1.0;
    pp.query("k_freq",k_freq);    
    Real epsilon = 0.25;//0.25;//0.25;
    Real kappa = 2.0*d*pow(k_freq,2.0); // This choice leads to cancellation analytically. Doesn't matter now.
    
    
    // Set BC based on Nprob:
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        if (Nprob<3) {
            bc_lo[idim] = 3;
            bc_hi[idim] = 3;
        }
        else if (Nprob>=3) {
            bc_lo[idim] = 0;
            bc_hi[idim] = 0;
        }
    }
    
    // determine whether boundary conditions are periodic
    Vector<int> is_periodic(AMREX_SPACEDIM,0);
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        if (bc_lo[idim] == INT_DIR && bc_hi[idim] == INT_DIR) {
            is_periodic[idim] = 1;
        }
        else {
            is_periodic[idim] = 0;  // Need this?
        }
    }
    // Check that we have external dirichlet conditions in place.
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        if (!(bc_lo[idim] == EXT_DIR && bc_hi[idim] == EXT_DIR)) {
            if (!(bc_lo[idim] == INT_DIR && bc_hi[idim] == INT_DIR)) {
                amrex::Abort("Not Dirichlet BC or periodic");
            }
        }
    }
    
    

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
        IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);
        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        ba.maxSize(max_grid_size);

       // This defines the physical box, [-1,1] in each direction.
        RealBox real_box({AMREX_D_DECL(-1.0,-1.0,-1.0)},
                         {AMREX_D_DECL( 1.0, 1.0, 1.0)});

        // This defines a Geometry object
        geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    }

    // Nghost = number of ghost cells for each array 
    int Nghost = 2;
    
    // Ncomp = number of components for each array
    int Ncomp  = 1;

    // time = starting time in the simulation
    Real time = 0.0;
  
    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // We allocate two phi multifabs; one will store the old state, the other the new.
    MultiFab phi_old(ba, dm, Ncomp, Nghost);
    MultiFab phi_new(ba, dm, Ncomp, Nghost);

    // Initialize phi_new by calling a Fortran routine (init_phi_2d.f90).
    // MFIter = MultiFab Iterator
    for ( MFIter mfi(phi_new); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.growntilebox();
        init_phi(BL_TO_FORTRAN_BOX(bx),
                 BL_TO_FORTRAN_ANYD(phi_new[mfi]),
                 geom.CellSize(), geom.ProbLo(), geom.ProbHi(), &k_freq, &Nprob);
    }
    amrex::Print() << "intial norm " << phi_new.norm0() << "\n";
    // Set up BCRec; see Src/Base/AMReX_BC_TYPES.H for supported types
    Vector<BCRec> bc(phi_old.nComp());
    for (int n = 0; n < phi_old.nComp(); ++n)
      {
	for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
	  {
	    // lo-side BCs
	    if (bc_lo[idim] == EXT_DIR) {
	      bc[n].setLo(idim, BCType::ext_dir);  // Dirichlet uses "external Dirichlet"
	    }
        else if (bc_lo[idim] == INT_DIR) {
            bc[n].setLo(idim, BCType::int_dir); // Periodic uses "internal Dirichlet"
        }
	    else {
	      amrex::Abort("Invalid bc_lo");
	    }
	    
	    // hi-side BCs
	    if (bc_hi[idim] == EXT_DIR) {
	      bc[n].setHi(idim, BCType::ext_dir);  // Dirichlet uses "external Dirichlet"
        }
        else if (bc_hi[idim] == INT_DIR) {
            bc[n].setHi(idim, BCType::int_dir);  // Periodic uses "internal Dirichlet"
        }
	    else {
	      amrex::Abort("Invalid bc_hi");
	    }
	  } 
      }

    // Build the flux multifabs
    std::array<MultiFab, AMREX_SPACEDIM> flux;
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        // flux(dir) has one component, zero ghost cells, and is nodal in direction dir
        BoxArray edge_ba = ba;
        edge_ba.surroundingNodes(dir);
        flux[dir].define(edge_ba, dm, 1, 2); // Made it have 2 ghost cells
    }

    // Make an SDC structure
    int Nnodes=5;  // Default to 8th order
    int Npieces=2; // Default is full MISDC
    int Nsweeps=2*Nnodes-2;  //  This will give highest formal accuracy for Lobatto nodes
    pp.get("Nnodes",Nnodes);
    pp.get("Npieces",Npieces);
    pp.get("Nsweeps",Nsweeps);  //  Uncomment to adjust Nsweeps          

    //  Build the structure
    SDCstruct SDCmats(Nnodes,Npieces,phi_old);
    SDCmats.Nsweeps =Nsweeps;  // Number of SDC sweeps per time step
    
    const Real* dx = geom.CellSize();
    
    // Write a plotfile of the initial data if plot_int > 0 (plot_int was defined in the inputs file)

    MultiFab::Copy(phi_old, phi_new, 0, 0, 1, 2);
    
    
    if (plot_int > 0)
      {
	if (plot_err == 1)  // Turn the solution into the error
	  {
	    MultiFab::Copy(phi_old, phi_new, 0, 0, 1, 2);
	    for ( MFIter mfi(phi_new); mfi.isValid(); ++mfi )
	      {
		const Box& bx = mfi.validbox();
		err_phi(BL_TO_FORTRAN_BOX(bx),
			BL_TO_FORTRAN_ANYD(phi_new[mfi]),
			geom.CellSize(), geom.ProbLo(), geom.ProbHi(),&a,&d,&r,&time, &epsilon,&k_freq, &kappa, &Nprob);
	      }
	  }
	amrex::Print() << "max error in phi " << phi_new.norm0() << "\n";
	int n = 0;
	const std::string& pltfile = amrex::Concatenate("plt",n,5);
	WriteSingleLevelPlotfile(pltfile, phi_new, {"phi"}, geom, time, 0);
	if (plot_err == 1)  // Put the solution back
	  MultiFab::Copy(phi_new, phi_old, 0, 0, 1, 2);	
      }

  // Set an assorment of solver and parallization options and parameters
  // see AMReX_MLLinOp.H for the defaults, accessors, and mutators
  LPInfo info;
  
  // Implicit solve using MLABecLaplacian class
  //  MLABecLaplacian mlabec({geom}, {ba}, {dm}, info);
  Kerrek mlabec({geom}, {ba}, {dm}, info);
  // order of stencil
  int linop_maxorder = 4; // Change this?
  mlabec.setMaxOrder(linop_maxorder);
  
  // build array of boundary conditions needed by MLABecLaplacian
  // see Src/Boundary/AMReX_LO_BCTYPES.H for supported types
  std::array<LinOpBCType,AMREX_SPACEDIM> mgbc_lo;
  std::array<LinOpBCType,AMREX_SPACEDIM> mgbc_hi;
  
  for (int n = 0; n < phi_old.nComp(); ++n) 
    {
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
	{
	  // lo-side BCs
	  if (bc[n].lo(idim) == BCType::ext_dir) {
	    mgbc_lo[idim] = LinOpBCType::Dirichlet;
	  }
      else if (bc[n].lo(idim) == BCType::int_dir) {
          mgbc_lo[idim] = LinOpBCType::Periodic;
      }
	  else {
	    amrex::Abort("Invalid bc_lo");
	  }
	  
	  // hi-side BCs
	  if (bc[n].hi(idim) == BCType::ext_dir) {
	    mgbc_hi[idim] = LinOpBCType::Dirichlet;
	  }
      else if (bc[n].hi(idim) == BCType::int_dir) {
          mgbc_hi[idim] = LinOpBCType::Periodic;
      }
	  else {
	    amrex::Abort("Invalid bc_hi");
	  }
	}
    }
  
  // tell the solver what the domain boundary conditions are
  mlabec.setDomainBC(mgbc_lo, mgbc_hi);
    
  // scaling factors
  Real ascalar = 1.0;
  Real bscalar = 1.0;
  mlabec.setScalars(ascalar, bscalar);
  // Set up coefficient matrices
  MultiFab acoef(ba, dm, 1, 0);
  
  // fill in the acoef MultiFab and load this into the solver
  acoef.setVal(1.0);
  mlabec.setACoeffs(0, acoef);
    /*
     // bcoef lives on faces so we make an array of face-centered MultiFabs
     //   then we will in face_bcoef MultiFabs and load them into the solver.
     std::array<MultiFab,AMREX_SPACEDIM> face_bcoef;
     for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
     {
     const BoxArray& bamg = amrex::convert(acoef.boxArray(),
     IntVect::TheDimensionVector(idim));
     face_bcoef[idim].define(bamg, acoef.DistributionMap(), 1, 0);
     face_bcoef[idim].setVal(1.0);
     }
     mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(face_bcoef));
     */
    
    // Cell centered Beta values
    MultiFab BccCoef(ba, dm, 1, 2);
    
    // Need to implement 4th order extrapolation; for now, we simply use that the use that this function is periodic.
    for ( MFIter mfi(BccCoef); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.growntilebox();
        init_beta(BL_TO_FORTRAN_BOX(bx),
                 BL_TO_FORTRAN_ANYD(BccCoef[mfi]),
                 geom.CellSize(), geom.ProbLo(), geom.ProbHi(),&epsilon, &k_freq, &Nprob);
    }
    BccCoef.FillBoundary();

    std::array<MultiFab,AMREX_SPACEDIM> face_bcoef;
   // Product Storage
    std::array<MultiFab,AMREX_SPACEDIM> prod_stor;
    
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    { // Assumes .boxArray() doesn't include ghostCells (validBoxArray command?)
        const BoxArray& bamg = amrex::convert(BccCoef.boxArray(),
                                              IntVect::TheDimensionVector(idim));
        face_bcoef[idim].define(bamg, BccCoef.DistributionMap(), 1, 2);  face_bcoef[idim].setBndry(0);
        prod_stor[idim].define(bamg, acoef.DistributionMap(), 1,2);
        // Apply fortran routine to find face valued b_coeffs using cell centered bcc.
        
        for ( MFIter mfi(BccCoef); mfi.isValid(); ++mfi )
        {          const Box& bx = mfi.validbox();
            cc_to_face_loc(BL_TO_FORTRAN_BOX(bx),
                       BL_TO_FORTRAN_ANYD(BccCoef[mfi]),
                       //  BL_TO_FORTRAN_BOX(bamg[mfi]),
                       BL_TO_FORTRAN_ANYD(face_bcoef[idim][mfi]),
                       &idim );
        }

        face_bcoef[idim].FillBoundary(geom.periodicity()); // Shouldn't need this?
    }

mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(face_bcoef));
    // Need to be able to find the boundary conditions at a given time. We do this intialization a quick and dirty way. Would be better to utilize maskvals and oitr.
    // Boundary values are stored in ghost cells in cross.
    MultiFab bdry_values(ba, dm, 1, 1); bdry_values.setBndry(0);
    
    for ( MFIter mfi(bdry_values); mfi.isValid(); ++mfi )
    {          const Box& bx = mfi.validbox();
        fill_bdry_values(BL_TO_FORTRAN_BOX(bx),
                         BL_TO_FORTRAN_ANYD(bdry_values[mfi]),
                         geom.CellSize(), geom.ProbLo(), geom.ProbHi(),&time, &epsilon,&k_freq, &kappa,&Nprob);
    }
/*    for ( MFIter mfi(bdry_values); mfi.isValid(); ++mfi )
    {          const Box& bx = mfi.validbox();
        print_multifab(BL_TO_FORTRAN_ANYD(phi_new[mfi]));
    }*/
    
    
    if (Nprob<3){ // Dirichlet case
    mlabec.fourthOrderBCFill(phi_new,bdry_values);
    }
    phi_new.FillBoundary(geom.periodicity());
    
    
  /*  for ( MFIter mfi(bdry_values); mfi.isValid(); ++mfi )
    {          const Box& bx = mfi.validbox();
        print_multifab(BL_TO_FORTRAN_ANYD(phi_new[mfi]));
    }*/
    
    // build an MLMG solver
    MLMG mlmg(mlabec);
  
  // set solver parameters
  int max_iter = 100;
  mlmg.setMaxIter(max_iter);
  int max_fmg_iter = 0;
  mlmg.setMaxFmgIter(max_fmg_iter);
  int verbose = 0;
  mlmg.setVerbose(verbose);
  int cg_verbose = 0;
  mlmg.setCGVerbose(cg_verbose);
  
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
 /*   /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // 1s and 0s TESTING SPACE /////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    MultiFab eval_storage(ba, dm, 1, 0);
    MultiFab rhs(ba, dm, 1, 0);
    MultiFab soln(ba, dm, 1, 2);
    MultiFab approx_soln(ba, dm, 1, 2);
    MultiFab bvalues(ba,dm,1,1);
    
    mlabec.setScalars(1.0, d);
    
    const Box& domain_bx = geom.Domain();
   // const Real* dx = geom.CellSize();
    Real zeroReal = 0.0;
    int oneInt = 1;
    
    
   /////////////////////////////////////////////////////
   // SET SOLN VALUE
   /////////////////////////////////////////////////////
    for ( MFIter mfi(phi_new); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();
        init_phi(BL_TO_FORTRAN_BOX(bx),
                 BL_TO_FORTRAN_ANYD(soln[mfi]),
                 geom.CellSize(), geom.ProbLo(), geom.ProbHi(), &k_freq);
    }
    
    
    
    ///////////////////////////////////////////////////////
    // Intialize Boundary Values
    ///////////////////////////////////////////////////////
    bvalues.setVal(0.0);
    
    for ( MFIter mfi(bdry_values); mfi.isValid(); ++mfi )
    {          const Box& bx = mfi.validbox();
        fill_bdry_values(BL_TO_FORTRAN_BOX(bx),
                         BL_TO_FORTRAN_ANYD(bvalues[mfi]),
                         geom.CellSize(), geom.ProbLo(), geom.ProbHi(),&time, &epsilon,&k_freq, &kappa);
    }
    
    
    
    /////////////////////////////////////////////////////
    // APPLY BC TO SOLN
    /////////////////////////////////////////////////////
    
    
    
    
    amrex::Print() << "Check 1"  << "\n";
    
    /////////////////////////////////////////////////////
    // Find L^2,4,4 eval of soln (make sure to comment out correct piece of feval)
    /////////////////////////////////////////////////////
    for ( MFIter mfi(soln); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();
        SDC_feval_F(BL_TO_FORTRAN_BOX(bx),
                    BL_TO_FORTRAN_BOX(domain_bx),
                    BL_TO_FORTRAN_ANYD(soln[mfi]),
                    BL_TO_FORTRAN_ANYD(flux[0][mfi]),
                    BL_TO_FORTRAN_ANYD(flux[1][mfi]),
#if (AMREX_SPACEDIM == 3)
                    BL_TO_FORTRAN_ANYD(flux[2][mfi]),
#endif
                    BL_TO_FORTRAN_ANYD(eval_storage[mfi]),
                    dx,&a,&d,&r,
                    BL_TO_FORTRAN_ANYD(face_bcoef[0][mfi]),
                    BL_TO_FORTRAN_ANYD(face_bcoef[1][mfi]),
                    BL_TO_FORTRAN_ANYD(prod_stor[0][mfi]),
                    BL_TO_FORTRAN_ANYD(prod_stor[1][mfi]),
                    &oneInt, &zeroReal, &zeroReal, &zeroReal, &zeroReal);
        
    }
    
    /////////////////////////////////////////////////////////////////
    // MAKE RHS = (I-L^2,4,4)soln
    /////////////////////////////////////////////////////////////////
    
    rhs.setVal(0.0);
    MultiFab::Saxpy(rhs,1.0,soln,0,0,1,0);
    MultiFab::Saxpy(rhs,-1.0,eval_storage,0,0,1,0);
    
    ////////////////////////////////////////////////////////////////
    // INITIAL GUESS FOR SOLN
    ////////////////////////////////////////////////////////////////
    
    approx_soln.setVal(0.0);
    
    ////////////////////////////////////////////////////////////////
    // SMOOTHING LOOP
    ////////////////////////////////////////////////////////////////
    for (int i =1;i<=10;i++){
        
    
        ///////////////////////////////////////////////////////////////
        // BC FOR APPROX
        ///////////////////////////////////////////////////////////////
      
        
        
        // SMOOTH
        mlabec.Fsmooth(0,0, approx_soln, rhs, 0);
    
    
    }
    
    
    
    
    ////////////////////////////////////////////////////
    // PRINT MULTIFAB to SEE RESULT
    ////////////////////////////////////////////////////
 //   for ( MFIter mfi(approx_soln); mfi.isValid(); ++mfi )
 //   {
 //       const Box& bx = mfi.validbox();
 //       print_multifab(BL_TO_FORTRAN_ANYD(approx_soln[mfi]));
 //   }
    
    
    
    //////////////////////////////////////////////////////
    // ERROR NORM
    //////////////////////////////////////////////////////
    
    MultiFab::Saxpy(approx_soln,-1.0,soln,0,0,1,0);
    amrex::Print() << "Error norm is " << approx_soln.norm0()  << "\n";
    
    
    // PAUSE AFTER RUNNING THIS AREA
    std::cin.get();
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // END 1s and 0s TESTING SPACE /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    */
   // mlabec.prepareForSolve();
    
    
  //  Do the time stepping
  // For now m_bcc assumed to just have spatial dependence. Time dependence should be easy to include.
  MultiFab::Copy(phi_old, phi_new, 0, 0, 1, 2);
  int totV=0;
  for (int n = 1; n <= Nsteps; ++n)
    {
      //amrex::Print() << "time" << time << "\n";
      // Do an SDC step
      
      SDC_advance(phi_old, phi_new,flux, dt, geom, bc, mlmg,mlabec,SDCmats,a,d,r ,face_bcoef, prod_stor,time, epsilon, k_freq, kappa, bdry_values, Nprob,Lord,totV);
       
      
      MultiFab::Copy(phi_old, phi_new, 0, 0, 1, 2);    
      time = time + dt;
      
      if (plot_err == 1)  // Turn the solution into the error
       // if(time >0.0999){
                for ( MFIter mfi(phi_new); mfi.isValid(); ++mfi )
                  {
                    const Box& bx = mfi.validbox();
                    err_phi(BL_TO_FORTRAN_BOX(bx),
                        BL_TO_FORTRAN_ANYD(phi_new[mfi]),
                        geom.CellSize(), geom.ProbLo(), geom.ProbHi(),&a,&d,&r,&time, &epsilon,&k_freq, &kappa, &Nprob);
                  }
            
                    // Tell the I/O Processor to write out which step we're doing
                  //  amrex::Print() << "Advanced step " << n << "\n";
                amrex::Print() << "max error in phi " << phi_new.norm0() << "\n";
                    // Write a plotfile of the current data (plot_int was defined in the inputs file)
                    if (plot_int > 0 && n%plot_int == 0)
                    {
                        const std::string& pltfile = amrex::Concatenate("plt",n,5);
                        WriteSingleLevelPlotfile(pltfile, phi_new, {"phi"}, geom, time, n);
                    }
            
       // }
    }
    // Call the timer again and compute the maximum difference between the start time and stop time
    //   over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time,IOProc);

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;
    amrex::Print() << "total V cycles = " << Nsteps << "  [ " << totV <<", "<< phi_new.norm0() <<"],\n";

}
