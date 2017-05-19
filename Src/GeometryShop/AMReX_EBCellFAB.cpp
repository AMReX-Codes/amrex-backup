/*
 *       {_       {__       {__{_______              {__      {__
 *      {_ __     {_ {__   {___{__    {__             {__   {__  
 *     {_  {__    {__ {__ { {__{__    {__     {__      {__ {__   
 *    {__   {__   {__  {__  {__{_ {__       {_   {__     {__     
 *   {______ {__  {__   {_  {__{__  {__    {_____ {__  {__ {__   
 *  {__       {__ {__       {__{__    {__  {_         {__   {__  
 * {__         {__{__       {__{__      {__  {____   {__      {__
 *
 */

#include "AMReX_EBCellFAB.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_BoxIterator.H"

namespace amrex
{
  /**********************/
  EBCellFAB::EBCellFAB():BaseEBCellFAB<Real>()
  {
  }
  /**********************/
  void
  EBCellFAB::define(const EBISBox& a_ebisBox,
                    const Box& a_region, int a_nVar)
  {
    BaseEBCellFAB<Real>::define(a_ebisBox, a_region, a_nVar);
  }
  /**********************/
  EBCellFAB::EBCellFAB(const EBISBox& a_ebisBox,
                       const Box& a_region, int a_nComp)
    :BaseEBCellFAB<Real>(a_ebisBox, a_region, a_nComp)
  {
  }
  /**********************/
  EBCellFAB::~EBCellFAB()
  {
  }
  /**********************/
  const FArrayBox&
  EBCellFAB::getFArrayBox() const
  {
    BL_ASSERT(isDefined());
    return (const FArrayBox&)m_regFAB;
  }
  /**********************/
  FArrayBox&
  EBCellFAB::getFArrayBox()
  {
    BL_ASSERT(isDefined());
    return (FArrayBox&)m_regFAB;
  }
  /**********************/
  EBCellFAB&
  EBCellFAB::negate(void)
  {
    BL_ASSERT(isDefined());
    (*this) *= -1.0;
    return *this;
  }

  /*******************************************************************************/
  // Routines to provide a consistent interface with FArrayBox
   FABio::Format
   EBCellFAB::getFormat ()
   {
       return getFArrayBox().getFormat();
   }

   void
   FArrayBox::writeOn (std::ostream& os) const
   {
       getFArrayBox().writeOn(os);
   }

   void
   FArrayBox::readFrom (std::ostream& os) const
   {
       getFArrayBox().readFrom(os);
   }

   void
   FArrayBox::readFrom (std::ostream& os, int compIndex) const
   {
       getFArrayBox().readFrom(os, int compIndex);
   }

/*
   *
   * Rework of minus to use explicitly specified src and dest boxes,
   * call BaseFab::minus for regular data. Still use BaseIVFAB<Real>::forall for irregular 
   * data
   *
   */

  EBCellFAB&
  EBCellFAB::minus(const EBCellFAB& a_src,
                   const Box&        srcbox,
                   const Box&        destbox,
                   int               a_srccomp,
                   int               a_dstcomp,
                   int               a_numcomp)
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(a_src.isDefined());
    BL_ASSERT(a_srccomp + a_numcomp <= a_src.nComp());
    BL_ASSERT(a_dstcomp + a_numcomp <= nComp());

    // Regular data
    m_regFAB.minus(a_src.m_regFAB, srcbox, destbox, a_srccomp, a_dstcomp, a_numcomp);

    // Irregular data - if both irrFABs do not contain all of the VOFS in srcbox, and that
    // box is not the same as the interection between the destbox and the box for the source 
    // we have a problem. Unlike regular data it isn't so easy to just rework the indexes to
    // map irregular data with a different index offset into the destination space.
         
    Box locRegion = ((a_src.getRegion() & getRegion()) & destbox) & srcbox;

    // Optimization - if regular boxes are the same as the destination box, 
    // the BaseIVFAB<Real>::forall will skip look up and operate on l[i], r[i]
    bool sameRegBox = ((a_src.m_regFAB.box() == m_regFAB.box() ) == destbox);
         
    if (!locRegion.isEmpty())
    {
      m_irrFAB.forall(a_src.m_irrFAB, locRegion, a_srccomp, a_dstcomp,
                      a_numcomp, sameRegBox, [](Real& dest, const Real& src){dest-=src;});
    }
    return *this;
  }
  
  EBCellFAB&
  EBCellFAB::minus(const EBCellFAB& a_src,
                   const Box&        bx,
                   int               a_srccomp,
                   int               a_dstcomp,
                   int               a_numcomp)
  {
    return minus( a_src, bx, bx, a_srccomp, a_dstcomp, a_numcomp);
   
  }
  
  EBCellFAB&
  EBCellFAB::invert (Real val
                    const Box& bx,
                    int        comp,
                    int        ncomp)
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(a_src.isDefined());
    BL_ASSERT(a_srccomp + a_numcomp <= a_src.nComp());
    BL_ASSERT(a_dstcomp + a_numcomp <= nComp());

    // Regular data
    m_regFAB.invert(val, bx, comp, ncomp);

    // Irregular data - if both irrFABs do not contain all of the VOFS in srcbox, and that
    // box is not the same as the interection between the destbox and the box for the source 
    // we have a problem. Unlike regular data it isn't so easy to just rework the indexes to
    // map irregular data with a different index offset into the destination space.
         
    Box locRegion =  getRegion() & bx;

    // Optimization - if regular boxes are the same as the destination box, 
    // the BaseIVFAB<Real>::forall will skip look up and operate on l[i], r[i]
    bool sameRegBox = ((a_src.m_regFAB.box() == m_regFAB.box() ) == destbox);

    if (!locRegion.isEmpty())
    {
      std::vector<VolIndex>  irrvofs = m_irrFAB.getVoFs();
      for(int ivof = 0; ivof < irrvofs.size(); ivof++)
      {
        const VolIndex& vof = irrvofs[ivof];
        if(locRegion.contains(vof.gridIndex()))
        {
          for (int icomp = comp; icomp < (comp+ncomp); ++icomp)
          {
            m_irrFAB(vof, icomp) = a / m_irrFAB(vof, icomp);
          }
        }
      }
      return *this;
    }
  }

  EBCellFAB&
  EBCellFAB::negate (const Box& bx,
                     int        comp,
                     int        ncomp)
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(a_src.isDefined());
    BL_ASSERT(a_srccomp + a_numcomp <= a_src.nComp());
    BL_ASSERT(a_dstcomp + a_numcomp <= nComp());

    // Regular data
    m_regFAB.negate(bx, comp, ncomp);

    // Irregular data - if both irrFABs do not contain all of the VOFS in srcbox, and that
    // box is not the same as the interection between the destbox and the box for the source 
    // we have a problem. Unlike regular data it isn't so easy to just rework the indexes to
    // map irregular data with a different index offset into the destination space.
         
    Box locRegion =  getRegion() & bx;

    // Optimization - if regular boxes are the same as the destination box, 
    // the BaseIVFAB<Real>::forall will skip look up and operate on l[i], r[i]
    bool sameRegBox = ((a_src.m_regFAB.box() == m_regFAB.box() ) == destbox);

    if (!locRegion.isEmpty())
    {
      std::vector<VolIndex>  irrvofs = m_irrFAB.getVoFs();
      for(int ivof = 0; ivof < irrvofs.size(); ivof++)
      {
        const VolIndex& vof = irrvofs[ivof];
        if(locRegion.contains(vof.gridIndex()))
        {
          for (int icomp = comp; icomp < (comp+ncomp); ++icomp)
          {
            m_irrFAB(vof, icomp) = -m_irrFAB(vof, icomp);
          }
        }
      }
      return *this;
    }
  }

  /*
   *
   * Scalar addition (a[i] <- a[i] + r) within a Box
   *
   */

  void
  EBCellFAB::plus (Real r, const Box& bx, int comp, int ncomp)
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(comp + ncomp <= nComp());


    m_regFAB.plus(r, bx, comp, ncomp);

    Box locRegion = getRegion() & bx;
    bool sameRegBox = (m_regFAB.box()  == bx);

    std::vector<VolIndex>  irrvofs = m_irrFAB.getVoFs();
    for(int ivof = 0; ivof < irrvofs.size(); ivof++)
      {
        const VolIndex& vof = irrvofs[ivof];
        if(locRegion.contains(vof.gridIndex()))
        {
          for (int icomp = comp; icomp < (comp+ncomp); ++icomp)
          {
            m_irrFAB(vof, icomp) += r;
          }
        }
      }
  }

  EBCellFAB::plus (Real r, const Box& bx)
  {
    plus(r, bx, 0, 1);
  }
  

  void
  EBCellFAB::mult (Real r, const Box& bx, int comp, int ncomp)
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(comp + ncomp <= nComp());


    m_regFAB.mult(r, bx, comp, ncomp);

    Box locRegion = getRegion() & bx;
    bool sameRegBox = (m_regFAB.box()  == bx);

    std::vector<VolIndex>  irrvofs = m_irrFAB.getVoFs();
    for(int ivof = 0; ivof < irrvofs.size(); ivof++)
      {
        const VolIndex& vof = irrvofs[ivof];
        if(locRegion.contains(vof.gridIndex()))
        {
          for (int icomp = comp; icomp < (comp+ncomp); ++icomp)
          {
            m_irrFAB(vof, icomp) *= r;
          }
        }
      }
  }

  

/**
   *
   * Rework of multiply to use explicitly specified src and dest boxes
   * call BaseFab::mult for regular data; still use BaseIVFab<Real>::foraall for irregular
   * data
   *
   */
  EBCellFAB&
  EBCellFAB::mult(const EBCellFAB& a_src,
                  const Box&       srcbox,
                  const Box&       destbox,
                  int a_srccomp,
                  int a_dstcomp,
                  int a_numcomp)
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(a_src.isDefined());
    BL_ASSERT(a_srccomp + a_numcomp <= a_src.nComp());
    BL_ASSERT(a_dstcomp + a_numcomp <= nComp());


    m_regFAB.mult(a_src.m_regFAB, srcbox, destbox, a_srccomp, a_dstcomp, a_numcomp);

    Box locRegion = ((a_src.getRegion() & getRegion()) & destbox) & srcbox;

    bool sameRegBox = (a_src.m_regFAB.box() == m_regFAB.box());
         
    if (!locRegion.isEmpty())
    {
      m_irrFAB.forall(a_src.m_irrFAB, locRegion, a_srccomp, a_dstcomp,
                      a_numcomp, sameRegBox, [](Real& dest, const Real& src){dest*=src;});
    }
         
    return *this;
  }

 /**
   *
   * Rework of divide to use explicitly specified src and dest boxes
   * call BaseFab::divide for regular data; still use BaseIVFab<Real>::foraall for irregular
   * data
   *
   */
  EBCellFAB&
  EBCellFAB::divide(const EBCellFAB& a_src,
                    const Box&       srcbox,
                    const Box&       destbox,
                    int              a_srccomp,
                    int              a_dstcomp,
                    int              a_numcomp)
  {
         
    BL_ASSERT(isDefined());
    BL_ASSERT(a_src.isDefined());
         
    BL_ASSERT(a_srccomp + a_numcomp <= a_src.nComp());
    BL_ASSERT(a_dstcomp + a_numcomp <= nComp());

    // regular data
    m_regFAB.divide(a_src.m_regFAB, srcbox, destbox,
                    a_srccomp, a_dstcomp, a_numcomp);
         
    // region for irregular data - only include VOFs in this box
    Box locRegion = ((a_src.getRegion() & getRegion()) & destbox) & srcbox;

    // Optimization to skip walking through VOFs
    bool sameRegBox = (a_src.m_regFAB.box() == m_regFAB.box());

    if (!locRegion.isEmpty())
    {
      m_irrFAB.forall(a_src.m_irrFAB, locRegion, a_srccomp, a_dstcomp,
                      a_numcomp, sameRegBox, [](Real& dest,
                      const Real& src){dest/=src;});
    }
    return *this;
  }
  
  
  EBCellFAB&
  EBCellFAB::divide(const EBCellFAB& a_src,
                    const Box&       bx,
                    int              a_srccomp,
                    int              a_dstcomp,
                    int              a_numcomp)
  {
      return divide(a_src, bx, bx, a_srccomp, a_dstcomp, a_numcomp);
  }

  /**
   * 
   * this = this + a*src
   * optimization through specialization.  Could do as this += a*src but
   * this version uses the same loop and may be able to use FMA as well as reduce
   * memory traffic
   *
   */
  EBCellFAB&
  EBCellFAB::
  saxpy(const Real& a_A,
        const EBCellFAB& src,
        const Box&       srcbox,
        const Box&       destbox,
        int              srccomp,
        int              dstcomp,
        int              numcomp)

{
    BL_ASSERT(isDefined());
    BL_ASSERT(a_src.isDefined());
         
    BL_ASSERT(a_srccomp + a_numcomp <= a_src.nComp());
    BL_ASSERT(a_dstcomp + a_numcomp <= nComp());

    // Regular data
    m_regFAB.saxpy( a_A, src.m_regFAB, srcbox, destbox, srccomp, dstcomp, numcomp);

   // region for irregular data - only include VOFs in this box
    Box locRegion = ((a_src.getRegion() & getRegion()) & destbox) & srcbox;

    // Optimization to skip walking through VOFs
    bool sameRegBox = (a_src.m_regFAB.box() == m_regFAB.box());

    if (!locRegion.isEmpty())
    {
      m_irrFAB.forall(a_src.m_irrFAB, locRegion, a_srccomp, a_dstcomp,
                      a_numcomp, sameRegBox, [](Real& dest,
                      const Real& src){dest+=a_A*src;});
    }
 
    return *this;
  }
  
    /**
   * 
   * this = a*this src
   * optimization through specialization.  Could do as this *= a and += src but
   * this version uses the same loop and may be able to use FMA as well as reduce
   * memory traffic
   *
   */
  EBCellFAB&
  EBCellFAB::
  xpay(const Real& a_A,
        const EBCellFAB& src,
        const Box&       srcbox,
        const Box&       destbox,
        int              srccomp,
        int              dstcomp,
        int              numcomp)

{
    BL_ASSERT(isDefined());
    BL_ASSERT(a_src.isDefined());
         
    BL_ASSERT(a_srccomp + a_numcomp <= a_src.nComp());
    BL_ASSERT(a_dstcomp + a_numcomp <= nComp());

    // Regular data
    m_regFAB.xpay( a_A, src.m_regFAB, srcbox, destbox, srccomp, dstcomp, numcomp);

   // region for irregular data - only include VOFs in this box
    Box locRegion = ((a_src.getRegion() & getRegion()) & destbox) & srcbox;

    // Optimization to skip walking through VOFs
    bool sameRegBox = (a_src.m_regFAB.box() == m_regFAB.box());

    if (!locRegion.isEmpty())
    {
      m_irrFAB.forall(a_src.m_irrFAB, locRegion, a_srccomp, a_dstcomp,
                      a_numcomp, sameRegBox, [](Real& dest,
                      const Real& src){dest=a_A*dest+src;});
    }
 
    return *this;
  }

/**
 *
 * this = a_A*X + a_B*Y
 *
 */
EBCellFAB&
EBCellFAB::linComb (const EBCellFAB& f1,
        			const Box&       b1,
        			int              comp1,
        			const EBCellFAB& f2,
        			const Box&       b2,
        			int              comp2,
        			Real             alpha,
        			Real             beta,
        			const Box&       b,
        			int              comp,
        			int              numcomp)
			
  {
    // Regular data
    m_regFAB.linComb (f1, b1,comp1,
			          f2,b2, comp2,
                      alpha, beta,
                        b,comp, numcomp);

    // region for irregular data - only include VOFs in this box
    Box locRegion = ((f1.getRegion() & f2.getRegion()) & b1) & b2 & b;
    
    std::vector<VolIndex>  irrvofs = m_irrFAB.getVoFs();
    for(int ivof = 0; ivof < irrvofs.size(); ivof++)
    {
        const VolIndex& vof = irrvofs[ivof];
        if(locRegion.contains(vof.gridIndex())){
            for (int icomp = comp; icomp < numcomp; icomp++)
            {
                m_irrFAB(vof, icomp) = alpha*f1.m_irrFAB(vof, icomp) 
                                       + beta*f2.m_irrFAB(vof, icomp);
            }
        }
    }
    return *this;
  }

    /*
     *
     * this += src1*src2
     *
     */

    EBCellFAB& addproduct (const Box&           destbox,
			               int                  destcomp,
			               int                  numcomp,
			               const EBCellFAB&     src1,
			               int                  comp1,
			               const EBCellFAB&     src2,
			               int                  comp2)

{
    // Regular data
    m_regFAB.addproduct (destbox, destcomp,
                         numcomp,
			             src1, comp1,
			             src2, comp2);
        

    // region for irregular data - only include VOFs in this box
    Box locRegion = ((getRegion() & src1.getRegion()) 
                      & src2.getRegion) & destbox;

    std::vector<VolIndex>  irrvofs = m_irrFAB.getVoFs();
    for(int ivof = 0; ivof < irrvofs.size(); ivof++)
    {
        const VolIndex& vof = irrvofs[ivof];
        if(locRegion.contains(vof.gridIndex())){
            for (int icomp = 0; icomp < numcomp; icomp++)
            {
                m_irrFAB(vof, destcomp+icomp) = src1.m_irrFAB(vof, comp1+icomp) 
                    + src2.m_irrFAB(vof, comp1+icomp);
            }
        }
    }
    return *this;
}


/**
 *
 * Dot product of x (i.e.,this) and y
 *
 */
  Real dot (const Box& xbx, int xcomp,
	   const EBCellFAB& y, const Box& ybx, int ycomp,
	   int numcomp) const
  {
    BL_ASSERT(isDefined());
    Real val = 0.0;
         
    // Dot product on intersection of two boxes
    const EBISBox& ebbox = getEBISBox();
    const IntVectSet validCells(xbx & ybx
    );
    for (VoFIterator vit(validCells, ebbox.getEBGraph()); vit.ok(); ++vit)
    {
      VolIndex vofi = vit();
      for (int icomp=0; icomp<numcomp; ++icomp)
      {
          val += (xbx(vofi, xcomp+icomp)*ybx(vofi,ycomp+icomp));
      } 

    }
    return val;
  }
}   

/**
 *
 * Set the initial value - either to initval or NaN
 *
 */
void
EBCellFAB::initVal ()
{
  BL_ASSERT(do_initval); // init as NaN not implemented
  setVal(initval);  
}

bool
EBCellFAB::set_do_initval (bool tf)
{
    bool o_tf = do_initval;
    do_initval = tf;
    return o_tf;
}

bool
EBCellFAB::get_do_initval ()
{
    return do_initval;
}

Real
EBCellFAB::set_initval (Real iv)
{
    Real o_iv = initval;
    initval = iv;
    return o_iv;
}

Real
EBCellFAB::get_initval ()
{
    return initval;
}

// Done by looping over vofs to exclude covered cells in the regular data
// from the check 
bool
EBCellFAB::contains_nan(const Box& bx, int scomp, int ncomp) const
{
    BL_ASSERT(isDefined());

    const IntVectSet validCells(bx);
    for (VoFIterator vit(validCells, ebbox.getEBGraph()); vit.ok(); ++vit)
    {
      VolIndex vofi = vit();
      for (int icomp = scomp; icomp < ncomp; icomp++)
       {
          if (std::isnan( (*this)(vofi, icomp)))
              return true;
        }
    }
    return false;
  }
}


// Done by looping over vofs to exclude covered cells in the regular data
// from the check 
bool
EBCellFAB::contains_inf(const Box& bx, int scomp, int ncomp) const
{
    BL_ASSERT(isDefined());

    const IntVectSet validCells(bx);
    for (VoFIterator vit(validCells, ebbox.getEBGraph()); vit.ok(); ++vit)
    {
      VolIndex vofi = vit();
      for (int icomp = scomp; icomp < ncomp; icomp++)
       {
          if (std::isinf( (*this)(vofi, icomp)))
              return true;
        }
    }
    return false;

}

// Done by looping over vofs to exclude covered cells in the regular data
// from the check 
IntVect
EBCellFAB::minIndex(const Box& bx, int comp) const
{
    BL_ASSERT(isDefined());
    Real val = 1.0e30;
    IntVect _min_loc(bx.smallEnd());

    const IntVectSet validCells(bx);
    for (VoFIterator vit(validCells, ebbox.getEBGraph()); vit.ok(); ++vit)
    {
        VolIndex vofi = vit();
        if ((*this)(vofi, a_comp) < val){
            val = (*this)(vofi, a_comp);
            _min_loc = vofi.gridIndex();
        }

    }
    return _min_loc;

}


// Done by looping over vofs to exclude covered cells in the regular data
// from the check 
IntVect
EBCellFAB::maxIndex(const Box& bx, int comp) const
{
    BL_ASSERT(isDefined());
    Real val = -1.0e30;
    IntVect _max_loc(bx.smallEnd());

    const IntVectSet validCells(bx);
    for (VoFIterator vit(validCells, ebbox.getEBGraph()); vit.ok(); ++vit)
    {
        VolIndex vofi = vit();
        if ((*this)(vofi, a_comp) > val){
            val = (*this)(vofi, a_comp);
            _max_loc = vofi.gridIndex();
        }

    }
    return _max_loc;

}

// Done by looping over vofs to exclude covered cells in the regular data
// from the check 
Real 
EBCellFAB::norm (const Box& bx, int p, int comp, int ncomp) const
{
    BL_ASSERT(isDefined());
    BL_ASSERT(comp >= 0 && comp + ncomp <= nvar);

    Real nrm;

    const IntVectSet validCells(bx);
    nrm = 0.0;
    if (p == 0 || p == 1)
    {
        for (VoFIterator vit(validCells, ebbox.getEBGraph()); vit.ok(); ++vit)
        {
            VolIndex vofi = vit();
            for (int icomp = comp; icomp < (comp+ncomp); ++icomp)
            {
                nrm = max(nrm, abs((*this)(vofi, icomp)));
            }
        }
    }
    else if (p ==1)
    {
        for (VoFIterator vit(validCells, ebbox.getEBGraph()); vit.ok(); ++vit)
        {
            VolIndex vofi = vit();
            for (int icomp = comp; icomp < (comp+ncomp); ++icomp)
            {
                nrm += abs((*this)(vofi, icomp));
            }
        }

    }
    else
    {
        amrex::Error("EBCellFAB::norm(): only p == 0 or p == 1 are supported");
    }
    return nrm;


}

// Done by looping over vofs to exclude covered cells in the regular data
// from the check 
Real 
EBCellFAB::sum (const Box& bx, int comp, int ncomp) const
{
    BL_ASSERT(isDefined());
    BL_ASSERT(comp >= 0 && comp + ncomp <= nvar);

    Real sum;

    const IntVectSet validCells(bx);
    sum = 0.0;

    for (VoFIterator vit(validCells, ebbox.getEBGraph()); vit.ok(); ++vit)
    {
        VolIndex vofi = vit();
        for (int icomp = comp; icomp < (comp+ncomp); ++icomp)
        {
            sum += abs((*this)(vofi, icomp));
        }
    }
  
    return sum;


}


  // End of routines to provide a consistent interface with FArrayBox
  /*******************************************************************************/
  /**********************/
  EBCellFAB&
  EBCellFAB::operator+=(const EBCellFAB& a_src)
  {
    BL_ASSERT(a_src.nComp() == nComp());
         
    plus(a_src, 0, 0, nComp());
         
    return *this;
  }
         
         
/**
 *
 * this = a_A*X + a_B*Y
 *
 */
  EBCellFAB&
  EBCellFAB::
  axby(const EBCellFAB& a_X, 
       const EBCellFAB& a_Y,
       const Real& a_A, const Real& a_B)
  {
    m_regFAB.linComb( a_X.m_regFAB, a_X.m_regFAB.box(), 0,
                      a_Y.m_regFAB, a_Y.m_regFAB.box(), 0,
                      a_A, b_B,
                      m_regFAB.box(), 0, nComp() )

    std::vector<VolIndex>  irrvofs = m_irrFAB.getVoFs();
    for(int ivof = 0; ivof < irrvofs.size(); ivof++)
    {
        for (int icomp = 0; icomp < nComp(); icomp++)
        {
          const VolIndex& vof = irrvofs[ivof];
          m_irrFAB(vof, icomp) = a_A*a_X.m_irrFAB(vof, icomp) + a_B*a_Y.m_irrFAB(vof, icomp);
        }
    }

    return *this;
  }

  void
  EBCellFAB::
  kappaWeight()
  {
    std::vector<VolIndex>  irrvofs = m_irrFAB.getVoFs();
    for(int ivof = 0; ivof < irrvofs.size(); ivof++)
    {
        for (int icomp = 0; icomp < nComp(); icomp++)
        {
          const VolIndex& vof = irrvofs[ivof];
          Real kappa = m_ebisBox.volFrac(vof);
          m_irrFAB(vof, icomp) *= kappa; 
        }
    }
  }

  EBCellFAB&
  EBCellFAB::plus(const EBCellFAB& a_src,
                  int a_srccomp,
                  int a_dstcomp,
                  int a_numcomp,
                  Real a_scale)
  {
    Box locRegion = a_src.getRegion() & getRegion();
    plus(a_src, locRegion, a_srccomp, a_dstcomp, a_numcomp, a_scale);
    return *this;
  }
         
  EBCellFAB& EBCellFAB::plus(const EBCellFAB& a_src,
                             const Box& a_region,
                             int  a_srccomp,
                             int  a_dstcomp,
                             int  a_numcomp,
                             Real a_scale)
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(a_src.isDefined());
    BL_ASSERT(a_srccomp + a_numcomp <= a_src.nComp());
    BL_ASSERT(a_dstcomp + a_numcomp <= nComp());
    const Box& locRegion = a_region;
         
    if (!locRegion.isEmpty())
    {
      Box region = locRegion & m_region;
      region &= a_src.m_region;
      for (BoxIterator boxit(region); boxit.ok(); ++boxit)
      {
        for (int icomp = 0; icomp < a_numcomp; icomp++)
        {
          int srcvar = a_srccomp + icomp;
          int dstvar = a_dstcomp + icomp;
          m_regFAB(boxit(), dstvar) += a_scale*a_src.m_regFAB(boxit(), srcvar);
        }
      }

      std::vector<VolIndex>  irrvofs = m_irrFAB.getVoFs();
      for(int ivof = 0; ivof < irrvofs.size(); ivof++)
      {
        const VolIndex& vof = irrvofs[ivof];
        if(locRegion.contains(vof.gridIndex()))
        {
          for (int icomp = 0; icomp < nComp(); icomp++)
          {
            m_irrFAB(vof, icomp) = a_scale*a_src.m_irrFAB(vof, icomp);
          }
        }
      }
    }
    return *this;
  }

  /**********************/
  void EBCellFAB::clone(const EBCellFAB& a_arg)
  {
    define(a_arg.m_ebisBox, a_arg.getRegion(), a_arg.nComp());
    copy(a_arg);
  }
  /**********************/
  EBCellFAB&
  EBCellFAB::operator-=(const EBCellFAB& a_src)
  {
    BL_ASSERT(a_src.nComp() == nComp());
         
    minus(a_src, 0, 0, nComp());
         
    return *this;
  }
  /**********************/
  EBCellFAB&
  EBCellFAB::minus(const EBCellFAB& a_src,
                   int a_srccomp,
                   int a_dstcomp,
                   int a_numcomp)
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(a_src.isDefined());
    BL_ASSERT(a_srccomp + a_numcomp <= a_src.nComp());
    BL_ASSERT(a_dstcomp + a_numcomp <= nComp());
         
    Box locRegion = a_src.getRegion() & getRegion();
    bool sameRegBox = (a_src.m_regFAB.box() == m_regFAB.box());
         
    if (!locRegion.isEmpty())
    {
      Box region = locRegion & m_region;
      region &= a_src.m_region;
      for (BoxIterator boxit(region); boxit.ok(); ++boxit)
      {
        for (int icomp = 0; icomp < a_numcomp; icomp++)
        {
          int srcvar = a_srccomp + icomp;
          int dstvar = a_dstcomp + icomp;
          m_regFAB(boxit(), dstvar) -= a_src.m_regFAB(boxit(), srcvar);
        }
      }
      m_irrFAB.forall(a_src.m_irrFAB, locRegion, a_srccomp, a_dstcomp, a_numcomp, sameRegBox, [](Real& dest, const Real& src){dest-=src;});
    }
    return *this;
  }
  

  
  
  /**********************/
  EBCellFAB&
  EBCellFAB::operator*=(const EBCellFAB& a_src)
  {
    BL_ASSERT(a_src.nComp() == nComp());
         
    mult(a_src, 0, 0, nComp());
         
    return *this;
  }
         
  EBCellFAB&
  EBCellFAB::mult(const EBCellFAB& a_src,
                  int a_srccomp,
                  int a_dstcomp,
                  int a_numcomp)
  {
    BL_ASSERT(isDefined());
    BL_ASSERT(a_src.isDefined());
    BL_ASSERT(a_srccomp + a_numcomp <= a_src.nComp());
    BL_ASSERT(a_dstcomp + a_numcomp <= nComp());

    Box locRegion = a_src.getRegion() & getRegion();
    bool sameRegBox = (a_src.m_regFAB.box() == m_regFAB.box());
         
    if (!locRegion.isEmpty())
    {
      Box region = locRegion & m_region;
      region &= a_src.m_region;
      for (BoxIterator boxit(region); boxit.ok(); ++boxit)
      {
        for (int icomp = 0; icomp < a_numcomp; icomp++)
        {
          int srcvar = a_srccomp + icomp;
          int dstvar = a_dstcomp + icomp;
          m_regFAB(boxit(), dstvar) *= a_src.m_regFAB(boxit(), srcvar);
        }
      }

      m_irrFAB.forall(a_src.m_irrFAB, locRegion, a_srccomp, a_dstcomp, a_numcomp, sameRegBox, [](Real& dest, const Real& src){dest*=src;});
    }
         
    return *this;
  }

  

  /**********************/
  EBCellFAB&
  EBCellFAB::operator/=(const EBCellFAB& a_src)
  {
    BL_ASSERT(a_src.nComp() == nComp());
         
    divide(a_src, 0, 0, nComp());
         
    return *this;
  }
  /**********************/
  EBCellFAB&
  EBCellFAB::divide(const EBCellFAB& a_src,
                    int a_srccomp,
                    int a_dstcomp,
                    int a_numcomp)
  {
         
    BL_ASSERT(isDefined());
    BL_ASSERT(a_src.isDefined());
    // Dan G. feels strongly that the assert below should NOT be commented out
    // Brian feels that a weaker version of the BL_ASSERT (if possible) is needed
    // Terry is just trying to get his code to work
    //BL_ASSERT(m_ebisBox == a_src.m_ebisBox);
         
    BL_ASSERT(a_srccomp + a_numcomp <= a_src.nComp());
    BL_ASSERT(a_dstcomp + a_numcomp <= nComp());
    bool sameRegBox = (a_src.m_regFAB.box() == m_regFAB.box());
         
    Box locRegion = a_src.getRegion() & getRegion();
    if (!locRegion.isEmpty())
    {
      Box region = locRegion & m_region;
      region &= a_src.m_region;
      for (BoxIterator boxit(region); boxit.ok(); ++boxit)
      {
        for (int icomp = 0; icomp < a_numcomp; icomp++)
        {
          int srcvar = a_srccomp + icomp;
          int dstvar = a_dstcomp + icomp;
          m_regFAB(boxit(), dstvar) /= a_src.m_regFAB(boxit(), srcvar);
        }
      }

      m_irrFAB.forall(a_src.m_irrFAB, locRegion, a_srccomp, a_dstcomp, a_numcomp, sameRegBox, [](Real& dest, const Real& src){dest/=src;});
    }
    return *this;
  }

 
  /**********************/
  EBCellFAB&
  EBCellFAB::operator+=(const Real& a_src)
  {
    BL_ASSERT(isDefined());
    for (BoxIterator boxit(m_region); boxit.ok(); ++boxit)
    {
      for (int icomp = 0; icomp < nComp(); icomp++)
      {
        m_regFAB(boxit(), icomp) += a_src;
      }
    }
         
    Real* l = m_irrFAB.dataPtr(0);
    int nvof = m_irrFAB.numVoFs();
    for (int i=0; i<nComp()*nvof; i++)
      l[i] += a_src;
         
    return *this;
  }
  /**********************/
  EBCellFAB&
  EBCellFAB::operator-=(const Real& a_src)
  {
    BL_ASSERT(isDefined());
    for (BoxIterator boxit(m_region); boxit.ok(); ++boxit)
    {
      for (int icomp = 0; icomp < nComp(); icomp++)
      {
        m_regFAB(boxit(), icomp) -= a_src;
      }
    }
         
    Real* l = m_irrFAB.dataPtr(0);
    int nvof = m_irrFAB.numVoFs();
    for (int i=0; i<nComp()*nvof; i++)
      l[i] -= a_src;
         
    return *this;
  }
  /**********************/
  EBCellFAB&
  EBCellFAB::operator*=(const Real& a_src)
  {
    BL_ASSERT(isDefined());
    for (BoxIterator boxit(m_region); boxit.ok(); ++boxit)
    {
      for (int icomp = 0; icomp < nComp(); icomp++)
      {
        m_regFAB(boxit(), icomp) *= a_src;
      }
    }
         
    Real* l = m_irrFAB.dataPtr(0);
    int nvof = m_irrFAB.numVoFs();
    for (int i=0; i<nComp()*nvof; i++)
      l[i] *= a_src;
         
    return *this;
  }
  /**********************/
  EBCellFAB&
  EBCellFAB::mult(Real a_src)
  {
    *this *= a_src;
         
    return *this;
  }
  /**********************/
  EBCellFAB&
  EBCellFAB::operator/=(const Real& a_src)
  {
    BL_ASSERT(isDefined());
    for (BoxIterator boxit(m_region); boxit.ok(); ++boxit)
    {
      for (int icomp = 0; icomp < nComp(); icomp++)
      {
        m_regFAB(boxit(), icomp) /= a_src;
      }
    }
         
    Real* l = m_irrFAB.dataPtr(0);
    int nvof = m_irrFAB.numVoFs();
    for (int i=0; i<nComp()*nvof; i++)
      l[i] /= a_src;
         
    return *this;
  }
         
  /**********************/
         
  //-----------------------------------------------------------------------
  Real
  EBCellFAB::max(int a_comp) const
  {
    BL_ASSERT(isDefined());
    Real val = -1.0e30;
         
    // Find the max on irregular cells.
    const EBISBox& ebbox = getEBISBox();
    const Box& box = BaseEBCellFAB<Real>::box();
    const IntVectSet validCells(box);
    for (VoFIterator vit(validCells, ebbox.getEBGraph()); vit.ok(); ++vit)
    {
      VolIndex vofi = vit();
      val = std::max(val, (*this)(vofi, a_comp));
    }
    return val;
  }

  Real
  EBCellFAB::max(const Box& subbox, int a_comp) const
  {
    BL_ASSERT(isDefined());
    Real val = -1.0e30;
         
    // Find the max on irregular cells.
    const EBISBox& ebbox = getEBISBox();
    const IntVectSet validCells(subbox);
    for (VoFIterator vit(validCells, ebbox.getEBGraph()); vit.ok(); ++vit)
    {
      VolIndex vofi = vit();
      val = std::max(val, (*this)(vofi, a_comp));
    }
    return val;
  }
  //-----------------------------------------------------------------------
         
  //-----------------------------------------------------------------------
  Real
  EBCellFAB::min(int a_comp) const
  {
    BL_ASSERT(isDefined());
    Real val = 1.0e30;
         
    // Find the min on irregular cells.
    const EBISBox& ebbox = getEBISBox();
    const Box& box = BaseEBCellFAB<Real>::box();
    const IntVectSet validCells(box);
    for (VoFIterator vit(validCells, ebbox.getEBGraph()); vit.ok(); ++vit)
    {
      VolIndex vofi = vit();
      val = std::min(val, (*this)(vofi, a_comp));
    }
    return val;
  }

  Real
  EBCellFAB::min(const Box& subbox, int a_comp) const
  {
    BL_ASSERT(isDefined());
    Real val = 1.0e30;
         
    // Find the min on irregular cells.
    const EBISBox& ebbox = getEBISBox();
    const IntVectSet validCells(subbox);
    for (VoFIterator vit(validCells, ebbox.getEBGraph()); vit.ok(); ++vit)
    {
      VolIndex vofi = vit();
      val = std::min(val, (*this)(vofi, a_comp));
    }
    return val;
  }
}         

const FABio&
EBCellFAB::getFABio ()
{
    return getFArrayBox().getFABio();
}
