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
         
  EBCellFAB&
  EBCellFAB::
  axby(const EBCellFAB& a_X, 
       const EBCellFAB& a_Y,
       const Real& a_A, const Real& a_B)
  {
    for (BoxIterator boxit(m_region); boxit.ok(); ++boxit)
      {
        for (int icomp = 0; icomp < nComp(); icomp++)
        {
          m_regFAB(boxit(), icomp) = a_A*a_X.m_regFAB(boxit(), icomp) + a_B*a_Y.m_regFAB(boxit(), icomp); // 
        }
      }

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
