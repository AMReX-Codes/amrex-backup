
#ifndef BL_MULTIPHAB_H

#define BL_MULTIPHAB_H


#include <cstdint>

#include <AMReX_BLassert.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FabArray.H>

namespace amrex
{

template <class T>
class PhData
{
    public:
    //! Construct an empty PhData, which must be resized  before use.
    PhData ();

    //!  Make BaseFab with desired domain (box) and number of components.
    // explicit BaseFab (const Box& bx,
    //                   int        n = 1,
	// 	      bool       alloc = true,
	// 	      bool       shared = false);
    //! The destructor deletes the array memory.
    // ~BaseFab ();
    void resize(N);
    void clear ();
    void define();
    
    //! Returns the number of components
    int nComp () const { return nvar; }

    //! Returns the number of points 
    long nPts () const { return numpts; }

    /**
    * \brief Returns a pointer to an object of type T that is the
    * value of the Nth component associated with the cell at the
    * low end of the domain.  This is commonly used to get a pointer
    * to data in the array which is then handed off to a Fortran
    * subroutine.  Remember that data is stored in Fortran array
    * order, with the component index coming last.   In other words,
    * dataPtr returns a pointer to all the Nth components.
    */
    T* dataPtr (int n = 0) { BL_ASSERT(!(dptr == 0)); return &dptr[n*numpts]; }
    //
    //! Same as above except works on const FABs.
    const T* dataPtr (int n = 0) const { BL_ASSERT(!(dptr == 0)); return &dptr[n*numpts]; }

    //! Returns true if the data for the FAB has been allocated.
    bool isAllocated () const { return dptr != 0; }

    /**
    * \brief This function puts numcomp component values, starting at
    * component N, from position pos in the domain into array data,
    * that must be allocated by the user.
    */
    void getVal (T*             data,
                 const IntVect& pos,
                 int            N,
                 int            numcomp) const;

    protected:

    private:

}


/**
 * \brief 
 * A collection (stored as an array) of EBArrayBox objects. 
 * An EBArrayBox object is like a FArrayBox, except that it addition
 * to a BaseFAB<Real> it also includes a BaseFAB<int> object that 
 * is a flag for the EB/multivalued content of the cell.
 * In the future it may also include additional BaseFAB<int> objects that
 * describe more complex connectivity between the cells
 *
 * Although this is called PhabArray, at the moment it is not templated
 * the way FabArray is, but rather extends FabArray<EBArrayBox>, where
 * EBArrayBox is derived form FArrayBox, which extends BaseFAB<Real>. 
 * In that sense it is more like MultiFab, but doesn't yet have the
 * capabilities that MultiFab adds to FabArray. In the future perhaps 
 * this can become a templated class that is the parent for a
 * MultiPhab drived from PhabArray<FEBArrayBox>, but not yet.
 *
 *
 * This class does NOT provide a copy constructor or assignment operator.
 */

class PhabArray
    :
    public FabArray<FArrayBox>
{
public:
    //
    // Constructs an empty MultiFab.  Data can be defined at a later
    // time using the define member functions inherited
    // from FabArray.
    //
    PhabArray ():FabArray(){
        std::cout << "Using multiphab" << std::endl;
    
    }


    PhabArray (const BoxArray&            bxs,
              const DistributionMapping& dm,
              int                        nvar,
              int                        ngrow,
	      const MFInfo&              info = MFInfo())
        : FabArray(bxs, dm, nvar, ngrow, info) {
        std::cout << "Using multiphab" << std::endl;
        
        }


protected:
    //
    // The EB data.
    //
    std::vector<BaseFab<int>*> m_eb_flag_fabs_v;

    // Map i,j,k,m to array of Real data
    // This will store the 2nd, 3rd, etc value of cells that
    // do not fit in the BaseFab<Real> pointed to by dptr
    //
    // Need to make this a class that holds the data in probably
    // a hash table so that e.g. extradata.DataPtr can be passed
    // to fortran
    std::vector<PhData> m_extra_data;
        
        std::map<std::tuple<int, int, int, int>, Real>> m_extra_data;

private:
};

// template <class T>
// void
// PhData<T>::define ()
// {
//     BL_ASSERT(nvar > 0);
//     BL_ASSERT(dptr == 0);
//     BL_ASSERT(numpts > 0);
//     BL_ASSERT(std::numeric_limits<long>::max()/nvar > numpts);
// 
//     truesize  = nvar*numpts;
//     dptr      = static_cast<T*>(amrex::The_Arena()->alloc(truesize*sizeof(T)));
//     ptr_owner = true;
//     //
//     // Now call T::T() on the raw memory so we have valid Ts.
//     //
//     T* ptr = dptr;
//     //
//     // Note this must be long not int for very large (e.g.,1024^3) boxes.
//     //
//     for (long i = 0; i < truesize; i++, ptr++)
//     {
//         new (ptr) T;
//     }
// 
//     amrex::update_fab_stats(numpts, truesize, sizeof(T));
// }
// 
// };

#endif /*BL_PHABARRAY*/
