
#ifndef BL_FABFACTORY_H_
#define BL_FABFACTORY_H_

#include <AMReX_Box.H>

namespace amrex
{

struct FabInfo
{
    bool alloc = true;
    bool shared = false;

    FabInfo& SetAlloc (bool a) {
        alloc = a;
        return *this;
    }

    FabInfo& SetShared (bool s) {
        shared = s;
        return *this;
    }
};

template <class FAB>
class FabFactory
{
public:
    virtual ~FabFactory () {}
    virtual FAB* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const = 0;
    virtual FAB* create (const FAB& rhs, MakeType make_type, int scomp, int ncomp) const;
    virtual FabFactory<FAB>* make_copy() const;
    FabFactory( ) { }
    FabFactory( const FabFactory& rhs) { }

};

template <class FAB>
class DefaultFabFactory
    : public FabFactory<FAB>
{
public:
    virtual FAB* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const override
    {
        return new FAB(box, ncomps, info.alloc, info.shared);
    }

    virtual FAB* create (const FAB& rhs, MakeType make_type, int scomp, int ncomp)
    {
        return new FAB(*rhs, make_type, scomp, ncomp);
    }

    virtual FabFactory<FAB>* make_copy() {return new DefaultFabFactory();}
};

}

#endif
