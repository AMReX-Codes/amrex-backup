#ifndef BL_TYPETRAITS_H_
#define BL_TYPETRAITS_H_

#include <type_traits>

#include <AMReX_BaseFab.H>
#include <AMReX_BaseEBCellFAB.H>
namespace amrex
{

    // BaseFab 

    // Always exists
    template <class A, class Enable = void> struct IsBaseFab : std::false_type {};
    // Only exists if D is a BaseFab or child of (otherwise removed by SFINAE)
    template <class D>
    struct IsBaseFab<D, typename std::enable_if<std::is_base_of<BaseFab<typename D::value_type>,D>::value>::type> : std::true_type {};


    // BaseEBCellFAB

    // Always exists
    template <class A, class Enable = void> struct  IsBaseEBCellFAB : std::false_type {};
    // Only exists if D is a BaseEBCellFab or child of (otherwise removed by SFINAE)
    template <class D>
    struct IsBaseEBCellFAB<D, typename std::enable_if<std::is_base_of<BaseEBCellFAB<typename D::value_type>,D>::value>::type> : std::true_type {};


    // Always exists
    template <class A, class Enable = void> struct  IsContainer : std::false_type {};

    // Exists if IsBaseEBCEllFAB<D> exists
    template <class D>
    struct IsContainer<D, typename std::enable_if<IsBaseEBCellFAB<D>::value>::type > : std::true_type {};

    //// Exists if IsBaseFab<D> exists - but fails with 
    template <class D>   
    struct IsContainer<D, typename std::enable_if<IsBaseFab<D>::value>::type > : std::true_type {};

}

#endif
 