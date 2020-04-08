#include <AssignPeopleToPlaces.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

void assignPeopleToPlaces (PersonContainer& people, PlaceContainer& places)
{
    BL_PROFILE("assignPeopleToPlaces");

    const int lev = 0;

    ParmParse pp;
    int max_passes = 100;
    pp.query("max_passes", max_passes);
    
    for (int i = 0; i < max_passes; ++i)
    {
        amrex::Print() << "pass " << i << "\n";
        people.jumpToRandomCell();
        people.Redistribute();

#ifdef _OPENMP
#pragma omp parallel
#endif        
        for (PersonIterator pti(people, lev); pti.isValid(); ++pti)
        {
            auto& local_people = people.ParticlesAt(lev, pti);
            auto& local_places = places.ParticlesAt(lev, pti);
            
            const auto npeople = local_people.numParticles();
            const auto nplaces = local_places.numParticles();
            
            if (nplaces == 0) continue;
            
            auto people_data = local_people.GetArrayOfStructs()().data();
            auto places_data = local_places.GetArrayOfStructs()().data();
            for (int ip = 0; ip < npeople; ++ip)
            {
                auto& person = people_data[ip];
                if (person.idata(PersonIntID::place) >= 0) continue;  // already assigned

                const int m = amrex::Random_int(nplaces);
                auto& place = places_data[m];
                if (amrex::Random() < 0.25) // this should depend on distance.
                {
                    if (place.idata(PlaceIntID::capacity) > 0)
                    {
                        person.idata(PersonIntID::place) = place.idata(PlaceIntID::global_id);
                        place.idata(PlaceIntID::capacity) -= 1;
                    }
                }
            }
        }

        long num_vacancies = places.numVacancies();
        amrex::Print() << "number of vacancies: " << num_vacancies << "\n";

        long num_unassigned = people.numUnassigned();
        amrex::Print() << "number of unassigned: " << num_unassigned << "\n";
        if ( (num_unassigned == 0) or (num_vacancies == 0) ) break;
    }

    people.restoreOriginal();
    people.Redistribute();
}
