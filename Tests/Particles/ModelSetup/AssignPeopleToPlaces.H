#ifndef ASSIGN_PEOPLE_TO_PLACES_H_
#define ASSIGN_PEOPLE_TO_PLACES_H_

#include <AMReX.H>

#include <PersonContainer.H>
#include <PlaceContainer.H>

/**
   /brief Assigns people to places in a distributed fashion. 

   Each place type has a parameter to set the number of neighbors to consider
   for assignment. Setting '0' means no limit, and there is some probability
   that a person might get assigned to any place in the simulation, no matter
   how far away. This function implements the no limit option.

   The algorithm is as follows.
   
   First, each person jumps to any other random cell in the simulation. 
   We exclude people that have already been assigned. This assumes an 
   even probability distribution, but in principle it should be weighted by 
   the number of vacancies. In SpatialSim, this information is lagged and 
   recomputed every so many passes.
   
   Then, we redistribute to send the person to the right task.
   
   Next, each person selects a place in the current cell at random, also 
   using a uniform probability distribution. This should be modified to 
   pick only non-full places.
   
   If the place has room, and if the person passes a probability check, 
   we assign the person to the place. This check in principle should be 
   a function of distance between the person and the place.
   
   We repeat the above in a loop until all places are full. Finally, we
   reset all the people back to their original locations and redistribute.

   This approach seems fairly robust to modifying the the number of vacancies
   and number of people, except in the case that the two are equal, where
   convergence is very slow.   
*/
void assignPeopleToPlaces (PersonContainer& people, PlaceContainer& places);

#endif
