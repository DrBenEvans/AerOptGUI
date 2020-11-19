/*********************************************
**
**	Created on: 	11/05/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:			 Enumerations.h
**
**********************************************/

#ifndef ENUMERATIONS
#define ENUMERATIONS

/**
 * Enums used in the AerOpt Project including ObjectiveFunction, Optimisation Method, and Mesh Densities.
 */
namespace Enum
{

    /**
     * @brief The ObjFunc enum represents the different criteria that can be optimised.
     */
	enum ObjFunc
	{
        FUNCNOTSET,         // No setting selected
        LIFTDRAG,           // Lift-Drag Ratio
        DISTORTION,         //
        MAXLIFT,            // Maximise Lift
        MINDRAG,            // Minimise Drag
        MAXDOWNFORCE,       // Maximise Downforce
        MINLIFT,            // Minimise Lift

	};

    /**
     * @brief The OptMethod enum represents the different Optimisation Algorithms that can be used.
     */
    enum OptMethod
    {
        MCS,            // Modified Cuckoo search
        DE,             //Differential Evolution
        PSO,            // Particle Swarm Optimisation
        METHODNOTSET
    };

    /**
     * @brief The Mesh enum represents different Mesh density settings.
     */
	enum Mesh
	{
        COARSE,
		MEDIUM,
		FINE
	};

}

#endif // ENUMERATIONS

