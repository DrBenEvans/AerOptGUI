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
 * All enums used in the AerOpt Project
 */
namespace Enum
{

	//1 - Lift/Drag
	//2 - Distortion
	//3 - max Lift
	//4 - min Drag
	//5 - max Downforce
	//6 - min Lift

    /**
     * @brief The ObjFunc enum
     * Optimisation criteria
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
     * @brief The OptMethod enum
     * Optimisation Algorithms
     */
    enum OptMethod
    {
        METHODNOTSET = 0,
        MCS,            // Modified Cuckoo search
        DE,             // Differential Evolution
        PSO             // Particle Swarm Optimisation
    };

    /**
     * @brief The Mesh enum
     * Mesh density settings
     */
	enum Mesh
	{
        COARSE,
		MEDIUM,
		FINE
	};

}

#endif // ENUMERATIONS

