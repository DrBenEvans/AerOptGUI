/*********************************************
**
**	Created on: 	11/05/2015 2015
**	Author: 	matt - Matt Edmunds
**	File:			 Enumerations.h
**
**********************************************/

#ifndef ENUMERATIONS
#define ENUMERATIONS

namespace Enum
{

	//1 - Lift/Drag
	//2 - Distortion
	//3 - max Lift
	//4 - min Drag
	//5 - max Downforce
	//6 - min Lift
	enum ObjFunc
	{
		LIFTDRAG = 1,
		DISTORTION,
		MAXLIFT,
		MINDRAG,
		MAXDOWNFORCE,
        MINLIFT,
        FUNCNOTSET
	};

    enum OptMethod
    {
        MCS = 0,
        DE,
        PSO,
        METHODNOTSET
    };

	enum TreeType
	{
		DEFAULT = 1000,
		ROOT,
		DATA,
		MESH,
		FUNCTION,
		BOUNDARY,
		OPTIMISER,
		RUNTIME,
		NODE
	};

	enum Mesh
	{
		COURSE,
		MEDIUM,
		FINE
	};
}

#endif // ENUMERATIONS

