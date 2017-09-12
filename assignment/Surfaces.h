#pragma once

#include "Vector3.h"

/*
	The abstract class for all the parametric surfaces
*/
class ParametricSurface
{

public:

	// Minimum and maximum values for the parametric values
	double u_min;
	double u_max;
	double v_min;
	double v_max;


	/*
		Constuctor of the class
	*/
	explicit
		ParametricSurface(double u_min, double u_max, double v_min, double v_max)
		: u_min(u_min), u_max(u_max), v_min(v_min), v_max(v_max) {}

	/*
		𝐫(𝑢,𝑣)
	*/
	virtual Vector3 r(double u, double v) = 0;

	/*
		∂𝐫/∂𝑢(𝑢,𝑣)
	*/
	virtual Vector3 pdr_pdu(double u, double v) = 0;

	/*
		∂𝐫/∂𝑣(𝑢,𝑣)
	*/
	virtual Vector3 pdr_pdv(double u, double v) = 0;

	/*
		∂2𝐫/∂𝑢∂𝑢(𝑢,𝑣)
	*/
	virtual Vector3 pdr2_pdupdu(double u, double v) = 0;

	/*
		∂2𝐫/∂𝑢∂𝑣(𝑢,𝑣)
	*/
	virtual Vector3 pdr2_pdupdv(double u, double v) = 0;

	/*
		∂2𝐫/∂𝑣∂𝑣(𝑢,𝑣)
	*/
	virtual Vector3 pdr2_pdvpdv(double u, double v) = 0;

};

/*
	The class representing the cylinder
*/
class Cylinder : public ParametricSurface
{
public:

	/*
		Constructor of the class
	*/
	Cylinder(double u_min, double u_max, double v_min, double v_max) : ParametricSurface(u_min, u_max, v_min, v_max) {}

	/*
		𝐫(𝑢,𝑣)
	*/
	Vector3 r(double u, double v);

	/*
		∂𝐫/∂𝑢(𝑢,𝑣)
	*/
	Vector3 pdr_pdu(double u, double v);

	/*
		∂𝐫/∂𝑣(𝑢,𝑣)
	*/
	Vector3 pdr_pdv(double u, double v);

	/*
		∂2𝐫/∂𝑢∂𝑢(𝑢,𝑣)
	*/
	Vector3 pdr2_pdupdu(double u, double v);

	/*
		∂2𝐫/∂𝑢∂𝑣(𝑢,𝑣)
	*/
	Vector3 pdr2_pdupdv(double u, double v);

	/*
		∂2𝐫/∂𝑣∂𝑣(𝑢,𝑣)
	*/
	Vector3 pdr2_pdvpdv(double u, double v);

};


/*
	The class representing the cone
*/
class Cone : public ParametricSurface
{
public:

	/*
		Constructor of the class
	*/
	Cone(double u_min, double u_max, double v_min, double v_max) : ParametricSurface(u_min, u_max, v_min, v_max) {}

	/*
		𝐫(𝑢,𝑣)
	*/
	Vector3 r(double u, double v);

	/*
		∂𝐫/∂𝑢(𝑢,𝑣)
	*/
	Vector3 pdr_pdu(double u, double v);

	/*
		∂𝐫/∂𝑣(𝑢,𝑣)
	*/
	Vector3 pdr_pdv(double u, double v);

	/*
		∂2𝐫/∂𝑢∂𝑢(𝑢,𝑣)
	*/
	Vector3 pdr2_pdupdu(double u, double v);

	/*
		∂2𝐫/∂𝑢∂𝑣(𝑢,𝑣)
	*/
	Vector3 pdr2_pdupdv(double u, double v);

	/*
		∂2𝐫/∂𝑣∂𝑣(𝑢,𝑣)
	*/
	Vector3 pdr2_pdvpdv(double u, double v);

};


/*
	The class representing the ellipsoid
*/
class Ellipsoid : public ParametricSurface
{
public:

	/*
		Constructor of the class
	*/
	Ellipsoid(double u_min, double u_max, double v_min, double v_max) : ParametricSurface(u_min, u_max, v_min, v_max) {}

	/*
		𝐫(𝑢,𝑣)
	*/
	Vector3 r(double u, double v);

	/*
		∂𝐫/∂𝑢(𝑢,𝑣)
	*/
	Vector3 pdr_pdu(double u, double v);

	/*
		∂𝐫/∂𝑣(𝑢,𝑣)
	*/
	Vector3 pdr_pdv(double u, double v);

	/*
		∂2𝐫/∂𝑢∂𝑢(𝑢,𝑣)
	*/
	Vector3 pdr2_pdupdu(double u, double v);

	/*
		∂2𝐫/∂𝑢∂𝑣(𝑢,𝑣)
	*/
	Vector3 pdr2_pdupdv(double u, double v);

	/*
		∂2𝐫/∂𝑣∂𝑣(𝑢,𝑣)
	*/
	Vector3 pdr2_pdvpdv(double u, double v);

};


/*
	The class representing the torus
*/
class Torus : public ParametricSurface
{
public:

	/*
		Constructor of the class
	*/
	Torus(double u_min, double u_max, double v_min, double v_max) : ParametricSurface(u_min, u_max, v_min, v_max) {}

	/*
		𝐫(𝑢,𝑣)
	*/
	Vector3 r(double u, double v);

	/*
		∂𝐫/∂𝑢(𝑢,𝑣)
	*/
	Vector3 pdr_pdu(double u, double v);

	/*
		∂𝐫/∂𝑣(𝑢,𝑣)
	*/
	Vector3 pdr_pdv(double u, double v);

	/*
		∂2𝐫/∂𝑢∂𝑢(𝑢,𝑣)
	*/
	Vector3 pdr2_pdupdu(double u, double v);

	/*
		∂2𝐫/∂𝑢∂𝑣(𝑢,𝑣)
	*/
	Vector3 pdr2_pdupdv(double u, double v);

	/*
		∂2𝐫/∂𝑣∂𝑣(𝑢,𝑣)
	*/
	Vector3 pdr2_pdvpdv(double u, double v);

};
