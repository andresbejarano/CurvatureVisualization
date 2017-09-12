#include "Surfaces.h"
#include <math.h>

Vector3 Cylinder::r(double u, double v)
{
	return Vector3(cos(v), sin(v), u);
}

Vector3 Cylinder::pdr_pdu(double u, double v)
{
	return Vector3(0, 0, 1);
}

Vector3 Cylinder::pdr_pdv(double u, double v)
{
	return Vector3(-sin(v), cos(v), 0);
}

Vector3 Cylinder::pdr2_pdupdu(double u, double v)
{
	return Vector3(0, 0, 0);
}

Vector3 Cylinder::pdr2_pdupdv(double u, double v)
{
	return Vector3(0, 0, 0);
}

Vector3 Cylinder::pdr2_pdvpdv(double u, double v)
{
	return Vector3(-cos(v), -sin(v), 0);
}

Vector3 Cone::r(double u, double v)
{
	return Vector3(0.5 * u * cos(v), 0.5 * u * sin(v), (sqrt(3.0) / 2.0) * u);
}

Vector3 Cone::pdr_pdu(double u, double v)
{
	return Vector3(0.5 * cos(v), 0.5 * sin(v), sqrt(3.0) / 2.0);
}

Vector3 Cone::pdr_pdv(double u, double v)
{
	return Vector3(-0.5 * u * sin(v), 0.5 * u * cos(v), 0);
}

Vector3 Cone::pdr2_pdupdu(double u, double v)
{
	return Vector3(0, 0, 0);
}

Vector3 Cone::pdr2_pdupdv(double u, double v)
{
	return Vector3(-0.5 * sin(v), 0.5 * cos(v), 0);
}

Vector3 Cone::pdr2_pdvpdv(double u, double v)
{
	return Vector3(-0.5 * u * cos(v), -0.5 * u * sin(v), 0);
}

Vector3 Ellipsoid::r(double u, double v)
{
	return Vector3(cos(u) * cos(v), cos(u) * sin(v), 2.0 * sin(u));
}

Vector3 Ellipsoid::pdr_pdu(double u, double v)
{
	return Vector3(-sin(u) * cos(v), -sin(u) * sin(v), 2.0 * cos(u));
}

Vector3 Ellipsoid::pdr_pdv(double u, double v)
{
	return Vector3(-cos(u) * sin(v), cos(u) * cos(v), 0);
}

Vector3 Ellipsoid::pdr2_pdupdu(double u, double v)
{
	return Vector3(-cos(u) * cos(v), -cos(u) * sin(v), -2.0 * sin(u));
}

Vector3 Ellipsoid::pdr2_pdupdv(double u, double v)
{
	return Vector3(sin(u) * sin(v), -sin(u) * cos(v), 0);
}

Vector3 Ellipsoid::pdr2_pdvpdv(double u, double v)
{
	return Vector3(-cos(u) * cos(v), -cos(u) * sin(v), 0);
}

Vector3 Torus::r(double u, double v)
{
	return Vector3((2.0 + cos(v)) * cos(u), (2.0 + cos(v)) * sin(u), sin(v));
}

Vector3 Torus::pdr_pdu(double u, double v)
{
	return Vector3(-(2.0 + cos(v)) * sin(u), (2.0 + cos(v)) * cos(u), 0);
}

Vector3 Torus::pdr_pdv(double u, double v)
{
	return Vector3(-sin(v) * cos(u), -sin(v) * sin(u), cos(v));
}

Vector3 Torus::pdr2_pdupdu(double u, double v)
{
	return Vector3(-(2.0 + cos(v)) * cos(u), -(2.0 + cos(v)) * sin(u), 0);
}

Vector3 Torus::pdr2_pdupdv(double u, double v)
{
	return Vector3(sin(v) * sin(u), -sin(v) * cos(u), 0);
}

Vector3 Torus::pdr2_pdvpdv(double u, double v)
{
	return Vector3(-cos(v) * cos(u), -cos(v) * sin(u), -sin(v));
}
