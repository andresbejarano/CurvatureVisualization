#pragma once

/*
	The class representing a 3D vector
*/
class Vector3
{

public:

	// The element of the vector
	double x;
	double y;
	double z;

	/*
		Constructor of the class
	*/
	Vector3(double x, double y, double z) : x(x), y(y), z(z) {}

	/*
		Constructor of the class
	*/
	Vector3(double v[3]) : x(v[0]), y(v[1]), z(v[2]) {}

	/*
		Returns the magnitude of the vector
	*/
	double magnitude();

	/*
		Normalizes the vector
	*/
	void normalize();

	/*
		Multuply the vector by the given scalar
	*/
	void multiply(double s);

	/*
		Clones the vector
	*/
	Vector3 clone();

	/*
		Sums the value of the given vector
	*/
	void sum(Vector3 B);

	/*
		Subtracts the values of the given vector
	*/
	void sub(Vector3 B);

	/*
		Calculates the dot product with vector B
	*/
	double dot(Vector3 B);

	/*
		Calculates the cross product with vector B
	*/
	Vector3 cross(Vector3 B);

	/*
		Indicates whether the current vector has the same values as vector B
	*/
	bool equal(Vector3 B);

};
