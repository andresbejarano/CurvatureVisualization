#include "Vector3.h"
#include <math.h>

double Vector3::magnitude()
{
	// Calculate the magnitude of the vector and return it
	return sqrt((this->x * this->x) + (this->y * this->y) + (this->z * this->z));
}

void Vector3::normalize()
{
	// Get the magnitude of the vector
	double mag = this->magnitude();

	// If it has a magnitude then normalize it
	if (mag > 0) 
	{
		this->x = this->x / mag;
		this->y = this->y / mag;
		this->z = this->z / mag;
	}
}

void Vector3::multiply(double s)
{
	// Multiply each component by the scalar
	this->x = this->x * s;
	this->y = this->y * s;
	this->z = this->z * s;
}

Vector3 Vector3::clone()
{
	return Vector3(this->x, this->y, this->z);
}

void Vector3::sum(Vector3 B)
{
	this->x += B.x;
	this->y += B.y;
	this->z += B.z;
}

void Vector3::sub(Vector3 B)
{
	this->x -= B.x;
	this->y -= B.y;
	this->z -= B.z;
}

double Vector3::dot(Vector3 B)
{
	return (this->x * B.x) + (this->y * B.y) + (this->z * B.z);
}

Vector3 Vector3::cross(Vector3 B)
{
	// Calculate the cross product values
	double i = (this->y * B.z) - (this->z * B.y);
	double j = (this->z * B.x) - (this->x * B.z);
	double k = (this->x * B.y) - (this->y * B.x);

	// Return the new vector
	return Vector3(i, j, k);
}

bool Vector3::equal(Vector3 B)
{
	return (this->x == B.x && this->y == B.y && this->z == B.z);
}
