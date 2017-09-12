#define vtkRenderingCore_AUTOINIT 3(vtkInteractionStyle, vtkRenderingFreeType, vtkRenderingOpenGL2)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL2)

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkTriangle.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkColorTransferFunction.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkLine.h>
#include <vtkAppendPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkInteractorStyleTrackballCamera.h>

#include "Surfaces.h"
#include <iostream>
#include <vector>
#include <math.h>

// Set constant PI
const double PI = 3.14159265358979323846;

// Define the pointer to the different geometries
ParametricSurface* geometry;
Cylinder* cylinder = new Cylinder(0.0, 4.0, 0.0, 2.0 * PI);
Cone* cone = new Cone(0.0, 2.0, 0.0, 2.0 * PI);
Ellipsoid* ellipsoid = new Ellipsoid(-PI / 2.0, PI / 2.0, 0.0, 2.0 * PI);
Torus* torus = new Torus(0.0, 2.0 * PI, 0.0, 2.0 * PI);

// Define the pointers to the surface and line actors for each geometry
vtkSmartPointer<vtkActor> current_surface_actor = NULL;
vtkSmartPointer<vtkActor> current_lines_actor = NULL;
vtkSmartPointer<vtkActor> cylinder_surface_actor = NULL;
vtkSmartPointer<vtkActor> cone_surface_actor = NULL;
vtkSmartPointer<vtkActor> ellipsoid_surface_actor = NULL;
vtkSmartPointer<vtkActor> torus_surface_actor = NULL;
vtkSmartPointer<vtkActor> cylinder_lines_actor = NULL;
vtkSmartPointer<vtkActor> cone_lines_actor = NULL;
vtkSmartPointer<vtkActor> ellipsoid_lines_actor = NULL;
vtkSmartPointer<vtkActor> torus_lines_actor = NULL;

// The pointer to the renderer object
vtkSmartPointer<vtkRenderer> ren = NULL;

// Indicates if the curvature lines are being drawn
bool drawing_lines = true;


/*
	Loads the requested geometry
*/
void loadGeometry(int type)
{
	// +-------------------------+
	// | Initialize the geometry |
	// +-------------------------+

	// Indicate which geometry is to be loaded
	// 1 = cylinder
	// 2 = cone
	// 3 = ellipsoid
	// 4 = torus
	

	// The pointer to the parametric surface
	ParametricSurface* geometry = NULL;

	// Initialize the geometry based on the option value
	switch (type)
	{
		// Load the cylinder
	case 1:
		cylinder = new Cylinder(0.0, 4.0, 0.0, 2.0 * PI);
		geometry = cylinder;
		break;

		// Load the cone
	case 2:
		cone = new Cone(0.0, 2.0, 0.0, 2.0 * PI);
		geometry = cone;
		break;

		// Load the ellipsoid
	case 3:
		ellipsoid = new Ellipsoid(-PI / 2.0, PI / 2.0, 0.0, 2.0 * PI);
		geometry = ellipsoid;
		break;

		// Load the torus
	case 4:
		torus = new Torus(0.0, 2.0 * PI, 0.0, 2.0 * PI);
		geometry = torus;
		break;
	}


	// +---------------------------------------------+
	// | Generate the surface points and line points |
	// +---------------------------------------------+

	// Define the parameters increment for each rectangle
	int u_div = 50;
	int v_div = 50;
	double u_inc = (geometry->u_max - geometry->u_min) / double(u_div);
	double v_inc = (geometry->v_max - geometry->v_min) / double(v_div);

	// Calculate some initial values that are obtained from the geometric parameters
	// NOTE: Boundering points are calculated as well. It is possible the geometry is not continuous in the boundary
	int n_points = (u_div + 1) * (v_div + 1);

	// Initialize the array where point color values will be stored
	vtkSmartPointer<vtkDoubleArray> surface_colors = vtkSmartPointer<vtkDoubleArray>::New();
	surface_colors->SetNumberOfValues(n_points);

	// Set the point index (it is kept manually for taking advantage of the next for loops
	vtkIdType i = 0;

	// Initialize the object where surface points will be stored
	vtkSmartPointer<vtkPoints> surface_points = vtkSmartPointer<vtkPoints>::New();

	// Initialize the object where line points will be stored
	vtkSmartPointer<vtkPoints> line_points = vtkSmartPointer<vtkPoints>::New();

	// Initialize the v value
	double v_value = geometry->v_min;

	// Generate the points by traversing the parametric space
	// NOTE: Parametric space traversal is bottom - up
	for (int v = 0; v <= v_div; v += 1)
	{
		// Initialize the u value
		double u_value = geometry->u_min;

		for (int u = 0; u <= u_div; u += 1)
		{
			// Calculate the point and insert it
			Vector3 p = geometry->r(u_value, v_value);
			//std::cout << "(" << u_value << ", " << v_value << ") = [" << p.x << ", " << p.y << ", " << p.z << "]" << std::endl;
			surface_points->InsertNextPoint(p.x, p.y, p.z);

			// Calculate the partial derivatives at the point
			Vector3 pdu = geometry->pdr_pdu(u_value, v_value);
			Vector3 pdv = geometry->pdr_pdv(u_value, v_value);

			// Calculate the unit surface normal vector at the point
			Vector3 NORMAL = pdu.cross(pdv);
			NORMAL.normalize();
			//normal_vector.push_back(NORMAL);

			// Calculate the elements of the first fundamental form
			double E = pdu.dot(pdu);
			double F = pdu.dot(pdv);
			double G = pdv.dot(pdv);

			// Calculate the elements of the second fundamental form
			double L = geometry->pdr2_pdupdu(u_value, v_value).dot(NORMAL);
			double M = geometry->pdr2_pdupdv(u_value, v_value).dot(NORMAL);
			double N = geometry->pdr2_pdvpdv(u_value, v_value).dot(NORMAL);

			// Calculate the gaussian curvature at the point
			double K = ((L * N) - (M * M)) / ((E * G) - (F * F));

			// The directions of the principal curvatures are obtained from equation
			// (FN - GM)lambda^2 + (EN - GL)lambda + (EM - FL) = 0   (3.53 in hyperbook)

			// Calculate the attributes of the quadratic equation and the discriminant
			double a = (F * N) - (G * M);
			double b = (E * N) - (G * L);
			double c = (E * M) - (F * L);
			double discr = (b * b) - (4.0 * a * c);

			// If discriminant is greater or equal to zero then it means there exist at 
			// least one lambda. Therefore we can calculate the directions of the 
			// principal curvatures
			double lambda_min = NAN;
			double lambda_max = NAN;
			if (discr >= 0)
			{
				// Calculate the directions of the principal curvatures based on the values
				// of the quadratic equation
				if (a == 0)
				{
					if (b == 0)
					{
						lambda_min = 0.0;
						lambda_max = 0.0;
					}
					else
					{
						lambda_min = -c / b;
						lambda_max = -c / b;
					}
				}
				else
				{
					lambda_min = (-b - sqrt(discr)) / (2 * a);
					lambda_max = (-b + sqrt(discr)) / (2 * a);
				}

				// Put the minimum lambda value in lambda_min. Do the same for the maximum value and lambda_max
				if (lambda_min > lambda_max)
				{
					double tmp = lambda_min;
					lambda_min = lambda_max;
					lambda_max = tmp;
				}

				// Define the components of the parametric surface equation. Make du = 1
				// so lambda = dv
				Vector3 ru_min = pdu.clone();
				Vector3 rv_min = pdv.clone();
				rv_min.multiply(lambda_min);
				ru_min.sum(rv_min);
				ru_min.normalize();
				ru_min.multiply(0.1);

				Vector3 ru_max = pdu.clone();
				Vector3 rv_max = pdv.clone();
				ru_max.multiply(lambda_max);
				rv_max.sum(ru_max);
				rv_max.normalize();
				rv_max.multiply(0.1);

				// Insert the current points of the surface
				//line_points->InsertNextPoint(p.x, p.y, p.z);
				//line_points->InsertNextPoint(p.x + ru_min.x, p.y + ru_min.y, p.z + ru_min.z);
				line_points->InsertNextPoint(p.x - ru_min.x, p.y - ru_min.y, p.z - ru_min.z);
				line_points->InsertNextPoint(p.x + ru_min.x, p.y + ru_min.y, p.z + ru_min.z);

				//line_points->InsertNextPoint(p.x, p.y, p.z);
				//line_points->InsertNextPoint(p.x + rv_max.x, p.y + rv_max.y, p.z + rv_max.z);
				line_points->InsertNextPoint(p.x - rv_max.x, p.y - rv_max.y, p.z - rv_max.z);
				line_points->InsertNextPoint(p.x + rv_max.x, p.y + rv_max.y, p.z + rv_max.z);

			}

			// Set the value of the gaussian curvature as the color of the point
			surface_colors->SetValue(i, (K != NAN) ? K : 0.0);

			// Increment the point index
			i += 1;

			// Increment the u value
			u_value += u_inc;
		}

		// Increment the v value
		v_value += v_inc;
	}

	// +------------------------+
	// | Generate the triangles |
	// +------------------------+

	// Define the array where the triangles will be stored
	vtkSmartPointer<vtkCellArray> surface_triangles = vtkSmartPointer<vtkCellArray>::New();

	// Define the width of the grid (in terms of points)
	double width = u_div + 1;

	// Build the rectangles in the the parametric space
	for (int v = 0; v < v_div; v += 1)
	{
		for (int u = 0; u < u_div; u += 1)
		{
			// Initialize the points of the rectangle;
			int p0 = u + (v * width);
			int p1 = (u + 1) + (v * width);
			int p2 = (u + 1) + ((v + 1) * width);
			int p3 = u + ((v + 1) * width);

			// Define the bottom triangle
			vtkSmartPointer<vtkTriangle> triangle1 = vtkSmartPointer<vtkTriangle>::New();
			triangle1->GetPointIds()->SetId(0, p0);		// bottom left corner
			triangle1->GetPointIds()->SetId(1, p1);		// bottom right corner
			triangle1->GetPointIds()->SetId(2, p2);		// top right corner
			surface_triangles->InsertNextCell(triangle1);

			// Define the upper triangle
			vtkSmartPointer<vtkTriangle> triangle2 = vtkSmartPointer<vtkTriangle>::New();
			triangle2->GetPointIds()->SetId(0, p0);		// bottom left corner
			triangle2->GetPointIds()->SetId(1, p2);		// top right corner
			triangle2->GetPointIds()->SetId(2, p3);		// top left corner
			surface_triangles->InsertNextCell(triangle2);
		}
	}


	// +--------------------------------------+
	// | Add elements to the surface polydata |
	// +--------------------------------------+

	// Set surface points and color values
	vtkSmartPointer<vtkPolyData> surface_polydata = vtkSmartPointer<vtkPolyData>::New();
	surface_polydata->SetPoints(surface_points);
	surface_polydata->SetPolys(surface_triangles);
	surface_polydata->GetPointData()->SetScalars(surface_colors);


	// +-------------------------------------+
	// | Define the surface mapper and actor |
	// +-------------------------------------+

	// Create a color table
	vtkSmartPointer<vtkColorTransferFunction> surface_ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
	surface_ctf->AddRGBPoint(-1, 1.0, 0.0, 0.0);	// red for -1
	surface_ctf->AddRGBPoint(0, 1.0, 1.0, 1.0); 	// white for 0
	surface_ctf->AddRGBPoint(1, 0.0, 0.0, 1.0);		// blue for 1

													// Create a mapper for the surface
	vtkSmartPointer<vtkPolyDataMapper> surface_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	surface_mapper->SetInputData(surface_polydata);
	surface_mapper->SetLookupTable(surface_ctf);

	// Create an actor for the surface
	vtkSmartPointer<vtkActor> surface_actor = vtkSmartPointer<vtkActor>::New();
	surface_actor->SetMapper(surface_mapper);


	// +--------------------+
	// | Generate the lines |
	// +--------------------+

	// Get the number of generated line points
	int n_line_points = line_points->GetNumberOfPoints();

	// Initialize the array where line point color values will be stored
	vtkSmartPointer<vtkDoubleArray> line_colors = vtkSmartPointer<vtkDoubleArray>::New();
	line_colors->SetNumberOfValues(n_line_points);

	// Initialize the array where the curvature lines will be stored
	vtkSmartPointer<vtkCellArray> curvature_lines = vtkSmartPointer<vtkCellArray>::New();

	// Traverse the line points
	for (int i = 0; i < n_line_points; i += 4)
	{
		// Generate the lines
		vtkSmartPointer<vtkLine> line1 = vtkSmartPointer<vtkLine>::New();
		line1->GetPointIds()->SetId(0, i);
		line1->GetPointIds()->SetId(1, i + 1);

		vtkSmartPointer<vtkLine> line2 = vtkSmartPointer<vtkLine>::New();
		line2->GetPointIds()->SetId(0, i + 2);
		line2->GetPointIds()->SetId(1, i + 3);

		// Set the value of the gaussian curvature as the color of the point
		line_colors->SetValue(i, -1);
		line_colors->SetValue(i + 1, -1);
		line_colors->SetValue(i + 2, 1);
		line_colors->SetValue(i + 3, 1);

		// Insert the line in the lines array
		curvature_lines->InsertNextCell(line1);
		curvature_lines->InsertNextCell(line2);
	}


	// +------------------------------------------+
	// | Add elements to the line points polydata |
	// +------------------------------------------+

	// Set surface points and color values
	vtkSmartPointer<vtkPolyData> lines_polydata = vtkSmartPointer<vtkPolyData>::New();
	lines_polydata->SetPoints(line_points);
	lines_polydata->SetLines(curvature_lines);
	lines_polydata->GetPointData()->SetScalars(line_colors);


	// +-----------------------------------+
	// | Define the lines mapper and actor |
	// +-----------------------------------+

	// Create a color table for the lines
	vtkSmartPointer<vtkColorTransferFunction> lines_ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
	lines_ctf->AddRGBPoint(-1, 1.0, 0.0, 0.0);	// red for -1
	lines_ctf->AddRGBPoint(0, 1.0, 1.0, 1.0); 	// white for 0
	lines_ctf->AddRGBPoint(1, 0.0, 0.0, 1.0);	// blue for 1

												// Create a mapper for the lines
	vtkSmartPointer<vtkPolyDataMapper> lines_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	lines_mapper->SetInputData(lines_polydata);
	lines_mapper->SetLookupTable(lines_ctf);

	// Create an actor for the lines
	vtkSmartPointer<vtkActor> lines_actor = vtkSmartPointer<vtkActor>::New();
	lines_actor->SetMapper(lines_mapper);

	switch (type) 
	{
		case 1:
			cylinder_surface_actor = surface_actor;
			cylinder_lines_actor = lines_actor;
			break;

		case 2:
			cone_surface_actor = surface_actor;
			cone_lines_actor = lines_actor;
			break;

		case 3:
			ellipsoid_surface_actor = surface_actor;
			ellipsoid_lines_actor = lines_actor;
			break;

		case 4:
			torus_surface_actor = surface_actor;
			torus_lines_actor = lines_actor;
			break;
	}
}



/*
	The class required for keyboard interaction
*/
class KeyPressInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:

	static KeyPressInteractorStyle* New();
	vtkTypeMacro(KeyPressInteractorStyle, vtkInteractorStyleTrackballCamera);

	/*
		Function called when a keyboard key is pressed
	*/
	virtual void OnKeyPress()
	{
		// Get the keypress
		vtkRenderWindowInteractor *rwi = this->Interactor;
		std::string key = rwi->GetKeySym();

		// If key "1" is pressed then load the cylinder
		if (key == "1") 
		{
			// If the cylinder surface actors hasn't been generated then generate the geometry
			if (cylinder_surface_actor == NULL) 
			{
				// Load the geometry into
				loadGeometry(1);
			}
			
			// Remove the current actors
			ren->RemoveActor(current_surface_actor);
			ren->RemoveActor(current_lines_actor);

			// Set the cylinder actors as the current ones
			current_surface_actor = cylinder_surface_actor;
			current_lines_actor = cylinder_lines_actor;

			// Add the actors
			ren->AddActor(current_surface_actor);
			if (drawing_lines)
			{
				ren->AddActor(current_lines_actor);
			}

			// Reset the camera and display the new geometry
			ren->ResetCamera();
			ren->Modified();
		}
		
		// If key "2" is pressed then load the cone
		if (key == "2") 
		{
			// If the cone surface actors hasn't been generated then generate the geometry
			if (cone_surface_actor == NULL)
			{
				// Load the geometry into
				loadGeometry(2);
			}

			// Remove the current actors
			ren->RemoveActor(current_surface_actor);
			ren->RemoveActor(current_lines_actor);

			// Set the cone actors as the current ones
			current_surface_actor = cone_surface_actor;
			current_lines_actor = cone_lines_actor;

			// Add the actors
			ren->AddActor(current_surface_actor);
			if (drawing_lines)
			{
				ren->AddActor(current_lines_actor);
			}

			// Reset the camera and display the new geometry
			ren->ResetCamera();
			ren->Modified();
		}
		
		// If key "3" is pressed then load the ellipsoid
		if (key == "3")
		{
			// If the ellipsoid surface actors hasn't been generated then generate the geometry
			if (ellipsoid_surface_actor == NULL)
			{
				// Load the geometry into
				loadGeometry(3);
			}

			// Remove the current actors
			ren->RemoveActor(current_surface_actor);
			ren->RemoveActor(current_lines_actor);

			// Set the ellipsoid actors as the current ones
			current_surface_actor = ellipsoid_surface_actor;
			current_lines_actor = ellipsoid_lines_actor;

			// Add the actors
			ren->AddActor(current_surface_actor);
			if (drawing_lines)
			{
				ren->AddActor(current_lines_actor);
			}

			// Reset the camera and display the new geometry
			ren->ResetCamera();
			ren->Modified();
		}
		
		// If key "4" is pressed then load the torus
		if (key == "4")
		{
			// If the torus surface actors hasn't been generated then generate the geometry
			if (torus_surface_actor == NULL)
			{
				// Load the geometry into
				loadGeometry(4);
			}

			// Remove the current actors
			ren->RemoveActor(current_surface_actor);
			ren->RemoveActor(current_lines_actor);

			// Set the torus actors as the current ones
			current_surface_actor = torus_surface_actor;
			current_lines_actor = torus_lines_actor;

			// Add the actors
			ren->AddActor(current_surface_actor);
			if (drawing_lines)
			{
				ren->AddActor(current_lines_actor);
			}

			// Reset the camera and display the new geometry
			ren->ResetCamera();
			ren->Modified();
		}

		// If key "l" is pressed then toggle the principal curvatures
		if (key == "l") 
		{
			// Check whether the lines are being drawn
			if (drawing_lines) 
			{
				// Remove the lines actor
				ren->RemoveActor(current_lines_actor);
				drawing_lines = false;
			}
			else 
			{
				// Add the current lines actor
				ren->AddActor(current_lines_actor);
				drawing_lines = true;
			}

			// Reset the camera and display the new geometry
			ren->ResetCamera();
			ren->Modified();
		}

		

		// Forward events
		vtkInteractorStyleTrackballCamera::OnKeyPress();
	}

};
vtkStandardNewMacro(KeyPressInteractorStyle);


/*
	The main function
*/
int main() 
{
	// Initiate the program by generating the cylinder
	loadGeometry(1);

	// Set the cylinder actors as the current actors
	current_surface_actor = cylinder_surface_actor;
	current_lines_actor = cylinder_lines_actor;

	// Create a renderer, render window, and interactor
	ren = vtkSmartPointer<vtkRenderer>::New();

	vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
	renWin->AddRenderer(ren);
	renWin->SetSize(800, 800);

	// Define the object for keyboard events
	vtkSmartPointer<KeyPressInteractorStyle> keyboard = vtkSmartPointer<KeyPressInteractorStyle>::New();

	// Initialize the renderer window interactor
	vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	iren->SetRenderWindow(renWin);
	iren->SetInteractorStyle(keyboard);

	// Add the actors to the scene
	ren->AddActor(current_surface_actor);
	if (drawing_lines) 
	{
		ren->AddActor(current_lines_actor);
	}

	// Set the background color
	ren->SetBackground(.1, .1, .1);

	// Start the renderer window interactor
	iren->Initialize();

	// Start the event loop.
	iren->Start();

	return 0;
}
