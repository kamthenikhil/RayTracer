#ifdef _WIN32
#include <windows.h>
#endif

#ifndef __APPLE__
#include <GL/gl.h>
#include <GL/glut.h>
#else
#include <OpenGL/gl.h>
#include <GLUT/glut.h>
#endif

#include "../Headers/ray_tracer.h"
#include <iostream>
#include <map>
#include <fstream>
#include <string>

//--------------------------------------------------------------------------------
// Change these to test ray tracer at different resolutions
#define WIDTH  640
#define HEIGHT 480

clock_t startTime = clock();
clock_t endTime = 0;
bool printed = false;
//--------------------------------------------------------------------------------
Render_World world;
Driver driver(world);
//--------------------------------------------------------------------------------
void Display() {
	glClear(GL_COLOR_BUFFER_BIT);
	if (world.camera.film.colors)
		glDrawPixels(WIDTH, HEIGHT, GL_RGBA, GL_UNSIGNED_INT_8_8_8_8,
				(GLvoid*) world.camera.film.colors);
	glFlush();
}
//--------------------------------------------------------------------------------
void Handle_Idle() {
	driver.Render_More(startTime, endTime, printed);
	Display();
}
//--------------------------------------------------------------------------------
void Initialize_Opengl_And_Glut(int argc, char** argv) {
	glutInit(&argc, (char**) argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
	glutInitWindowSize(WIDTH, HEIGHT);
	glutCreateWindow("Ray Tracer");
	glutDisplayFunc(Display);
	glutIdleFunc(Handle_Idle);
	glClearColor(0, 0, 0, 0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(0.0, WIDTH, 0.0, HEIGHT, -1.0, 1.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}
//--------------------------------------------------------------------------------
void Usage(const std::string& exec) {
	std::cout << "Usage: " << exec << " <test number>"
			<< "\n\t <test number> 1-4" << std::endl;
	exit(1);
}
//--------------------------------------------------------------------------------
// Initialize_World
//    Sets up objects, materials, and lights for tests 1 - 4
//--------------------------------------------------------------------------------
void Initialize_World(Render_World& world, const int width, const int height,
		const int test_number, char* filename) {
	Vector_3D<double> color1(1, 0, 0);
	Vector_3D<double> color2(.2, .2, .8);
	Vector_3D<double> plane_color(1, 0.775, 0.5431);

	Vector_3D<double> color[4];
	color[0] = Vector_3D<double>(1, 0, 0);
	color[1] = Vector_3D<double>(0.2, 0.2, 0.8);
	color[2] = Vector_3D<double>(0.2, 1, 0.2);
	color[3] = Vector_3D<double>(0.7, 1, 0.1);

	ifstream input(filename);

	int counter = 2;
	switch (test_number) {
	case 1:
		while (!input.eof()) {
			int i, j, k;
			char comma;
			input >> i >> comma >> j >> comma >> k;
			Sphere* obj = new Sphere(Vector_3D<double>(i, j, k), 1.0);
			int c = counter % 4;
			obj->material_shader = new Flat_Shader(world, color[c]);
			world.objects.push_back(obj);
			counter += 1;
		}
		world.enable_shadows = false;
		world.recursion_depth_limit = 0;
		world.camera.Position_And_Aim_Camera(Vector_3D<double>(0, 4, -6),
				Vector_3D<double>(0, 1, 0), Vector_3D<double>(0, 1, 0));
		world.camera.Focus_Camera(1, (double) width / (double) height,
				(double) 70 / (double) 180 * PI);
		break;
	case 2:

		while (!input.eof()) {
			int i, j, k;
			char comma;
			input >> i >> comma >> j >> comma >> k;
			Sphere* obj = new Sphere(Vector_3D<double>(i, j, k), 1.0);
			int c = counter % 4;
			obj->material_shader = new Phong_Shader(world, color[c], color[c]);
			world.objects.push_back(obj);
			counter += 1;
		}
		world.enable_shadows = false;
		world.recursion_depth_limit = 0;
		world.camera.Position_And_Aim_Camera(Vector_3D<double>(0, 4, -6),
				Vector_3D<double>(0, 1, 0), Vector_3D<double>(0, 1, 0));
		world.camera.Focus_Camera(1, (double) width / (double) height,
				(double) 70 / (double) 180 * PI);
		break;
	case 3:
		while (!input.eof()) {
			int i, j, k;
			char comma;
			input >> i >> comma >> j >> comma >> k;
			Sphere* obj = new Sphere(Vector_3D<double>(i, j, k), 1.0);
			int c = counter % 4;
			obj->material_shader = new Reflective_Shader(world, color[c],
					color[c], Vector_3D<double>(1, 1, 1), 50, .5);
			world.objects.push_back(obj);
			counter += 1;
		}
		world.enable_shadows = true;
		world.recursion_depth_limit = 0;
		world.camera.Position_And_Aim_Camera(Vector_3D<double>(0, 4, -6),
				Vector_3D<double>(0, 1, 0), Vector_3D<double>(0, 1, 0));
		world.camera.Focus_Camera(1, (double) width / (double) height,
				(double) 70 / (double) 180 * PI);
		break;
	case 4:
		while (!input.eof()) {
			int i, j, k;
			char comma;
			input >> i >> comma >> j >> comma >> k;
			Sphere* obj = new Sphere(Vector_3D<double>(i, j, k), 1.0);
			int c = counter % 4;
			obj->material_shader = new Reflective_Shader(world, color[c],
					color[c], Vector_3D<double>(1, 1, 1), 50, .5);
			world.objects.push_back(obj);
			counter += 1;
		}

		world.enable_shadows = true;
		world.recursion_depth_limit = 4;
		world.camera.Position_And_Aim_Camera(Vector_3D<double>(0, 4, -6),
				Vector_3D<double>(-1, 1, 0), Vector_3D<double>(0, 1, 0));
		world.camera.Focus_Camera(1, (double) width / (double) height,
				(double) 70 / (double) 180 * PI);
		break;
	case 5:

		while (!input.eof()) {
			int i, j, k;
			char comma;
			input >> i >> comma >> j >> comma >> k;
			Sphere* obj = new Sphere(Vector_3D<double>(i, j, k), 1.0);
			int c = counter % 4;
			obj->material_shader = new Refractive_Shader(world, color[c],
					color[c], Vector_3D<double>(1, 1, 1), 50, .2, 1.005);
			world.objects.push_back(obj);
			counter += 1;
		}

		world.enable_shadows = true;
		world.recursion_depth_limit = 4;
		world.camera.Position_And_Aim_Camera(Vector_3D<double>(0, 4, -6),
				Vector_3D<double>(-1, 1, 0), Vector_3D<double>(0, 1, 0));
		world.camera.Focus_Camera(1, (double) width / (double) height,
				(double) 70 / (double) 180 * PI);
		break;
	default:
		std::cout << "Unrecognized test number" << std::endl;
		exit(1);
	}
	input.close();

	// lights
	Light* point_light = new Point_Light(Vector_3D<double>(-2, 7, -3),
			Vector_3D<double>(1, 1, 1), .25);
	Light* point_light2 = new Point_Light(Vector_3D<double>(3, 5, -3),
			Vector_3D<double>(1, 1, 1), .25);

	if (world.useIndex) {
		clock_t startTimeTreeConstruction = clock();
		world.kdtree = new KdTree(world.objects);
		cout << "Time taken for Index creation: "
				<< (double) (clock() - startTimeTreeConstruction)
						/ CLOCKS_PER_SEC << "s" << endl;
	}

	world.lights.push_back(point_light);
	world.lights.push_back(point_light2);

	// camera
	world.camera.film.Set_Resolution(width, height);
}
//--------------------------------------------------------------------------------
int main(int argc, char** argv) {
	int test_number = atoi(argv[2]);
	if (test_number < 1 || test_number > 5)
		Usage(argv[0]);

	Initialize_Opengl_And_Glut(argc, argv);
	std::string arg = argv[3];
	if (arg == "true") {
		world.useIndex = true;
	} else {
		world.useIndex = false;
	}
	Initialize_World(world, WIDTH, HEIGHT, test_number, argv[1]);

	glutMainLoop();
	return 0;
}
