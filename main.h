// Shawn Halayka -- shalayka@gmail.com
// June 26, 2010
//
// This code and data is in the public domain.


#ifndef main_H
#define main_H
	

#include "uv_camera.h"

#include "custom_math.h"
using namespace custom_math;

#include "primitives.h"

#include "marching_squares.h"

#include <fstream>
using std::ofstream;

#include <ios>
using std::ios;

#include <vector>
#include <sstream>
using namespace std;

#include <random>
using std::mt19937;

#include <complex>
using namespace std;

#include <set>
using std::set;


size_t point_res = 10;


float grid_max = 1.5;
complex<float> C(0.2, 0.5);
unsigned short int max_iterations = 8;
float threshold = 4.0;
float beta = 2.0f;
bool mandelbrot_mode = FALSE;

vector_3 background_colour(1.0, 1.0, 1.0);

float outline_width = 3.0;
static const float outline_colour[] = {0.0, 0.0, 0.0};

bool draw_curves = true;
bool draw_mesh = true;
bool draw_outline = true;
bool draw_axis = true;
bool draw_control_list = true;
bool screenshot_mode = false;

uv_camera main_camera;

GLint win_id = 0;
GLint win_x = 800, win_y = 600;
float camera_w = 4;
float camera_fov = 45;
float camera_x_transform = 0;
float camera_y_transform = 0;
double u_spacer = 0.01;
double v_spacer = 0.5*u_spacer;
double w_spacer = 0.1;
double camera_near = 0.1;
double camera_far = 1000.0;

GLUquadricObj* glu_obj = gluNewQuadric(); // Probably should delete this before app exit... :)

bool lmb_down = false;
bool mmb_down = false;
bool rmb_down = false;
int mouse_x = 0;
int mouse_y = 0;

void idle_func(void);
void init_opengl(const int &width, const int &height);
void reshape_func(int width, int height);
void display_func(void);
void keyboard_func(unsigned char key, int x, int y);
void mouse_func(int button, int state, int x, int y);
void motion_func(int x, int y);
void passive_motion_func(int x, int y);
void draw_objects(bool disable_colouring = false);


	

vector<vector<vector_4>> all_4d_points;
vector<vector<vector_4>> pos;




vector<triangle> tris;
vector<vertex_3> face_normals;
vector<vertex_3> vertices;
vector<vertex_3> vertex_normals;

float mesh_transparent[] = { 0.0f, 0.5f, 1.0f, 0.2f };
float mesh_solid[] = { 0.0f, 0.5f, 1.0f, 1.0f };



complex<float> pow_complex(const complex<float>& in, const float beta)
{
	float fabs_beta = fabsf(beta);

	float self_dot = in.real() * in.real() + in.imag() * in.imag();

	if (self_dot == 0)
	{
		//qOut->x = 0;
//		qOut->y = 0;
	//	qOut->z = 0;
		//qOut->w = 0;

		return complex<float>(0, 0);
	}

	float len = std::sqrtf(self_dot);
	float self_dot_beta = std::powf(self_dot, fabs_beta / 2.0f);

	complex<float> out = complex<float>(
		self_dot_beta * std::cos(fabs_beta * std::acos(in.real() / len)),
		in.imag() * self_dot_beta * std::sin(fabs_beta * std::acos(in.real() / len)) / sqrtf(in.imag() * in.imag()));

	if (beta < 0)
		out = conj(out) / powf(abs(out), 2.0f);

	return out;
}

float iterate_mandelbrot_2d(vector< complex<float> >& trajectory_points,
	complex<float> Z,
	complex<float> C,
	const short unsigned int max_iterations,
	const float threshold,
	const float exponent)
{
	C = Z;
	Z = complex<float>(0, 0);

	trajectory_points.clear();
	trajectory_points.push_back(Z);

	for (short unsigned int i = 0; i < max_iterations; i++)
	{
		Z = pow_complex(Z, exponent);
		Z += C;

		trajectory_points.push_back(Z);

		if (abs(Z) >= threshold)
			break;
	}

	return abs(Z);
}





float iterate_julia_2d(vector< complex<float> >& trajectory_points,
	complex<float> Z,
	const complex<float> C,
	const short unsigned int max_iterations,
	const float threshold,
	const float exponent)
{
	trajectory_points.clear();
	trajectory_points.push_back(Z);

	for (short unsigned int i = 0; i < max_iterations; i++)
	{
		Z = pow_complex(Z, exponent);
		Z += C;

		trajectory_points.push_back(Z);

		if (abs(Z) >= threshold)
			break;
	}

	return abs(Z);
}


float iterate_2d(bool mandelbrot,
	vector< complex<float> >& trajectory_points,
	complex<float> Z,
	const complex<float> C,
	const short unsigned int max_iterations,
	const float threshold,
	const float exponent)
{
	if (mandelbrot)
	{
		return iterate_mandelbrot_2d(
			trajectory_points,
			Z,
			C,
			max_iterations,
			threshold,
			exponent);
	}
	else
	{
		return iterate_julia_2d(
			trajectory_points,
			Z,
			C,
			max_iterations,
			threshold,
			exponent);
	}

}

void get_vertices_and_normals_from_triangles(vector<triangle>& t, vector<vertex_3>& fn, vector<vertex_3>& v, vector<vertex_3>& vn)
{
	fn.clear();
	v.clear();
	vn.clear();

	if (0 == t.size())
		return;

	cout << "Triangles: " << t.size() << endl;

	cout << "Welding vertices" << endl;

	// Insert unique vertices into set.
	set<vertex_3> vertex_set;

	for (vector<triangle>::const_iterator i = t.begin(); i != t.end(); i++)
	{
		vertex_set.insert(i->vertex[0]);
		vertex_set.insert(i->vertex[1]);
		vertex_set.insert(i->vertex[2]);
	}

	cout << "Vertices: " << vertex_set.size() << endl;

	cout << "Generating vertex indices" << endl;

	// Add indices to the vertices.
	for (set<vertex_3>::const_iterator i = vertex_set.begin(); i != vertex_set.end(); i++)
	{
		size_t index = v.size();
		v.push_back(*i);
		v[index].index = index;
	}

	vertex_set.clear();

	// Re-insert modifies vertices into set.
	for (vector<vertex_3>::const_iterator i = v.begin(); i != v.end(); i++)
		vertex_set.insert(*i);

	cout << "Assigning vertex indices to triangles" << endl;

	// Find the three vertices for each triangle, by index.
	set<vertex_3>::iterator find_iter;

	for (vector<triangle>::iterator i = t.begin(); i != t.end(); i++)
	{
		find_iter = vertex_set.find(i->vertex[0]);
		i->vertex[0].index = find_iter->index;

		find_iter = vertex_set.find(i->vertex[1]);
		i->vertex[1].index = find_iter->index;

		find_iter = vertex_set.find(i->vertex[2]);
		i->vertex[2].index = find_iter->index;
	}

	vertex_set.clear();

	cout << "Calculating normals" << endl;
	fn.resize(t.size());
	vn.resize(v.size());

	for (size_t i = 0; i < t.size(); i++)
	{
		vertex_3 v0 = t[i].vertex[1] - t[i].vertex[0];
		vertex_3 v1 = t[i].vertex[2] - t[i].vertex[0];
		fn[i] = v0.cross(v1);
		fn[i].normalize();

		vn[t[i].vertex[0].index] = vn[t[i].vertex[0].index] + fn[i];
		vn[t[i].vertex[1].index] = vn[t[i].vertex[1].index] + fn[i];
		vn[t[i].vertex[2].index] = vn[t[i].vertex[2].index] + fn[i];
	}

	for (size_t i = 0; i < vn.size(); i++)
		vn[i].normalize();
}



void get_isosurface(
	const bool mandelbrot,
	const float grid_max,
	const size_t res,
	const complex<float> C,
	const unsigned short int max_iterations,
	const float threshold,
	const float exponent)
{
	tris.clear();
	face_normals.clear();
	vertices.clear();
	vertex_normals.clear();

	vector<float> image(res*res, 0);

	const float x_grid_min = -grid_max;
	const size_t x_res = res;
	const complex<float> x_step_size((grid_max - x_grid_min) / (x_res - 1), 0);

	const float y_grid_min = -grid_max;
	const size_t y_res = res;
	const complex<float> y_step_size(0, (grid_max - y_grid_min) / (y_res - 1));

	complex<float> Z(x_grid_min, y_grid_min);

	vector< complex<float> > trajectory_points;

	for (size_t x = 0; x < x_res; x++, Z += x_step_size)
	{
		Z = complex<float>(Z.real(), y_grid_min);

		for (size_t y = 0; y < y_res; y++, Z += y_step_size)
		{
			image[x_res * y + x] = iterate_2d(mandelbrot, trajectory_points, Z, C, max_iterations, threshold, exponent);

			if (image[x_res * y + x] > threshold*2.0f)
				image[x_res * y + x] = threshold*2.0f;

			image[x_res * y + x] /= 2.0f * threshold;
			image[x_res * y + x] = 1 - image[x_res * y + x];
		}
	}


	grid_square g;

	cout << "Generating geometric primitives..." << endl;
	cout << endl;

	double grid_x_pos = x_grid_min; // Start at minimum x.
	double grid_y_pos = y_grid_min; // Start at maximum y.

	float step_size = (grid_max - x_grid_min) / static_cast<double>(res - 1);


	// Begin march.
	for (short unsigned int y = 0; y < (y_res - 1); y++, grid_y_pos += step_size, grid_x_pos = x_grid_min)
	{
		for (short unsigned int x = 0; x < (x_res - 1); x++, grid_x_pos += step_size)
		{
			// Corner vertex order: 03
			//                      12
			// e.g.: clockwise, as in OpenGL
			g.vertex[0] = vertex_3(grid_x_pos, grid_y_pos, 0, 0);
			g.vertex[1] = vertex_3(grid_x_pos, grid_y_pos - step_size, 0, 0);
			g.vertex[2] = vertex_3(grid_x_pos + step_size, grid_y_pos - step_size, 0, 0);
			g.vertex[3] = vertex_3(grid_x_pos + step_size, grid_y_pos, 0, 0);

			g.value[0] = image[y * x_res + x];
			g.value[1] = image[(y + 1) * x_res + x];
			g.value[2] = image[(y + 1) * x_res + (x + 1)];
			g.value[3] = image[y * x_res + (x + 1)];

			g.generate_primitives(tris, 0.5F);
		}
	}

	get_vertices_and_normals_from_triangles(tris, face_normals, vertices, vertex_normals);

}


// https://stackoverflow.com/questions/785097/how-do-i-implement-a-bézier-curve-in-c
vector_4 getBezierPoint(vector<vector_4> points, float t)
{
	size_t i = points.size() - 1;

	while (i > 0)
	{
		for (size_t k = 0; k < i; k++)
		{
			points[k].x += t * (points[k + 1].x - points[k].x);
			points[k].y += t * (points[k + 1].y - points[k].y);
			points[k].z += t * (points[k + 1].z - points[k].z);
			points[k].w += t * (points[k + 1].w - points[k].w);
		}

		i--;
	}

	return points[0];
}



void get_points(
	const bool mandelbrot,
	const float grid_max,
	const size_t res,
	const complex<float> C,
	const unsigned short int max_iterations,
	const float threshold,
const float exponent)
{
	all_4d_points.clear();
	pos.clear();

	const float x_grid_max = grid_max;
	const float x_grid_min = -x_grid_max;
	const size_t x_res = res;
	const complex<float> x_step_size((x_grid_max - x_grid_min) / (x_res - 1), 0);

	const float y_grid_max = grid_max;
	const float y_grid_min = -y_grid_max;
	const size_t y_res = res;
	const complex<float> y_step_size(0, (y_grid_max - y_grid_min) / (y_res - 1));

	complex<float> Z(x_grid_min, y_grid_min);

	vector< complex<float> > trajectory_points;

	for (size_t x = 0; x < x_res; x++, Z += x_step_size)
	{
		Z = complex<float>(Z.real(), y_grid_min);

		for (size_t y = 0; y < y_res; y++, Z += y_step_size)
		{
			float magnitude = iterate_2d(mandelbrot, trajectory_points, Z, C, max_iterations, threshold, exponent);

			if (magnitude < threshold)
			{
				vector<vector_4> v;

				for (size_t i = 0; i < trajectory_points.size(); i++)
				{
					vector_4 p;
					p.x = trajectory_points[i].real();
					p.y = trajectory_points[i].imag();
					v.push_back(p);
				}

				all_4d_points.push_back(v);
			}
		}
	}


	cout << "trajectory count " << all_4d_points.size() << endl;

	size_t orbit_count = 0;

	for (size_t i = 0; i < all_4d_points.size(); i++)
	{
		set<vector_4> point_set;

		for (size_t j = 0; j < all_4d_points[i].size(); j++)
		{
			point_set.insert(all_4d_points[i][j]);
		}

		if (point_set.size() != all_4d_points[i].size())
		{
			cout << point_set.size() << endl;
			orbit_count++;
		}
	}

	cout << "orbit count " << orbit_count << endl;






	for (size_t i = 0; i < all_4d_points.size(); i++)
	{
		vector<vector_4> p;

		for (float t = 0; t <= 0.2f; t += 0.01f)
		//for (float t = 0; t <= 0.85f; t += 0.01f)
		{
			vector_4 v = getBezierPoint(all_4d_points[i], t);
			p.push_back(v);
		}

		pos.push_back(p);
	}


	 get_isosurface(
		mandelbrot,
		grid_max,
		1000,
		C,
		max_iterations,
		threshold,
		exponent);



}


// TODO: fix camera bug where portrait mode crashes.
void take_screenshot(size_t num_cams_wide, const char *filename, const bool reverse_rows = false)
{
	screenshot_mode = true;

	//get_points(res);

	// Set up Targa TGA image data.
	unsigned char  idlength = 0;
	unsigned char  colourmaptype = 0;
	unsigned char  datatypecode = 2;
	unsigned short int colourmaporigin = 0;
	unsigned short int colourmaplength = 0;
	unsigned char  colourmapdepth = 0;
	unsigned short int x_origin = 0;
	unsigned short int y_origin = 0;

	cout << "Image size: " << static_cast<size_t>(win_x)*num_cams_wide << "x" << static_cast<size_t>(win_y)*num_cams_wide << " pixels" << endl;

	if (static_cast<size_t>(win_x)*num_cams_wide > static_cast<unsigned short>(-1) ||
		static_cast<size_t>(win_y)*num_cams_wide > static_cast<unsigned short>(-1))
	{
		cout << "Image too large. Maximum width and height is " << static_cast<unsigned short>(-1) << endl;
		return;
	}

	unsigned short int px = win_x*static_cast<unsigned short>(num_cams_wide);
	unsigned short int py = win_y*static_cast<unsigned short>(num_cams_wide);
	unsigned char  bitsperpixel = 24;
	unsigned char  imagedescriptor = 0;
	vector<char> idstring;

	size_t num_bytes = 3*px*py;
	vector<unsigned char> pixel_data(num_bytes);

	// Adjust some parameters for large screen format.
	bool temp_draw_control_list = draw_control_list;
	draw_control_list = false;

	float temp_outline_width = outline_width;
	outline_width = 6;

	vector<unsigned char> fbpixels(3*win_x*win_y);

	const size_t total_cams = num_cams_wide * num_cams_wide;
	size_t cam_count = 0;
	// Loop through subcameras.
	for(size_t cam_num_x = 0; cam_num_x < num_cams_wide; cam_num_x++)
	{
		for(size_t cam_num_y = 0; cam_num_y < num_cams_wide; cam_num_y++)
		{
			cout << "Camera: " << cam_count + 1 << " of " << total_cams << endl;

			// Set up camera, draw, then copy the frame buffer.
			main_camera.Set_Large_Screenshot(num_cams_wide, cam_num_x, cam_num_y);
			display_func();
			glReadPixels(0, 0, win_x, win_y, GL_RGB, GL_UNSIGNED_BYTE, &fbpixels[0]);

			// Copy pixels to large image.
			for(GLint i = 0; i < win_x; i++)
			{
				for(GLint j = 0; j < win_y; j++)
				{
					size_t fb_index = 3*(j*win_x + i);

					size_t screenshot_x = cam_num_x*win_x + i;
					size_t screenshot_y = cam_num_y*win_y + j;
					size_t screenshot_index = 3*(screenshot_y*(win_x*num_cams_wide) + screenshot_x);

					pixel_data[screenshot_index] = fbpixels[fb_index + 2];
					pixel_data[screenshot_index + 1] = fbpixels[fb_index + 1];
					pixel_data[screenshot_index + 2] = fbpixels[fb_index ];
				}
			}

			cam_count++;
		}

	}

	screenshot_mode = false;

	// Restore the parameters.
	draw_control_list = temp_draw_control_list;
	outline_width = temp_outline_width;
	main_camera.Set();

	// Write Targa TGA file to disk.
	ofstream out(filename, ios::binary);

	if(!out.is_open())
	{
		cout << "Failed to open TGA file for writing: " << filename << endl;
		return;	
	}

	out.write(reinterpret_cast<char *>(&idlength), 1);
	out.write(reinterpret_cast<char *>(&colourmaptype), 1);
	out.write(reinterpret_cast<char *>(&datatypecode), 1);
	out.write(reinterpret_cast<char *>(&colourmaporigin), 2);
	out.write(reinterpret_cast<char *>(&colourmaplength), 2);
	out.write(reinterpret_cast<char *>(&colourmapdepth), 1);
	out.write(reinterpret_cast<char *>(&x_origin), 2);
	out.write(reinterpret_cast<char *>(&y_origin), 2);
	out.write(reinterpret_cast<char *>(&px), 2);
	out.write(reinterpret_cast<char *>(&py), 2);
	out.write(reinterpret_cast<char *>(&bitsperpixel), 1);
	out.write(reinterpret_cast<char *>(&imagedescriptor), 1);

	out.write(reinterpret_cast<char *>(&pixel_data[0]), num_bytes);

	//get_points(point_res);
}

void idle_func(void)
{
	glutPostRedisplay();
}




void init_opengl(const int &width, const int &height)
{
	win_x = width;
	win_y = height;

	if(win_x < 1)
		win_x = 1;

	if(win_y < 1)
		win_y = 1;

	glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_ALPHA|GLUT_DEPTH);
	glutInitWindowPosition(0, 0);
	glutInitWindowSize(win_x, win_y);
	win_id = glutCreateWindow("2D v4");

	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glDepthMask(GL_TRUE);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
	glShadeModel(GL_SMOOTH);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_LINE_SMOOTH);
	glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_POLYGON_SMOOTH);

	float light_colour[] = {0.9f, 0.9f, 0.9f, 1.0f};
	float light0_position[] = {0.0f, 0.0f, 1.0f, 0.0f };
	float light1_position[] = {0.0f, 0.0f, -1.0f, 0.0f };
	float light2_position[] = {1.0f, 0.0f, 0.0f, 0.0f };
	float light3_position[] = {-1.0f, 0.0f, 0.0f, 0.0f };

	glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_colour);
	glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_colour);
	glLightfv(GL_LIGHT2, GL_POSITION, light2_position);
	glLightfv(GL_LIGHT2, GL_DIFFUSE, light_colour);
	glLightfv(GL_LIGHT3, GL_POSITION, light3_position);
	glLightfv(GL_LIGHT3, GL_DIFFUSE, light_colour);

	float light_colour2[] = {0.5f, 0.5f, 0.5f, 1.0f};
	float light4_position[] = {0.0f, 1.0f, 0.0f, 0.0f };
	float light5_position[] = {0.0f, -1.0f, 0.0f, 0.0f };

	glLightfv(GL_LIGHT4, GL_POSITION, light4_position);
	glLightfv(GL_LIGHT4, GL_DIFFUSE, light_colour2);
	glLightfv(GL_LIGHT5, GL_POSITION, light5_position);
	glLightfv(GL_LIGHT5, GL_DIFFUSE, light_colour2);

	glClearColor(static_cast<float>(background_colour.x), static_cast<float>(background_colour.y), static_cast<float>(background_colour.z), 1.0f);
	glClearDepth(1.0f);

	main_camera.Set(0, 0, camera_w, camera_fov, win_x, win_y, camera_near, camera_far);

	get_points(mandelbrot_mode, grid_max, point_res, C, max_iterations, threshold, beta);

}

void reshape_func(int width, int height)
{
	win_x = width;
	win_y = height;

	if(win_x < 1)
		win_x = 1;

	if(win_y < 1)
		win_y = 1;

	glutSetWindow(win_id);
	glutReshapeWindow(win_x, win_y);
	glViewport(0, 0, win_x, win_y);

	main_camera.Set(main_camera.u, main_camera.v, main_camera.w, main_camera.fov, win_x, win_y, camera_near, camera_far);
}

// Text drawing code originally from "GLUT Tutorial -- Bitmap Fonts and Orthogonal Projections" by A R Fernandes
void render_string(int x, const int y, void *font, const string &text)
{
	for(size_t i = 0; i < text.length(); i++)
	{
		glRasterPos2i(x, y);
		glutBitmapCharacter(font, text[i]);
		x += glutBitmapWidth(font, text[i]) + 1;
	}
}
// End text drawing code.

void display_func(void)
{
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);



	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_LIGHT2);
	glEnable(GL_LIGHT3);
	glEnable(GL_LIGHT4);
	glEnable(GL_LIGHT5);

	glEnable(GL_ALPHA);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glMaterialfv(GL_FRONT, GL_DIFFUSE, mesh_transparent);

	glBegin(GL_TRIANGLES);

	for (size_t i = 0; i < tris.size(); i++)
	{
		size_t v_index0 = tris[i].vertex[0].index;
		size_t v_index1 = tris[i].vertex[1].index;
		size_t v_index2 = tris[i].vertex[2].index;

		glNormal3f(vertex_normals[v_index0].x, vertex_normals[v_index0].y, vertex_normals[v_index0].z);
		glVertex3f(vertices[v_index0].x, vertices[v_index0].y, vertices[v_index0].z);
		glNormal3f(vertex_normals[v_index1].x, vertex_normals[v_index1].y, vertex_normals[v_index1].z);
		glVertex3f(vertices[v_index1].x, vertices[v_index1].y, vertices[v_index1].z);
		glNormal3f(vertex_normals[v_index2].x, vertex_normals[v_index2].y, vertex_normals[v_index2].z);
		glVertex3f(vertices[v_index2].x, vertices[v_index2].y, vertices[v_index2].z);
	}

	glEnd();

	glDisable(GL_BLEND);

	if(true == draw_outline)
	{
		glDisable(GL_DEPTH_TEST);

		// Draw outline code from NeHe lesson 37:
		// http://nehe.gamedev.net/data/lessons/lesson.asp?lesson=37
		glPushAttrib(GL_ALL_ATTRIB_BITS);
		glLineWidth(outline_width);
		glCullFace(GL_BACK);
		glPolygonMode(GL_FRONT, GL_LINE);
		glColor3fv(&outline_colour[0]);

		draw_objects(true);


		glPopAttrib();
		// End draw outline code.

		glEnable(GL_DEPTH_TEST);
	}


	// Draw the model's components using OpenGL/GLUT primitives.



	draw_objects();




	if(true == draw_outline)
	{
		// Draw outline code from NeHe lesson 37:
		// http://nehe.gamedev.net/data/lessons/lesson.asp?lesson=37
		glPushAttrib(GL_ALL_ATTRIB_BITS);
		glLineWidth(outline_width);
		glCullFace(GL_FRONT);
		glPolygonMode(GL_BACK, GL_LINE);

		draw_objects(true);

		glPopAttrib();
		// End draw outline code.
	}

	if(true == draw_control_list)
	{
		// Text drawing code originally from "GLUT Tutorial -- Bitmap Fonts and Orthogonal Projections" by A R Fernandes
		// http://www.lighthouse3d.com/opengl/glut/index.php?bmpfontortho
		glMatrixMode(GL_PROJECTION);
		glPushMatrix();
		glLoadIdentity();
		gluOrtho2D(0, win_x, 0, win_y);
		glScalef(1, -1, 1); // Neat. :)
		glTranslatef(0, -static_cast<float>(win_y), 0); // Neat. :)
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
		glLoadIdentity();

		glColor3f(0, 0, 0);

		size_t break_size = 13;
		size_t start = 20;

		render_string(10, start, GLUT_BITMAP_HELVETICA_12, string("Keyboard controls:"));

		render_string(10, start + 1*break_size, GLUT_BITMAP_HELVETICA_10, string("H: Draw outlines"));
		render_string(10, start + 2*break_size, GLUT_BITMAP_HELVETICA_10, string("J: Draw axis"));
		render_string(10, start + 3*break_size, GLUT_BITMAP_HELVETICA_10, string("K: Draw this list"));
		render_string(10, start + 4*break_size, GLUT_BITMAP_HELVETICA_10, string("L: Take screenshot"));

		glPopMatrix();
		glMatrixMode(GL_PROJECTION);
		glPopMatrix();
		//glMatrixMode(GL_MODELVIEW);
		// End text drawing code.
	}


	glEnable(GL_DEPTH_TEST);







	if (false == screenshot_mode)
	{
		glFlush();
		glutSwapBuffers();
	}
}

void keyboard_func(unsigned char key, int x, int y)
{
	switch (tolower(key))
	{
	case 'y':
	{
		draw_mesh = !draw_mesh;
		break;
	}
	case 'u':
	{
		draw_curves = !draw_curves;
		break;
	}
	case 'h':
	{
		draw_outline = !draw_outline;
		break;
	}
	case 'j':
	{
		draw_axis = !draw_axis;
		break;
	}
	case 'k':
	{
		draw_control_list = !draw_control_list;
		break;
	}
	case 'n':
	{
		get_points(mandelbrot_mode, grid_max, point_res, C, max_iterations, threshold, beta);
		take_screenshot(8, "screenshot.tga");
		get_points(mandelbrot_mode, grid_max, point_res, C, max_iterations, threshold, beta);
		break;
	}
	case 'm':
	{
		get_points(mandelbrot_mode, grid_max, 50, C, max_iterations, threshold, beta);
		take_screenshot(8, "screenshot.tga");
		get_points(mandelbrot_mode, grid_max, point_res, C, max_iterations, threshold, beta);
		break;
	}
	case 'b':
	{
		vector<float> betas = { -8, -7, -6, -5, -4, -3, 2, 3, 4, 5, 6, 7, 8 };

		for (size_t i = 0; i < betas.size(); i++)
		{
			ostringstream oss;
			oss << betas[i];

			string filename = "julia2d_" + oss.str() + ".tga";

			cout << filename << endl;

			get_points(false, grid_max, 20, C, max_iterations, threshold, betas[i]);
			take_screenshot(8, filename.c_str());
		}
	
		for (size_t i = 0; i < betas.size(); i++)
		{
			ostringstream oss;
			oss << betas[i];

			string filename = "mandelbrot2d_" + oss.str() + ".tga";

			cout << filename << endl;

			get_points(true, grid_max, 20, C, max_iterations, threshold, betas[i]);
			take_screenshot(2, filename.c_str());
		}

		break;
	}
	default:
		break;
	}
}




void mouse_func(int button, int state, int x, int y)
{
	if(GLUT_LEFT_BUTTON == button)
	{
		if(GLUT_DOWN == state)
			lmb_down = true;
		else
			lmb_down = false;
	}
	else if(GLUT_MIDDLE_BUTTON == button)
	{
		if(GLUT_DOWN == state)
			mmb_down = true;
		else
			mmb_down = false;
	}
	else if(GLUT_RIGHT_BUTTON == button)
	{
		if(GLUT_DOWN == state)
			rmb_down = true;
		else
			rmb_down = false;
	}
}

void motion_func(int x, int y)
{
	int prev_mouse_x = mouse_x;
	int prev_mouse_y = mouse_y;

	mouse_x = x;
	mouse_y = y;

	int mouse_delta_x = mouse_x - prev_mouse_x;
	int mouse_delta_y = prev_mouse_y - mouse_y;

	if(true == lmb_down && (0 != mouse_delta_x || 0 != mouse_delta_y))
	{
		main_camera.u -= static_cast<float>(mouse_delta_y)*u_spacer;
		main_camera.v += static_cast<float>(mouse_delta_x)*v_spacer;
	}
	else if(true == rmb_down && (0 != mouse_delta_y))
	{
		main_camera.w -= static_cast<float>(mouse_delta_y)*w_spacer;
	}

	main_camera.Set(); // Calculate new camera vectors.
}

void passive_motion_func(int x, int y)
{
	mouse_x = x;
	mouse_y = y;
}



class RGB
{
public:
	unsigned char r, g, b;
};

RGB HSBtoRGB(unsigned short int hue_degree, unsigned char sat_percent, unsigned char bri_percent)
{
	float R = 0.0f;
	float G = 0.0f;
	float B = 0.0f;

	if (hue_degree > 359)
		hue_degree = 359;

	if (sat_percent > 100)
		sat_percent = 100;

	if (bri_percent > 100)
		bri_percent = 100;

	float hue_pos = 6.0f - ((static_cast<float>(hue_degree) / 359.0f) * 6.0f);

	if (hue_pos >= 0.0f && hue_pos < 1.0f)
	{
		R = 255.0f;
		G = 0.0f;
		B = 255.0f * hue_pos;
	}
	else if (hue_pos >= 1.0f && hue_pos < 2.0f)
	{
		hue_pos -= 1.0f;

		R = 255.0f - (255.0f * hue_pos);
		G = 0.0f;
		B = 255.0f;
	}
	else if (hue_pos >= 2.0f && hue_pos < 3.0f)
	{
		hue_pos -= 2.0f;

		R = 0.0f;
		G = 255.0f * hue_pos;
		B = 255.0f;
	}
	else if (hue_pos >= 3.0f && hue_pos < 4.0f)
	{
		hue_pos -= 3.0f;

		R = 0.0f;
		G = 255.0f;
		B = 255.0f - (255.0f * hue_pos);
	}
	else if (hue_pos >= 4.0f && hue_pos < 5.0f)
	{
		hue_pos -= 4.0f;

		R = 255.0f * hue_pos;
		G = 255.0f;
		B = 0.0f;
	}
	else
	{
		hue_pos -= 5.0f;

		R = 255.0f;
		G = 255.0f - (255.0f * hue_pos);
		B = 0.0f;
	}

	if (100 != sat_percent)
	{
		if (0 == sat_percent)
		{
			R = 255.0f;
			G = 255.0f;
			B = 255.0f;
		}
		else
		{
			if (255.0f != R)
				R += ((255.0f - R) / 100.0f) * (100.0f - sat_percent);
			if (255.0f != G)
				G += ((255.0f - G) / 100.0f) * (100.0f - sat_percent);
			if (255.0f != B)
				B += ((255.0f - B) / 100.0f) * (100.0f - sat_percent);
		}
	}

	if (100 != bri_percent)
	{
		if (0 == bri_percent)
		{
			R = 0.0f;
			G = 0.0f;
			B = 0.0f;
		}
		else
		{
			if (0.0f != R)
				R *= static_cast<float>(bri_percent) / 100.0f;
			if (0.0f != G)
				G *= static_cast<float>(bri_percent) / 100.0f;
			if (0.0f != B)
				B *= static_cast<float>(bri_percent) / 100.0f;
		}
	}

	if (R < 0.0f)
		R = 0.0f;
	else if (R > 255.0f)
		R = 255.0f;

	if (G < 0.0f)
		G = 0.0f;
	else if (G > 255.0f)
		G = 255.0f;

	if (B < 0.0f)
		B = 0.0f;
	else if (B > 255.0f)
		B = 255.0f;

	RGB rgb;

	rgb.r = static_cast<unsigned char>(R);
	rgb.g = static_cast<unsigned char>(G);
	rgb.b = static_cast<unsigned char>(B);

	return rgb;
}




// This render mode won't apply to a curved 3D space.
void draw_objects(bool disable_colouring)
{





	if (false == disable_colouring)
	{
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHT1);
		glEnable(GL_LIGHT2);
		glEnable(GL_LIGHT3);
		glEnable(GL_LIGHT4);
		glEnable(GL_LIGHT5);
	}
	else
	{
		glColor3f(0.0f, 0.0f, 0.0f);
		glDisable(GL_LIGHTING);
	}

	static const float rad_to_deg = 180.0f / static_cast<float>(pi);

	glPushMatrix();

	glTranslatef(camera_x_transform, camera_y_transform, 0);

	if (draw_curves)
	{


		for (size_t i = 0; i < pos.size(); i++)
		{
			for (size_t j = 0; j < pos[i].size() - 1; j++)
			{
				double t = j / static_cast<double>(pos[i].size() - 1);

				//set<vector_4> point_set;

				//for (size_t j = 0; j < all_4d_points[i].size(); j++)
				//{
				//	point_set.insert(all_4d_points[i][j]);
				//}

				//if (point_set.size() == all_4d_points[i].size())
				//	continue;

//				double t = static_cast<float>(point_set.size()) / static_cast<float>(all_4d_points[i].size());

				RGB rgb = HSBtoRGB(static_cast<unsigned short>(300.f * t), 75, 100);

				float colour[] = { rgb.r / 255.0f, rgb.g / 255.0f, rgb.b / 255.0f, 1.0f };

				glMaterialfv(GL_FRONT, GL_DIFFUSE, colour);

				vector_4 line = pos[i][j + 1] - pos[i][j];
				glPushMatrix();
				glTranslatef(static_cast<float>(pos[i][j].x), static_cast<float>(pos[i][j].y), 0);

				float line_len = static_cast<float>(line.length());
				line.normalize();

				float yaw = 0.0f;

				if (fabsf(static_cast<float>(line.x)) < 0.00001f && fabsf(static_cast<float>(line.z)) < 0.00001f)
					yaw = 0.0f;
				else
					yaw = atan2f(static_cast<float>(line.x), static_cast<float>(line.z));

				float pitch = -atan2f(static_cast<float>(line.y), static_cast<float>(sqrt(line.x * line.x + line.z * line.z)));

				glRotatef(yaw * rad_to_deg, 0.0f, 1.0f, 0.0f);
				glRotatef(pitch * rad_to_deg, 1.0f, 0.0f, 0.0f);

				if (j == 0)
					glutSolidSphere(0.005 * 1.5, 16, 16);

				if (j < pos[i].size() - 2)
					gluCylinder(glu_obj, 0.005, 0.005, line_len, 20, 2);
				else
					glutSolidCone(0.005 * 4, 0.005 * 8, 20, 20);

				glPopMatrix();
			}

		}

		



		//for (size_t i = 0; i < all_4d_points.size(); i++)
		//{
		//	for (size_t j = 0; j < all_4d_points[i].size() - 1; j++)
		//	{
		//		double t = j / static_cast<double>(all_4d_points[i].size() - 1);

		//		RGB rgb = HSBtoRGB(static_cast<unsigned short>(300.f * t), 75, 100);

		//		float colour[] = { rgb.r / 255.0f, rgb.g / 255.0f, rgb.b / 255.0f, 1.0f };

		//		glMaterialfv(GL_FRONT, GL_DIFFUSE, colour);

		//		vector_4 line = all_4d_points[i][j + 1] - all_4d_points[i][j];

		//		glPushMatrix();
		//		glTranslatef(static_cast<float>(all_4d_points[i][j].x), static_cast<float>(all_4d_points[i][j].y), 0);

		//		float line_len = static_cast<float>(line.length());
		//		//line.normalize();

		//		float yaw = 0.0f;

		//		if (fabsf(static_cast<float>(line.x)) < 0.00001f)
		//			yaw = 0.0f;
		//		else
		//			yaw = atan2f(static_cast<float>(line.y), 0);

		//		float pitch = -atan2f(static_cast<float>(line.y), static_cast<float>(sqrt(line.x * line.x)));

		//		glRotatef(yaw * rad_to_deg, 0.0f, 1.0f, 0.0f);
		//		glRotatef(pitch * rad_to_deg, 1.0f, 0.0f, 0.0f);

		//		if (j == 0)
		//			glutSolidSphere(0.005 * 1.5, 16, 16);

		//		if (j < all_4d_points[i].size() - 2)
		//			gluCylinder(glu_obj, 0.005, 0.005, line_len, 20, 2);
		//		else
		//			glutSolidCone(0.005 * 4, 0.005 * 8, 20, 20);
		//		
		//		glPopMatrix();
		//	}

		//}

	}










	glDisable(GL_LIGHTING);

	// If we do draw the axis at all, make sure not to draw its outline.
	if(draw_axis && false == disable_colouring)
	{
		glEnable(GL_ALPHA);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

		glLineWidth(outline_width);

		glBegin(GL_LINES);

		//glColor4f(0, 0, 0, 0.5);

		glColor4f(1, 0, 0, 0.25);
		glVertex3f(0, 0, 0);
		glVertex3f(1, 0, 0);
		glColor4f(0, 1, 0, 0.25);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 1, 0);
		glColor4f(0, 0, 1, 0.25);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 0, 1);

		//glColor4f(0, 0, 0, 0.25);
		//glVertex3f(0, 0, 0);
		//glVertex3f(-1, 0, 0);
		//glVertex3f(0, 0, 0);
		//glVertex3f(0, -1, 0);
		//glVertex3f(0, 0, 0);
		//glVertex3f(0, 0, -1);

		glEnd();

		glDisable(GL_BLEND);
		glDisable(GL_ALPHA);
	}

	glPopMatrix();
}








#endif