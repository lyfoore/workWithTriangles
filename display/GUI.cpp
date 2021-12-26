#include "GUI.h"
#include "model.h"
#include "distance.h"

#include <cmath>
#include <my_include/gl.h>
#include <my_include/glu.h>
#include <my_include/glut.h>
#include <iostream>

extern Model model;
const int N_CELLS = 100;
double Z_PLANE;
double distances2D[N_CELLS+1][N_CELLS+1][3];  // [3] == [x, y, distance]

void display();

void draw_distance();

void draw_model_by_edges();

void draw_triangle();
void draw_origin();

void color()
{
    glClearColor(0.0, 0.0, 0.0, 1.0);
}

void reshape(int, int);

void timer(int);

float rx = 0;
float ry = 0;
float rx0 = 0;
float ry0 = 0;
float mx0 = 0;
float my0 = 0;
float scale_k = 1;
float blur_k = 4;

float mouse_x=0.0;
float mouse_y=0.0;
int rotate=0;


//mouse_move
void m_m(int x,int y) //mouse move
{
    if (rotate==1)
    {
        rx=rx0+0.5*(x-mx0);
        ry=ry0+0.5*(y-my0);
    }
    glutPostRedisplay();
}


//mouse down
void m_d(int button, int state, int x, int y)  //mouse down
{

    if (state==GLUT_UP)
    {
        rotate=0;
        rx0=rx;
        ry0=ry;
    }
    if (state==GLUT_DOWN)
    {
        rotate=1;
        mx0=x;
        my0=y;

    }

//    double W_WIDTH=W_HEIGHT*(w_x1-w_x0)/(w_y1-w_y0);
//    mouse_x=(1.0*x)/W_WIDTH;

//    mouse_y=(W_HEIGHT-(1.0*y))/W_HEIGHT;

    glutPostRedisplay();
}

void keyboardFunc (int key, int x, int y)
{
    if (key == GLUT_KEY_UP)
        // zoom in
        scale_k += 0.1;
    if (key == GLUT_KEY_DOWN)
        // zoom out
        scale_k -= 0.1;
    if (key == GLUT_KEY_LEFT)
        blur_k += 0.1;
    if (key == GLUT_KEY_RIGHT)
        blur_k -= 0.1;

    glutPostRedisplay();
}


void getDistribution2D()
{
    double x_size = model.m_xMax - model.m_xMin;
    double y_size = model.m_yMax - model.m_yMin;
    for (int i = 0; i < N_CELLS; i++) {
        for (int j = 0; j < N_CELLS; j++) {
            distances2D[i][j][0] = model.m_xMin + (x_size * i) / N_CELLS;
            distances2D[i][j][1] = model.m_yMin + (y_size * j) / N_CELLS;
        }
    }
}

//int shetchik = 0;

void GUI::drawDistance()
{
    Z_PLANE = (model.m_zMax - model.m_zMin) / 2;
    getDistribution2D();
    for (int i = 0; i < N_CELLS; i++)
    {
        for (int j = 0; j < N_CELLS; j++)
        {
            distances2D[i][j][2] = get_distance(distances2D[i][j][0], distances2D[i][j][1], Z_PLANE);
        }
    }
//    std::cout << shetchik << std::endl;
}


GUI::GUI()
{

}

void GUI::init(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);

    glutInitWindowPosition(200, 100);
    glutInitWindowSize(500, 500);

    glutCreateWindow("Window 1");

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
//    glutTimerFunc(0, timer, 0);
    glutMotionFunc(m_m);
    glutMouseFunc(m_d);
    glutSpecialFunc(keyboardFunc);

    glEnable(GL_DEPTH_TEST);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    color();

    glutMainLoop();
}

//float x_pos = -10.;
//float state = 1.;
float angle = 0.;

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    // translation
    glTranslatef(10., 0., -110.);
    glScalef(scale_k, scale_k, scale_k);
    glRotatef(rx, 0., 1., 0.);
    glRotatef(ry, 1., 0., 0.);

//    draw

//    for (int i = 0; i < model.m_all_faces.size(); i++)
//    {
//        glColor3f(0.7, 0.7, 0.7);
//        glBegin(GL_TRIANGLES);
//        for (int j = 0; j < 3; j++)
//        {
//            glVertex3d(model.m_all_vertexes[ model.m_all_faces[i].m_vertexes[j] ].m_x,
//                       model.m_all_vertexes[ model.m_all_faces[i].m_vertexes[j] ].m_y,
//                       model.m_all_vertexes[ model.m_all_faces[i].m_vertexes[j] ].m_z);
//        }
//        glEnd();
//    }

    draw_model_by_edges();

    draw_distance();

    draw_origin();

//    glColor3f(0.7, 0.7, 0.7);
//    glBegin(GL_TRIANGLES);

//    glVertex3d(model.m_all_vertexes[ model.m_all_faces[10].m_vertexes[0] ].m_x,
//               model.m_all_vertexes[ model.m_all_faces[10].m_vertexes[0] ].m_y,
//               model.m_all_vertexes[ model.m_all_faces[10].m_vertexes[0] ].m_z);
//    glVertex3d(model.m_all_vertexes[ model.m_all_faces[10].m_vertexes[1] ].m_x,
//               model.m_all_vertexes[ model.m_all_faces[10].m_vertexes[1] ].m_y,
//               model.m_all_vertexes[ model.m_all_faces[10].m_vertexes[1] ].m_z);
//    glVertex3d(model.m_all_vertexes[ model.m_all_faces[10].m_vertexes[2] ].m_x,
//               model.m_all_vertexes[ model.m_all_faces[10].m_vertexes[2] ].m_y,
//               model.m_all_vertexes[ model.m_all_faces[10].m_vertexes[2] ].m_z);

//    glVertex3d(model.m_all_faces[10].m_point0_new.m_x,
//               model.m_all_faces[10].m_point0_new.m_y,
//               model.m_all_faces[10].m_point0_new.m_z);
//    glVertex3d(model.m_all_faces[10].m_point1_new.m_x,
//               model.m_all_faces[10].m_point1_new.m_y,
//               model.m_all_faces[10].m_point1_new.m_z);
//    glVertex3d(model.m_all_faces[10].m_point2_new.m_x,
//               model.m_all_faces[10].m_point2_new.m_y,
//               model.m_all_faces[10].m_point2_new.m_z);
//    glVertex3d(0.,
//               model.m_all_faces[10].lines[3][0],
//               model.m_all_faces[10].lines[3][1]);
//    glVertex3d(0.,
//               model.m_all_faces[10].lines[3][2],
//               model.m_all_faces[10].lines[3][3]);
//    glVertex3d(0.,
//               model.m_all_faces[10].lines[3][2],
//               model.m_all_faces[10].lines[3][3]);

//    glEnd();

//    draw_triangle();

//    glBegin(GL_TRIANGLES);

//    glColor3f(10., 0., 0.);
//    glVertex3d(0., 10., 0);
//    glVertex3d(-10., -10., 0.);
//    glVertex3d(-10., -10., 0.);

//    glEnd();

    glutSwapBuffers();
}

void reshape(int w, int h)
{
//    viewport
    glViewport(0, 0, w, h);
//    projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60, 1, 2, 150);
    glMatrixMode(GL_MODELVIEW);

}

void timer(int)
{
    glutPostRedisplay();
    glutTimerFunc(1000./60., timer, 0);

    angle += 1.;
    if (angle >= 360.)
        angle -= 360.;

//    if (x_pos >= -9. || x_pos <= -20.)
//        state *= -1;
//    x_pos += state * 0.2;
}

void draw_distance()
{
//    glBegin(GL_POINTS);
//    for (int i = 0; i < N_CELLS; i++)
//    {
//        for (int j = 0; j < N_CELLS; j++)
//        {
//            glColor3d(1. * distances2D[i][j][2] / model.m_minSize * blur_k, 0., 0.);
//            glVertex3f(distances2D[i][j][0], distances2D[i][j][1], Z_PLANE);
//        }
//    }
//    glEnd();


    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    for (int i = 0; i < N_CELLS - 1 ; i++)
    {
        glBegin(GL_TRIANGLE_STRIP);
        for (int j = 0; j < N_CELLS; j++)
        {
            glColor3d(1. * distances2D[i][j][2] / model.m_minSize * blur_k, 0., 0.);
            glVertex3f(distances2D[i][j][0], distances2D[i][j][1], Z_PLANE);
            glColor3d(1. * distances2D[i+1][j][2] / model.m_minSize * blur_k, 0., 0.);
            glVertex3f(distances2D[i+1][j][0], distances2D[i+1][j][1], Z_PLANE);
        }
        glEnd();
    }

}

void draw_model_by_edges()
{
    glColor3f(0.7, 0.7, 0.7);
    glBegin(GL_LINES);
    for (int i = 0; i < model.m_all_edges.size(); i++)
    {
        for (int j = 0; j < 2; j++)
        {
            glVertex3f( model.m_all_vertexes[model.m_all_edges[i].m_vertexes[j]].m_x,
                        model.m_all_vertexes[model.m_all_edges[i].m_vertexes[j]].m_y,
                        model.m_all_vertexes[model.m_all_edges[i].m_vertexes[j]].m_z);
        }
    }
    glEnd();
}


void draw_triangle()
{
    int i = 0;
    std::vector<std::vector<double>> rot_matrix;
    Point point0 = {model.m_all_vertexes[model.m_all_faces[i].m_vertexes[0]].m_x,
                model.m_all_vertexes[model.m_all_faces[i].m_vertexes[0]].m_y,
                model.m_all_vertexes[model.m_all_faces[i].m_vertexes[0]].m_z};

    Point point1 = {model.m_all_vertexes[model.m_all_faces[i].m_vertexes[1]].m_x,
                model.m_all_vertexes[model.m_all_faces[i].m_vertexes[1]].m_y,
                model.m_all_vertexes[model.m_all_faces[i].m_vertexes[1]].m_z};

    Point point2 = {model.m_all_vertexes[model.m_all_faces[i].m_vertexes[2]].m_x,
                model.m_all_vertexes[model.m_all_faces[i].m_vertexes[2]].m_y,
                model.m_all_vertexes[model.m_all_faces[i].m_vertexes[2]].m_z};
    Point xyz = {0., 1., 0.};

    // calc the angle between two points and Ox
    double alpha = atan((point1.m_y - point0.m_y)/(point1.m_z - point0.m_z));
    if (point1.m_z < point0.m_z) alpha += M_PI;

    // rotate around Ox
    Point point0_new =  {point0.m_x - point0.m_x, point0.m_y - point0.m_y, point0.m_z - point0.m_z};
    Point point1_new =  {point1.m_x - point0.m_x, point1.m_y - point0.m_y, point1.m_z - point0.m_z};
    Point point2_new =  {point2.m_x - point0.m_x, point2.m_y - point0.m_y, point2.m_z - point0.m_z};
    Point xyz_new = {xyz.m_x - point0.m_x, xyz.m_y - point0.m_y, xyz.m_z - point0.m_z};

    point0_new = multiply_matrix_point( get_rotation_matrix(alpha, 1), point0_new );
    point1_new = multiply_matrix_point( get_rotation_matrix(alpha, 1), point1_new );
    point2_new = multiply_matrix_point( get_rotation_matrix(alpha, 1), point2_new );
    xyz_new = multiply_matrix_point( get_rotation_matrix(alpha, 1), xyz_new );

    // calc the angle between two points and Oy
    double gamma = atan((point1_new.m_x - point0_new.m_x)/(point1_new.m_z - point0_new.m_z));
    gamma *= -1;

    // rotate around Oy
    point0_new = multiply_matrix_point( get_rotation_matrix(gamma, 2), point0_new );
    point1_new = multiply_matrix_point( get_rotation_matrix(gamma, 2), point1_new );
    point2_new = multiply_matrix_point( get_rotation_matrix(gamma, 2), point2_new );
    xyz_new = multiply_matrix_point( get_rotation_matrix(gamma, 2), xyz_new );

    // calc the angle between two points and Oz
    double beta = atan((point2_new.m_x - point1_new.m_x)/(point2_new.m_y - point1_new.m_y));
    if (point2_new.m_y < 0) beta += M_PI;

    // rotate around Oz
    point0_new = multiply_matrix_point( get_rotation_matrix(beta, 3), point0_new );
    point1_new = multiply_matrix_point( get_rotation_matrix(beta, 3), point1_new );
    point2_new = multiply_matrix_point( get_rotation_matrix(beta, 3), point2_new );
    xyz_new = multiply_matrix_point( get_rotation_matrix(beta, 3), xyz_new );

    rot_matrix = multiply_matrixes(get_rotation_matrix(gamma, 2), get_rotation_matrix(alpha, 1));
    rot_matrix = multiply_matrixes(get_rotation_matrix(beta, 3), rot_matrix);

    point0_new =  {point0.m_x - point0.m_x, point0.m_y - point0.m_y, point0.m_z - point0.m_z};
    point1_new =  {point1.m_x - point0.m_x, point1.m_y - point0.m_y, point1.m_z - point0.m_z};
    point2_new =  {point2.m_x - point0.m_x, point2.m_y - point0.m_y, point2.m_z - point0.m_z};
    xyz_new = {xyz.m_x - point0.m_x, xyz.m_y - point0.m_y, xyz.m_z - point0.m_z};

    point0_new = multiply_matrix_point( rot_matrix, point0_new );
    point1_new = multiply_matrix_point( rot_matrix, point1_new );
    point2_new = multiply_matrix_point( rot_matrix, point2_new );
    xyz_new = multiply_matrix_point( rot_matrix, xyz_new );

    // points for normals in vertexes
    Point p01, p02, p10, p12, p20, p21; // first number - vertex index, second number - neighbour index

    // get normal equations
    double k02 = (point2_new.m_y - point0_new.m_y)/(point2_new.m_z - point0_new.m_z);
    k02 = 1/k02 * (-1);

    p02 = {0., -0.5 * k02, -0.5};
    p20 = {0., (point2_new.m_z - 0.5) * k02 + point2_new.m_y - k02 * point2_new.m_z, point2_new.m_z - 0.5};

    double k12 = (point2_new.m_y - point1_new.m_y)/(point2_new.m_z - point1_new.m_z);
    k12 = 1/k12 * (-1);

    p21 = {0., (point2_new.m_z + 0.5) * k12 + point2_new.m_y - k12 * point2_new.m_z, point2_new.m_z + 0.5};
    p12 = {0., (point1_new.m_z + 0.5) * k12 + point1_new.m_y - k12 * point1_new.m_z, point1_new.m_z + 0.5};

    // for z-axis k = 0
    p01 = {0., point0_new.m_y - 0.5, point0_new.m_z};
    p10 = {0., point1_new.m_y - 0.5, point1_new.m_z};

    double lines[9][4] = {{point0_new.m_y, point0_new.m_z, point2_new.m_y, point2_new.m_z}, // 0
                          {point2_new.m_y, point2_new.m_z, point1_new.m_y, point1_new.m_z}, // 1
                          {point1_new.m_y, point1_new.m_z, point0_new.m_y, point0_new.m_z}, // 2
                          {point0_new.m_y, point0_new.m_z, p01.m_y, p01.m_z},               // 3
                          {point0_new.m_y, point0_new.m_z, p02.m_y, p02.m_z},               // 4
                          {point2_new.m_y, point2_new.m_z, p21.m_y, p21.m_z},               // 5
                          {point2_new.m_y, point2_new.m_z, p20.m_y, p20.m_z},               // 6
                          {point1_new.m_y, point1_new.m_z, p10.m_y, p10.m_z},               // 7
                          {point1_new.m_y, point1_new.m_z, p12.m_y, p12.m_z}};              // 8

    std::vector<int> results;
    for (int j = 0; j < 9; j++)
    {
        results.push_back(edge_equation(xyz_new.m_z, xyz_new.m_y, lines[j][1], lines[j][0], lines[j][3], lines[j][2]));
    }

    if (results[0] > 0 && results[1] > 0 && results[2] > 0) std::cout << abs(xyz_new.m_x) << std::endl;
    if (results[0] < 0 && results[4] > 0 && results[6] < 0) std::cout << sqrt(pow(get_distance_point_line_2D(xyz_new, point0_new, point2_new), 2) + pow(xyz_new.m_x, 2)) << std::endl;
    if (results[1] < 0 && results[5] > 0 && results[8] < 0) std::cout << sqrt(pow(get_distance_point_line_2D(xyz_new, point2_new, point1_new), 2) + pow(xyz_new.m_x, 2)) << std::endl;
    if (results[2] < 0 && results[7] > 0 && results[3] < 0) std::cout << sqrt(pow(get_distance_point_line_2D(xyz_new, point1_new, point0_new), 2) + pow(xyz_new.m_x, 2)) << std::endl;
    if (results[4] < 0 && results[3] > 0) std::cout << get_distance_point_point(point0_new, xyz_new) << std::endl;
    if (results[5] < 0 && results[6] > 0) std::cout << get_distance_point_point(point2_new, xyz_new) << std::endl;
    if (results[7] < 0 && results[8] > 0) std::cout << get_distance_point_point(point1_new, xyz_new) << std::endl;


    glColor3f(0.7, 0.7, 0.7);
    glBegin(GL_POINTS);
    glVertex3d(xyz.m_x, xyz.m_y, xyz.m_z);
    glVertex3d(xyz_new.m_x, xyz_new.m_y, xyz_new.m_z);
    glEnd();
    glBegin(GL_LINES);
    glVertex3d(point0_new.m_x, point0_new.m_y, point0_new.m_z);
    glVertex3d(p02.m_x, p02.m_y, p02.m_z);
    glVertex3d(point2_new.m_x, point2_new.m_y, point2_new.m_z);
    glVertex3d(p20.m_x, p20.m_y, p20.m_z);
    glVertex3d(point2_new.m_x, point2_new.m_y, point2_new.m_z);
    glVertex3d(p21.m_x, p21.m_y, p21.m_z);
    glVertex3d(point1_new.m_x, point1_new.m_y, point1_new.m_z);
    glVertex3d(p12.m_x, p12.m_y, p12.m_z);
    glVertex3d(point0_new.m_x, point0_new.m_y, point0_new.m_z);
    glVertex3d(p01.m_x, p01.m_y, p01.m_z);
    glVertex3d(point1_new.m_x, point1_new.m_y, point1_new.m_z);
    glVertex3d(p10.m_x, p10.m_y, p10.m_z);
    glEnd();


    glColor3f(0.7, 0.7, 0.7);
    glBegin(GL_TRIANGLES);
    glVertex3f(point0_new.m_x, point0_new.m_y, point0_new.m_z);
    glVertex3f(point1_new.m_x, point1_new.m_y, point1_new.m_z);
    glVertex3f(point2_new.m_x, point2_new.m_y, point2_new.m_z);

    glVertex3f(point0.m_x, point0.m_y, point0.m_z);
    glVertex3f(point1.m_x, point1.m_y, point1.m_z);
    glVertex3f(point2.m_x, point2.m_y, point2.m_z);
    glEnd();
}

void draw_origin()
{

    glBegin(GL_LINES);
    glColor3f(1., 0., 0.);
    glVertex3f(0., 0., 0.);
    glVertex3f(1., 0., 0.);
    glColor3f(0., 1., 0.);
    glVertex3f(0., 0., 0.);
    glVertex3f(0., 1., 0.);
    glColor3f(0., 0., 1.);
    glVertex3f(0., 0., 0.);
    glVertex3f(0., 0., 1.);
    glEnd();
}
