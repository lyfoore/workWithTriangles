#include "GUI.h"
#include "model.h"

#include <my_include/gl.h>
#include <my_include/glu.h>
#include <my_include/glut.h>
#include <iostream>

extern Model model;
const int N_CELLS = 100;
extern double distances2D[N_CELLS][N_CELLS][3];
extern const double Z_PLANE;

void display();

void draw_distance();

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
    glutPostRedisplay();
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

    for (int i = 0; i < model.m_all_faces.size(); i++)
    {
        glColor3f(0.7, 0.7, 0.7);
        glBegin(GL_TRIANGLES);
        for (int j = 0; j < 3; j++)
        {
            glVertex3d(model.m_all_vertexes[ model.m_all_faces[i].m_vertexes[j] ].m_x,
                       model.m_all_vertexes[ model.m_all_faces[i].m_vertexes[j] ].m_y,
                       model.m_all_vertexes[ model.m_all_faces[i].m_vertexes[j] ].m_z);
        }
        glEnd();
    }

    draw_distance();

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
    glBegin(GL_POINTS);
    for (int i = 0; i < N_CELLS; i++)
    {
        for (int j = 0; j < N_CELLS; j++)
        {
            glColor3f(1. * distances2D[i][j][2] / model.m_minSize * 5, 0., 0.);
            glVertex3f(distances2D[i][j][0], distances2D[i][j][1], Z_PLANE);
        }
    }
    glEnd();
}
