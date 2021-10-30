#include "GUI.h"

#include <my_include/gl.h>
#include <my_include/glu.h>
#include <my_include/glut.h>
//#include <iostream>

void display();

void color()
{
    glClearColor(0.0, 0.0, 0.0, 1.0);
}

void reshape(int, int);

GUI::GUI()
{

}

void GUI::init(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);

    glutInitWindowPosition(200, 100);
    glutInitWindowSize(500, 500);

    glutCreateWindow("Window 1");

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);

    color();

    glutMainLoop();
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();

//    draw
    glBegin(GL_TRIANGLES);

    glVertex2f(0.0, 5.0);
    glVertex2f(4.0, -3.0);
    glVertex2f(-4.0, -3.0);


    glEnd();

    glutSwapBuffers();

}

void reshape(int w, int h)
{
//    viewport
    glViewport(0, 0, w, h);
//    projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(-10, 10, -10, 10);
    glMatrixMode(GL_MODELVIEW);

}
