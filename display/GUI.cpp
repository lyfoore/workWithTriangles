#include "GUI.h"
#include "main.h"

#include <my_include/gl.h>
#include <my_include/glu.h>
#include <my_include/glut.h>
#include <iostream>

extern Model model;

void display();

void color()
{
    glClearColor(0.0, 0.0, 0.0, 1.0);
}

void reshape(int, int);

void timer(int);

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
    glutTimerFunc(0, timer, 0);

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
    glTranslatef(0., 0., -90.);
    glRotatef(angle, 1., 1., 1.);

//    draw

    for (int i = 0; i < model.m_all_faces.size(); i++)
    {
        glBegin(GL_TRIANGLES);
        for (int j = 0; j < 3; j++)
        {
            glVertex3d(model.m_all_vertexes[ model.m_all_faces[i].m_vertexes[j] ].m_x,
                       model.m_all_vertexes[ model.m_all_faces[i].m_vertexes[j] ].m_y,
                       model.m_all_vertexes[ model.m_all_faces[i].m_vertexes[j] ].m_z);
        }
        glEnd();
    }

//    glBegin(GL_TRIANGLES);

//    glColor3f(1., 0., 0.);
//    glVertex3d(0., 1., 0);
//    glVertex3d(-1., -1., 0.);
//    glVertex3d(1., -1., 0.);

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
