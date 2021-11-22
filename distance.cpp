#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string.h>

#include "model.h"
#include "distance.h"


extern Model model;
const int N_CELLS = 100;
extern double distances2D[N_CELLS][N_CELLS][3];
extern const double Z_PLANE;

double distance_point_vertex(double x, double y, double z, Vertex vertex)
{
    return sqrt( pow( vertex.m_x - x, 2 ) + pow( vertex.m_y - y, 2 ) + pow( vertex.m_z - z, 2 ));
}

void get_distance_vertex()            // may be optimized
{
    int n_vertexes = model.m_all_vertexes.size();
    double dist_to_vertexes[n_vertexes];
    for (int i = 0; i < N_CELLS; i++)
    {
        for (int j = 0; j < N_CELLS; j++)
        {
            distances2D[i][j][2] = distance_point_vertex(distances2D[i][j][0], distances2D[i][j][1], Z_PLANE, model.m_all_vertexes[0]);
            for (int k = 1; k < n_vertexes; k++)
            {
//                std::cout << distance_point_vertex(distances2D[i][j][0], distances2D[i][j][1], model.m_all_vertexes[k]) << std::endl;
                if (distance_point_vertex(distances2D[i][j][0], distances2D[i][j][1], Z_PLANE, model.m_all_vertexes[k]) < distances2D[i][j][2])
                    distances2D[i][j][2] = distance_point_vertex(distances2D[i][j][0], distances2D[i][j][1], Z_PLANE, model.m_all_vertexes[k]);
            }
        }
    }
}


Point multiply_matrix_point(std::vector<std::vector<double>> A, Point point)
{
    Point temp;

    temp.m_x = A[0][0] * point.m_x + A[0][1] * point.m_y + A[0][2] * point.m_z;
    temp.m_y = A[1][0] * point.m_x + A[1][1] * point.m_y + A[1][2] * point.m_z;
    temp.m_z = A[2][0] * point.m_x + A[2][1] * point.m_y + A[2][2] * point.m_z;

    return temp;
}


std::vector<std::vector<double>> get_rotation_matrix(double alpha, int axis)
{
    std::vector<std::vector<double>> temp (3, std::vector<double> (3, 0.));

    switch (axis) {
    case 1:
        temp[0][0] = 1;  temp[0][1] = 0;           temp[0][2] = 0;
        temp[1][0] = 0;  temp[1][1] = cos(alpha);  temp[1][2] = (-1)*sin(alpha);
        temp[2][0] = 0;  temp[2][1] = sin(alpha);  temp[2][2] = cos(alpha);
        break;
    case 2:
        temp[0][0] = cos(alpha);       temp[0][1] = 0;  temp[0][2] = sin(alpha);
        temp[1][0] = 0;                temp[1][1] = 1;  temp[1][2] = 0;
        temp[2][0] = (-1)*sin(alpha);  temp[2][1] = 0;  temp[2][2] = cos(alpha);
        break;

    case 3:
        temp[0][0] = cos(alpha);  temp[0][1] = (-1)*sin(alpha);  temp[0][2] = 0;
        temp[1][0] = sin(alpha);  temp[1][1] = cos(alpha);       temp[1][2] = 0;
        temp[2][0] = 0;           temp[2][1] = 0;                temp[2][2] = 1;
        break;
    }

    return temp;
}


double get_distance(double x, double y)
{
    // go over at all triangles
    for (int i = 0; i < model.m_all_faces.size(); i++)
    {
        Point p1 = {model.m_all_vertexes[model.m_all_faces[i].m_vertexes[0]].m_x,
                    model.m_all_vertexes[model.m_all_faces[i].m_vertexes[0]].m_y,
                    model.m_all_vertexes[model.m_all_faces[i].m_vertexes[0]].m_z};

        // get mean of two ponts of triangle
        Point p2 = {(model.m_all_vertexes[model.m_all_faces[i].m_vertexes[1]].m_x + model.m_all_vertexes[model.m_all_faces[i].m_vertexes[2]].m_x)/2,
                    (model.m_all_vertexes[model.m_all_faces[i].m_vertexes[1]].m_y + model.m_all_vertexes[model.m_all_faces[i].m_vertexes[2]].m_y)/2,
                    (model.m_all_vertexes[model.m_all_faces[i].m_vertexes[1]].m_z + model.m_all_vertexes[model.m_all_faces[i].m_vertexes[2]].m_z)/2};

        // calc the angle between two points and Oy

        // rotate around Oy

        // calc the angle between two points and Oz

        // rotate around Oz

        // calc the angle between two points and Ox

        // rotate around Ox
    }

}
