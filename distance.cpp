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

double get_distance_point_point(Point point1, Point point2)
{
    return pow(point1.m_x - point2.m_x, 2) +  pow(point1.m_y - point2.m_y, 2) + pow(point1.m_z - point2.m_z, 2);
}

void get_distance_vertex()            // may be optimized
{
    int n_vertexes = model.m_all_vertexes.size();
//    double dist_to_vertexes[n_vertexes];
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

void get_distance_field()
{
    for (int i = 0; i < N_CELLS; i++)
    {
        for (int j = 0; j < N_CELLS; j++)
        {
            distances2D[i][j][2] = get_distance(distances2D[i][j][0], distances2D[i][j][1], Z_PLANE);
//            std::cout << (N_CELLS * i + j) * 100. / (N_CELLS*N_CELLS) << std::endl;
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


std::vector<std::vector<double>> multiply_matrixes(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B)
{
    std::vector<std::vector<double>> temp (3, std::vector<double> (3, 0.));

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                temp[i][j] += A[i][k] * B[k][j];
            }
        }
    }

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

double edge_equation(double z, double y, double Z1, double Y1, double Z2, double Y2)
{
    // E(x, y) = (z - Z)dY - (y - Y)dZ
    double dZ = Z2 - Z1;
    double dY = Y2 - Y1;
    double l = sqrt(dZ*dZ + dY*dY);
    dZ /= l;
    dY /= l;
    double res = (z - Z1) * dY - (y - Y1) * dZ;

    return res;
}

double get_distance_point_line_2D(Point point, Point line1, Point line2)
{
    // a*z + b*y + c = 0
    double a = (line2.m_y - line1.m_y)/(line2.m_z - line1.m_z);
    double c = line1.m_y - a * line1.m_z;

    return abs(a * point.m_z + point.m_y + c)/sqrt(1. + a*a);
}

double get_distance(double x, double y, double z)
{
    double dist = 1e15;
    double dist_temp;
    Point xyz = {x, y, z};
    Point xyz_new;

    // go over the all triangles
    for (int i = 0; i < model.m_all_faces.size(); i++)
    {
        xyz_new = {x - model.m_all_vertexes[model.m_all_faces[i].m_vertexes[0]].m_x, y - model.m_all_vertexes[model.m_all_faces[i].m_vertexes[0]].m_y, z - model.m_all_vertexes[model.m_all_faces[i].m_vertexes[0]].m_z};
        xyz_new = multiply_matrix_point(model.m_all_faces[i].m_matrix, xyz_new);
        double results[9];
        for (int j = 0; j < 9; j++)
        {
            results[j] = edge_equation(xyz_new.m_z, xyz_new.m_y, model.m_all_faces[i].lines[j][1], model.m_all_faces[i].lines[j][0], model.m_all_faces[i].lines[j][3], model.m_all_faces[i].lines[j][2]);
        }

        if (results[0] > 0 && results[1] > 0 && results[2] > 0) dist_temp = xyz_new.m_x * xyz_new.m_x;
        else if (results[0] < 0 && results[4] > 0 && results[6] < 0) dist_temp = pow(results[0], 2) + pow(xyz_new.m_x, 2);
        else if (results[1] < 0 && results[5] > 0 && results[8] < 0) dist_temp = pow(results[1], 2) + pow(xyz_new.m_x, 2);
        else if (results[2] < 0 && results[7] > 0 && results[3] < 0) dist_temp = pow(results[2], 2) + pow(xyz_new.m_x, 2);
        else if (results[4] < 0 && results[3] > 0) dist_temp = get_distance_point_point(model.m_all_faces[i].m_point0_new, xyz_new);
        else if (results[5] < 0 && results[6] > 0) dist_temp = get_distance_point_point(model.m_all_faces[i].m_point2_new, xyz_new);
        else if (results[7] < 0 && results[8] > 0) dist_temp = get_distance_point_point(model.m_all_faces[i].m_point1_new, xyz_new);

        if (dist_temp < dist) dist = dist_temp;
    }
    return sqrt(dist);
}
