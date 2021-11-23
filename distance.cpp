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
    return sqrt( pow(point1.m_x - point2.m_x, 2) +  pow(point1.m_y - point2.m_y, 2) + pow(point1.m_z - point2.m_z, 2));
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

int edge_equation(double z, double y, double Z1, double Y1, double Z2, double Y2)
{
    // E(x, y) = (z - Z)dY - (y - Y)dZ
    double dZ = Z2 - Z1;
    double dY = Y2 - Y1;
    double res = (z - Z1) * dY - (y - Y1) * dZ;
//    std::cout << res << std::endl;

    if (res > 0) return 1; // on the right side
    if (res < 0) return -1; // on the left side
    if (res == 0) return 0; // on line
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
    Point point0;
    Point point1;
    Point point2;
    double alpha;
    Point point0_new;
    Point point1_new;
    Point point2_new;
    Point xyz_new;
    double gamma;
    double beta;

    // points for normals in vertexes
    Point p01, p02, p10, p12, p20, p21; // first number - vertex index, second number - neighbour index

    double k02;
    double k12;
    double dist_temp;

    // go over the all triangles
    for (int i = 0; i < model.m_all_faces.size(); i++)
    {
        point0 = {model.m_all_vertexes[model.m_all_faces[i].m_vertexes[0]].m_x,
              model.m_all_vertexes[model.m_all_faces[i].m_vertexes[0]].m_y,
              model.m_all_vertexes[model.m_all_faces[i].m_vertexes[0]].m_z};

        point1 = {model.m_all_vertexes[model.m_all_faces[i].m_vertexes[1]].m_x,
              model.m_all_vertexes[model.m_all_faces[i].m_vertexes[1]].m_y,
              model.m_all_vertexes[model.m_all_faces[i].m_vertexes[1]].m_z};

        point2 = {model.m_all_vertexes[model.m_all_faces[i].m_vertexes[2]].m_x,
              model.m_all_vertexes[model.m_all_faces[i].m_vertexes[2]].m_y,
              model.m_all_vertexes[model.m_all_faces[i].m_vertexes[2]].m_z};

        // calc the angle between two points and Ox
        alpha = atan((point1.m_y - point0.m_y)/(point1.m_z - point0.m_z));
        if (point1.m_z < point0.m_z) alpha += M_PI;

        // rotate around Ox
        point0_new =  {point0.m_x - point0.m_x, point0.m_y - point0.m_y, point0.m_z - point0.m_z};
        point1_new =  {point1.m_x - point0.m_x, point1.m_y - point0.m_y, point1.m_z - point0.m_z};
        point2_new =  {point2.m_x - point0.m_x, point2.m_y - point0.m_y, point2.m_z - point0.m_z};
        xyz_new = {x - point0.m_x, y - point0.m_y, z - point0.m_z};

        point0_new = multiply_matrix_point( get_rotation_matrix(alpha, 1), point0_new );
        point1_new = multiply_matrix_point( get_rotation_matrix(alpha, 1), point1_new );
        point2_new = multiply_matrix_point( get_rotation_matrix(alpha, 1), point2_new );
        xyz_new = multiply_matrix_point( get_rotation_matrix(alpha, 1), xyz_new );

        // calc the angle between two points and Oy
        gamma = atan((point1_new.m_x - point0_new.m_x)/(point1_new.m_z - point0_new.m_z));
        gamma *= -1;

        // rotate around Oy
        point0_new = multiply_matrix_point( get_rotation_matrix(gamma, 2), point0_new );
        point1_new = multiply_matrix_point( get_rotation_matrix(gamma, 2), point1_new );
        point2_new = multiply_matrix_point( get_rotation_matrix(gamma, 2), point2_new );
        xyz_new = multiply_matrix_point( get_rotation_matrix(gamma, 2), xyz_new );

        // calc the angle between two points and Oz
        beta = atan((point2_new.m_x - point1_new.m_x)/(point2_new.m_y - point1_new.m_y));
        if (point2_new.m_y < 0) beta += M_PI;

        // rotate around Oz
        point0_new = multiply_matrix_point( get_rotation_matrix(beta, 3), point0_new );
        point1_new = multiply_matrix_point( get_rotation_matrix(beta, 3), point1_new );
        point2_new = multiply_matrix_point( get_rotation_matrix(beta, 3), point2_new );
        xyz_new = multiply_matrix_point( get_rotation_matrix(beta, 3), xyz_new );


        // get normal equations
        k02 = (point2_new.m_y - point0_new.m_y)/(point2_new.m_z - point0_new.m_z);
        k02 = 1/k02 * (-1);

        p02 = {0., -0.5 * k02, -0.5};
        p20 = {0., (point2_new.m_z - 0.5) * k02 + point2_new.m_y - k02 * point2_new.m_z, point2_new.m_z - 0.5};

        k12 = (point2_new.m_y - point1_new.m_y)/(point2_new.m_z - point1_new.m_z);
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
        for (int i = 0; i < 9; i++)
        {
            results.push_back(edge_equation(xyz_new.m_z, xyz_new.m_y, lines[i][1], lines[i][0], lines[i][3], lines[i][2]));
    //            std::cout << results[i] << std::endl;
        }


        if (results[0] > 0 && results[1] > 0 && results[2] > 0) dist_temp = abs(xyz_new.m_x);
        if (results[0] < 0 && results[4] > 0 && results[6] < 0) dist_temp = sqrt(pow(get_distance_point_line_2D(xyz_new, point0_new, point2_new), 2) + pow(xyz_new.m_x, 2));
        if (results[1] < 0 && results[5] > 0 && results[8] < 0) dist_temp = sqrt(pow(get_distance_point_line_2D(xyz_new, point2_new, point1_new), 2) + pow(xyz_new.m_x, 2));
        if (results[2] < 0 && results[7] > 0 && results[3] < 0) dist_temp = sqrt(pow(get_distance_point_line_2D(xyz_new, point1_new, point0_new), 2) + pow(xyz_new.m_x, 2));
        if (results[4] < 0 && results[3] > 0) dist_temp = get_distance_point_point(point0_new, xyz_new);
        if (results[5] < 0 && results[6] > 0) dist_temp = get_distance_point_point(point2_new, xyz_new);
        if (results[7] < 0 && results[8] > 0) dist_temp = get_distance_point_point(point1_new, xyz_new);

        if(dist_temp < dist) dist = dist_temp;
    }
    return dist;
}
