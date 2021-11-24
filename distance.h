#ifndef DISTANCE_H
#define DISTANCE_H

#include <vector>
#include "model.h"

double distance_point_vertex(double x, double y, double z, Vertex vertex);

void get_distance_vertex();

double get_distance(double x, double y, double z);

Point multiply_matrix_point(std::vector<std::vector<double>> A, Point point);

std::vector<std::vector<double>> get_rotation_matrix(double alpha, int axis);

int edge_equation(double z, double y, double Z1, double Y1, double Z2, double Y2);

double get_distance_point_line_2D(Point point, Point line1, Point line2);

double get_distance_point_point(Point point1, Point point2);

void get_distance_field();

std::vector<std::vector<double>> multiply_matrixes(std::vector<std::vector<double>> A, std::vector<std::vector<double>> B);

#endif // DISTANCE_H
