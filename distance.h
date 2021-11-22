#ifndef DISTANCE_H
#define DISTANCE_H

#include <vector>
#include "model.h"

double distance_point_vertex(double x, double y, double z, Vertex vertex);

void get_distance_vertex();

double get_distance(double x, double y);

Point multiply_matrix_point(std::vector<std::vector<double>> A, Point point);

std::vector<std::vector<double>> get_rotation_matrix(double alpha, int axis);

#endif // DISTANCE_H
