#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string.h>

#include "distance.h"
#include "model.h"

extern Model model;
const int N_CELLS = 100;
extern double distances2D[N_CELLS][N_CELLS][3];

double distance_point_vertex(double x, double y, Vertex vertex)
{
    return pow( pow( vertex.m_x - x, 2 ) + pow( vertex.m_y - y, 2 ) , 1/2 );
}

void get_distance_vertex()            // may be optimized
{
    int n_vertexes = model.m_all_vertexes.size();
    double z_plane = (model.m_zMax - model.m_zMin) / 2;
    double dist_to_vertexes[n_vertexes];
    for (int i = 0; i < N_CELLS; i++)
    {
        for (int j = 0; j < N_CELLS; j++)
        {
            distances2D[i][j][2] = distance_point_vertex(distances2D[i][j][0], distances2D[i][j][1], model.m_all_vertexes[0]);
            for (int k = 1; k < n_vertexes; k++)
            {
                if (distance_point_vertex(distances2D[i][j][0], distances2D[i][j][1], model.m_all_vertexes[k]) < distances2D[i][j][2])
                    distances2D[i][j][2] = distance_point_vertex(distances2D[i][j][0], distances2D[i][j][1], model.m_all_vertexes[k]);
            }
        }
    }
}
