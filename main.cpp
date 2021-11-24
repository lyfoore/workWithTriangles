#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string.h>
#include <time.h>

#include "model.h"
#include "display/GUI.h"
#include "distance.h"

// using namespace std;

GUI gui;
Model model;
const int N_CELLS = 100;
double Z_PLANE;
double distances2D[N_CELLS][N_CELLS][3];  // [3] == [x, y, distance]

int main(int argc, char **argv)
{
    time_t start, end;
    time(&start);

    char dir[1024],file_mod[1024];
    sprintf(dir,"%s",__FILE__);
    dir[strlen(dir)-9]=0; //hack to get source directory
    sprintf(file_mod,"%s/stl_files/dodeca_half_a.stl",dir);

    model.load(file_mod);
    model.getMinMax();
    std::cout << "x_min - " << model.m_xMin << std::endl;
    std::cout << "before deleting twins - " << model.m_all_vertexes.size() << std::endl;
    model.deleting_twins();
    std::cout << "after deleting twins - " << model.m_all_vertexes.size() << std::endl;
    model.init_matrixes();
    std::cout << "init finished" << std:: endl;
    model.getEdges();
    std::cout << "number of edges - "<< model.m_all_edges.size() << std:: endl;

//    for (int l = 0; l < model.m_all_edges.size(); l++)
//    {
//        std::cout << model.m_all_edges[l].m_vertexes[0] << ' ' << model.m_all_edges[l].m_vertexes[1] << std::endl;
//    }

    Z_PLANE = (model.m_zMax - model.m_zMin) / 2;
    model.distribution2D();
    get_distance_field();
//    get_distance_vertex();

//    for (int i = 0; i < N_CELLS; i++)
//    {
//        for (int j = 0; j < N_CELLS; j++)
//        {
//            std::cout << distances2D[i][j][0] << ' ' << distances2D[i][j][1] << ' ' <<  distances2D[i][j][2] << std::endl;
//        }
//    }

//    double dist = get_distance(22.3233, 28.7363, Z_PLANE);

//    std::cout << "dist = " << dist << std::endl;

//    std::vector<std::vector<double>> temp;

    time(&end);
    double seconds = difftime(end, start);
    std::cout << seconds << " seconds" << std::endl;

//    std::vector<std::vector<double>> resulttt = multiply_matrixes(get_rotation_matrix(3.1415, 1), get_rotation_matrix(3.1415, 1));

//    std::cout << resulttt[0][0] << " " << resulttt[0][1] << " " << resulttt[0][2] << std::endl;
//    std::cout << resulttt[1][0] << " " << resulttt[1][1] << " " << resulttt[1][2] << std::endl;
//    std::cout << resulttt[2][0] << " " << resulttt[2][1] << " " << resulttt[2][2] << std::endl;

    gui.init(argc,argv);

    return 0;
}
