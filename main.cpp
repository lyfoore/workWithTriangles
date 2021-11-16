#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string.h>

#include "model.h"
#include "display/GUI.h"
#include "distance.h"

// using namespace std;

GUI gui;
Model model;
const int N_CELLS = 100;
double distances2D[N_CELLS][N_CELLS][3];  // [3] == [x, y, distance]

int main(int argc, char **argv)
{
    char dir[1024],file_mod[1024];
    sprintf(dir,"%s",__FILE__);
    dir[strlen(dir)-9]=0; //hack to get source directory
    sprintf(file_mod,"%s/stl_files/teapot.stl",dir);

    model.load(file_mod);
    model.getMinMax();
    std::cout << "x_min - " << model.m_xMin << std::endl;
    std::cout << "before deleting twins - " << model.m_all_vertexes.size() << std::endl;
    model.deleting_twins();
    std::cout << "after deleting twins - " << model.m_all_vertexes.size() << std::endl;

    get_distance_vertex();

    gui.init(argc,argv);

    return 0;
}
