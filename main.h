#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string.h>

class Vertex
{
public:
    double m_x;
    double m_y;
    double m_z;
    int number;
    std::vector<int> m_faces;
    double getDistFromInit();
    Vertex(char *facet, int i);
};

class Face
{
public:
    std::vector<int> m_vertexes;
};

class Model
{
public:
    unsigned long m_faceN;
    char m_faceNchar[4];
    char header[80];
    std::vector<Vertex> all_vertexes;
    std::vector<Face> all_faces;
    void load(char *fname);
};



