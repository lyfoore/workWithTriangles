#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <vector>
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
    double m_xMax, m_yMax, m_zMax, m_xMin, m_yMin, m_zMin;
    double m_minSize;
    unsigned long m_faceN;
    char m_faceNchar[4];
    char header[80];
    std::vector<Vertex> m_all_vertexes;
    std::vector<Face> m_all_faces;
    void load(char *fname);
    void getMinMax();
    void deleting_twins();
    void distribution2D();
};

double difference(Vertex a, Vertex b);


#endif // MODEL_H
