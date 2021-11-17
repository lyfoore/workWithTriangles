#ifndef MODEL_H
#define MODEL_H

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string.h>

//extern const int N_CELLS;

class Vertex
{
public:
    double m_x;
    double m_y;
    double m_z;
    int number;
    std::vector<int> m_faces;
    std::vector<int> m_edges;
    double getDistFromInit();
    Vertex(char *facet, int i);
};

class Face
{
public:
    std::vector<int> m_vertexes;
    std::vector<int> m_edges;
};

class Edge
{
public:
    std::vector<int> m_vertexes;
    std::vector<int> m_faces;
};

class Model
{
public:
    double m_xMax, m_yMax, m_zMax, m_xMin, m_yMin, m_zMin;
    double m_minSize, m_maxSize;
    unsigned long m_faceN;
    char m_faceNchar[4];
    char header[80];
    std::vector<Vertex> m_all_vertexes;
    std::vector<Face> m_all_faces;
    std::vector<Edge> m_all_edges;
    void load(char *fname);
    void getMinMax();
    void deleting_twins();
    void distribution2D();
    void getEdges();
};

double difference(Vertex a, Vertex b);


#endif // MODEL_H
