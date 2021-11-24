#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string.h>

#include "model.h"
#include "display/GUI.h"
#include "distance.h"

const double PRECIZE = 0.01 / 100; //0.01%
double SIZE;
const int N_CELLS = 100;
extern double distances2D[N_CELLS][N_CELLS][3];

void Model::load(char *fname)
{
    // creating binary file object for reading
    std::ifstream file(fname, std::ios::in | std::ios::binary);
    // reading header from STL
    file.read(header, 80);
    std::cout << "Header - " << header << std::endl;
    // reading number of faces from STL
    file.read(m_faceNchar, 4);
    // converting char to unsigned long
    m_faceN = *(unsigned long *)m_faceNchar;
    std::cout << "Number of faces - " << m_faceN << std::endl;
    // creating buffer for reading face by face (12 4-bytes numbers & 2-byte "attribute byte count")
    char buffer[50];
    for (int i = 0; i < m_faceN; i++)
    {
        // writing to the buffer
        file.read(buffer, 50);
        // getting vertixes coordinates
        Vertex v1(buffer + 12, i);
        Vertex v2(buffer + 24, i);
        Vertex v3(buffer + 36, i);

        // adding vertixes to the all_vertixes vector
        m_all_vertexes.push_back(v1);
        m_all_vertexes.push_back(v2);
        m_all_vertexes.push_back(v3);

        // adding indexes of vertexes to the face
        Face triangle;
        triangle.m_vertexes.push_back(3 * i);
        triangle.m_vertexes.push_back(3 * i + 1);
        triangle.m_vertexes.push_back(3 * i + 2);

        // adding face to the all_faces vector
        m_all_faces.push_back(triangle);
    }
    file.close();
}

void Model::getMinMax() {

    m_xMax = 1e-10;     m_xMin = 1e10;
    m_yMax = 1e-10;     m_yMin = 1e10;
    m_zMax = 1e-10;     m_zMin = 1e10;

    for (unsigned long i = 0; i < m_faceN * 3; i++) {
        if (m_all_vertexes[i].m_x < m_xMin) m_xMin = m_all_vertexes[i].m_x;
        if (m_all_vertexes[i].m_y < m_yMin) m_yMin = m_all_vertexes[i].m_y;
        if (m_all_vertexes[i].m_z < m_zMin) m_zMin = m_all_vertexes[i].m_z;

        if (m_all_vertexes[i].m_x > m_xMax) m_xMax = m_all_vertexes[i].m_x;
        if (m_all_vertexes[i].m_y > m_yMax) m_yMax = m_all_vertexes[i].m_y;
        if (m_all_vertexes[i].m_z > m_zMax) m_zMax = m_all_vertexes[i].m_z;
    }
    // finding the minimum between three dimensionals
    m_minSize = 1e10;
    for (int i = 0; i < 3; i++) {
        if ((m_xMax - m_xMin) < m_minSize) m_minSize = m_xMax - m_xMin;
        if ((m_yMax - m_yMin) < m_minSize) m_minSize = m_yMax - m_yMin;
        if ((m_zMax - m_zMin) < m_minSize) m_minSize = m_zMax - m_zMin;
    }

    // finding the maximun between three dimensionals
    m_maxSize = 1e-10;
    for (int i = 0; i < 3; i++) {
        if ((m_xMax - m_xMin) > m_maxSize) m_maxSize = m_xMax - m_xMin;
        if ((m_yMax - m_yMin) > m_maxSize) m_maxSize = m_yMax - m_yMin;
        if ((m_zMax - m_zMin) > m_maxSize) m_maxSize = m_zMax - m_zMin;
    }
}

void Model::deleting_twins()
{
    for (int i = m_all_vertexes.size() - 1; i > 0; i--)
    {
//        std::cout << (i * 100) / m_all_vertexes.size() << std::endl;
        for (int j = i - 1; j >= 0; j--)
        {
            if (difference(m_all_vertexes[i], m_all_vertexes[j]) < m_minSize * PRECIZE)
            {
                // editing faces of vertexes
                for(int k = 0; k < m_all_vertexes[i].m_faces.size(); k++)
                {
                    m_all_vertexes[j].m_faces.push_back(m_all_vertexes[i].m_faces[k]);

                    // editing vertexes of faces
                    for (int l = 0; l < 3; l++)
                    {
                        if (m_all_faces[m_all_vertexes[i].m_faces[k]].m_vertexes[l] == i)
                        {
                            m_all_faces[m_all_vertexes[i].m_faces[k]].m_vertexes[l] = j;
                            break;
                        }
                    }
                }

                // deleting element
                m_all_vertexes.erase(m_all_vertexes.begin() + i);

                // fixing trouble with deleting (transition indexes of vertexes by -1 after deleting)
                for (int m = 0; m < m_all_faces.size(); m++)
                {
                    for (int n = 0; n < 3; n++)
                    {
                        if (m_all_faces[m].m_vertexes[n] > i)
                            m_all_faces[m].m_vertexes[n] -= 1;
                    }
                }

                break;
            }
        }
    }
}

double Vertex::getDistFromInit()
{
    return sqrt(m_x * m_x + m_y * m_y + m_z * m_z);
}

Vertex::Vertex(char *facet, int i)
{
    char x[4] = {facet[0], facet[1], facet[2], facet[3]};
    char y[4] = {facet[4], facet[5], facet[6], facet[7]};
    char z[4] = {facet[8], facet[9], facet[10], facet[11]};

    float xx = *(float *)x;
    float yy = *(float *)y;
    float zz = *(float *)z;

    m_x = (double)xx;
    m_y = (double)yy;
    m_z = (double)zz;
    m_faces.push_back(i);
}


double difference(Vertex a, Vertex b)
{
    return sqrt(pow(b.m_x - a.m_x, 2) +
                pow(b.m_y - a.m_y, 2) +
                pow(b.m_z - a.m_z, 2));
}


void Model::distribution2D()
{
//    double distr[N_CELLS][N_CELLS][3]; // [3] == [x, y, distance]
    for (int i = 0; i < N_CELLS; i++) {
        for (int j = 0; j < N_CELLS; j++) {
            distances2D[i][j][0] = m_xMin + ((m_xMax - m_xMin) * i) / N_CELLS;
            distances2D[i][j][1] = m_yMin + ((m_yMax - m_yMin) * j) / N_CELLS;
//            std::cout << distances2D[i][j][0] << ' ' << distances2D[i][j][1] << std::endl;
        }
    }
}


bool consistsIn(int item, std::vector<int> vector)
{
    for (int i = 0; i < vector.size(); i++)
    {
        if (item == vector[i])
            return 1;
    }
    return 0;
}


void Model::getEdges()
{
    // go over the all vertixes
    for (int i = 0; i < m_all_vertexes.size(); i++)
    {
        // make temporary vector for each vertex to hold another vertexes for avoiding repetitions
        std::vector<int> temp;
        // go over all triangles that has "i" vertex
        for (int j = 0; j < m_all_vertexes[i].m_faces.size(); j++)
        {
            // go over the all vertixes inside of triangle
            for (int k = 0; k < 3; k++)
            {
                // ignore all vertexes with index that smaller than "i" for avoiding repetitions (all edges for vertexes with indexes < "i" already exists)
                if (m_all_faces[ m_all_vertexes[i].m_faces[j] ].m_vertexes[k] > i && !(consistsIn(m_all_faces[ m_all_vertexes[i].m_faces[j] ].m_vertexes[k], temp)))
                    temp.push_back(m_all_faces[ m_all_vertexes[i].m_faces[j] ].m_vertexes[k]);
            }
        }
        if (!(temp.empty()))
        {
            // go over the all temp == go over the all edges for "i" vertex
            for (int l = 0; l < temp.size(); l++)
            {
                Edge edge;
                // connect edge with vertex indexes
                edge.m_vertexes.push_back(i);
                edge.m_vertexes.push_back(temp[l]);

                m_all_edges.push_back(edge);

                // finding triangles with this edge
                for (int m = 0; m < m_all_vertexes[i].m_faces.size(); m++)
                {
                    // if temp[l] consists in the triangle with "i" vertex
                    if (consistsIn(temp[l], m_all_faces[ m_all_vertexes[i].m_faces[m] ].m_vertexes))
                    {
                        // connect edge with faces indexes
                        m_all_edges[m_all_edges.size() - 1].m_faces.push_back(m_all_vertexes[i].m_faces[m]);

                        // connect triangle with edge indexes
                        m_all_faces[ m_all_vertexes[i].m_faces[m] ].m_edges.push_back(m_all_edges.size() - 1);

                        // connect vertexes with edge indexes
                        m_all_vertexes[i].m_edges.push_back(m_all_edges.size() - 1);
                        m_all_vertexes[temp[l]].m_edges.push_back(m_all_edges.size() - 1);
                    }
                }
            }
        }
    }
}


void Model::init_matrixes()
{
    double alpha;
    double gamma;
    double beta;

    // points for normals in vertexes
    Point p01, p02, p10, p12, p20, p21; // first number - vertex index, second number - neighbour index

    double k02;
    double k12;

    for (int i = 0; i < m_all_faces.size(); i++)
    {
        Point point0;
        Point point1;
        Point point2;
        std::vector<std::vector<double>> rot_matrix;

        point0 = {m_all_vertexes[m_all_faces[i].m_vertexes[0]].m_x,
                  m_all_vertexes[m_all_faces[i].m_vertexes[0]].m_y,
                  m_all_vertexes[m_all_faces[i].m_vertexes[0]].m_z};

        point1 = {m_all_vertexes[m_all_faces[i].m_vertexes[1]].m_x,
                  m_all_vertexes[m_all_faces[i].m_vertexes[1]].m_y,
                  m_all_vertexes[m_all_faces[i].m_vertexes[1]].m_z};

        point2 = {m_all_vertexes[m_all_faces[i].m_vertexes[2]].m_x,
                  m_all_vertexes[m_all_faces[i].m_vertexes[2]].m_y,
                  m_all_vertexes[m_all_faces[i].m_vertexes[2]].m_z};

        // calc the angle between two points and Ox
        alpha = atan((point1.m_y - point0.m_y)/(point1.m_z - point0.m_z));
        if (point1.m_z < point0.m_z) alpha += M_PI;

        // rotate around Ox
        Point point0_new =  {point0.m_x - point0.m_x, point0.m_y - point0.m_y, point0.m_z - point0.m_z};
        Point point1_new =  {point1.m_x - point0.m_x, point1.m_y - point0.m_y, point1.m_z - point0.m_z};
        Point point2_new =  {point2.m_x - point0.m_x, point2.m_y - point0.m_y, point2.m_z - point0.m_z};

        point0_new = multiply_matrix_point( get_rotation_matrix(alpha, 1), point0_new );
        point1_new = multiply_matrix_point( get_rotation_matrix(alpha, 1), point1_new );
        point2_new = multiply_matrix_point( get_rotation_matrix(alpha, 1), point2_new );

        // calc the angle between two points and Oy
        gamma = atan((point1_new.m_x - point0_new.m_x)/(point1_new.m_z - point0_new.m_z));
        gamma *= -1;

        // rotate around Oy
        point0_new = multiply_matrix_point( get_rotation_matrix(gamma, 2), point0_new );
        point1_new = multiply_matrix_point( get_rotation_matrix(gamma, 2), point1_new );
        point2_new = multiply_matrix_point( get_rotation_matrix(gamma, 2), point2_new );

        // calc the angle between two points and Oz
        beta = atan((point2_new.m_x - point1_new.m_x)/(point2_new.m_y - point1_new.m_y));
        if (point2_new.m_y < 0) beta += M_PI;

        // rotate around Oz
        point0_new = multiply_matrix_point( get_rotation_matrix(beta, 3), point0_new );
        point1_new = multiply_matrix_point( get_rotation_matrix(beta, 3), point1_new );
        point2_new = multiply_matrix_point( get_rotation_matrix(beta, 3), point2_new );

        // get full rotation matrix
        rot_matrix = multiply_matrixes(get_rotation_matrix(gamma, 2), get_rotation_matrix(alpha, 1));
        rot_matrix = multiply_matrixes(get_rotation_matrix(beta, 3), rot_matrix);

        m_all_faces[i].m_matrix = rot_matrix;

        point0_new =  {point0.m_x - point0.m_x, point0.m_y - point0.m_y, point0.m_z - point0.m_z};
        point1_new =  {point1.m_x - point0.m_x, point1.m_y - point0.m_y, point1.m_z - point0.m_z};
        point2_new =  {point2.m_x - point0.m_x, point2.m_y - point0.m_y, point2.m_z - point0.m_z};

        point0_new = multiply_matrix_point(rot_matrix, point0_new);
        point1_new = multiply_matrix_point(rot_matrix, point1_new);
        point2_new = multiply_matrix_point(rot_matrix, point2_new);


        m_all_faces[i].m_point0_new = point0_new;
        m_all_faces[i].m_point1_new = point1_new;
        m_all_faces[i].m_point2_new = point2_new;

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

//        std::vector<std::vector<double>> m_all_faces[i].lines (9, std::vector<double> (4, 0.));

        m_all_faces[i].lines = {{point0_new.m_y, point0_new.m_z, point2_new.m_y, point2_new.m_z},  // 0
                               {point2_new.m_y, point2_new.m_z, point1_new.m_y, point1_new.m_z},   // 1
                               {point1_new.m_y, point1_new.m_z, point0_new.m_y, point0_new.m_z},   // 2
                               {point0_new.m_y, point0_new.m_z, p01.m_y, p01.m_z},         // 3
                               {point0_new.m_y, point0_new.m_z, p02.m_y, p02.m_z},         // 4
                               {point2_new.m_y, point2_new.m_z, p21.m_y, p21.m_z},         // 5
                               {point2_new.m_y, point2_new.m_z, p20.m_y, p20.m_z},         // 6
                               {point1_new.m_y, point1_new.m_z, p10.m_y, p10.m_z},         // 7
                               {point1_new.m_y, point1_new.m_z, p12.m_y, p12.m_z}};        // 8
    }

}
