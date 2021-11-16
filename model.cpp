#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string.h>

#include "model.h"
#include "display/GUI.h"

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


void Model::distribution2D() {
//    double distr[N_CELLS][N_CELLS][3]; // [3] == [x, y, distance]
    for (int i = 0; i < N_CELLS; i++) {
        for (int j = 0; j < N_CELLS; j++) {
            distances2D[i][j][0] = m_xMin + ((m_xMax - m_xMin) * i) / N_CELLS;
            distances2D[i][j][1] = m_yMin + ((m_yMax - m_yMin) * j) / N_CELLS;
//            std::cout << distances2D[i][j][0] << ' ' << distances2D[i][j][1] << std::endl;
        }
    }
}
