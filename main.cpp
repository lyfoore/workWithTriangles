#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string.h>

#include "main.h"

// using namespace std;

const double PRECIZE = 0.01 / 100; //0.01%
double SIZE;

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

        all_vertexes.push_back(v1);
        all_vertexes.push_back(v2);
        all_vertexes.push_back(v3);

        Face triangle;
        triangle.m_vertexes.push_back(3 * i);
        triangle.m_vertexes.push_back(3 * i + 1);
        triangle.m_vertexes.push_back(3 * i + 2);
    }
    file.close();
};

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

bool comparison(Vertex a, Vertex b)
{
    return (a.getDistFromInit() < b.getDistFromInit());
}

double difference(Vertex a, Vertex b)
{
    return sqrt(pow(b.m_x - a.m_x, 2) +
                pow(b.m_y - a.m_y, 2) +
                pow(b.m_z - a.m_z, 2));
}

std::vector<Vertex> merging(std::vector<Vertex> vector)
{
    for (int i = 0; i < vector.size(); i++)
    {
        for (int j = i; j < vector.size(); j++)
        {
            if (difference(vector[i], vector[j]) < SIZE * PRECIZE)
            {
                return vector; // не доделал
            }
        }
    }
}

int main()
{
    // Point a1{1, 1, 1, 1};
    // Point a2{10, 4, 3, 2};
    // Point a3{9, 5, 2, 3};
    // Point b1{0, 4, 5, 4};
    // Point b2{10, 0, 3, 5};
    // Point b3{6, 4, 6, 6};
    // int n = 5;
    // std::vector<Point> vec = {a1, a2, a3, b1, b2, b3};
    // for (auto now : vec)
    // {
    //     std::cout << now.number << " ";
    // }
    // std::cout << std::endl;
    // std::sort(vec.begin(), vec.end(), comparison);
    // for (auto now : vec)
    // {
    //     std::cout << now.number << " ";
    // }
    // std::cout << std::endl;
    // for (auto now : vec)
    // {
    //     std::cout << now.getDistFromInit() << " ";
    // }
    // SIZE = vec[vec.size() - 1].getDistFromInit() * 2;
    Model model;
    char *path = "D:\\projects\\workWithTriangles\\stl_files\\dodeca_half_a.stl";
    model.load(path);

    return 0;
}
