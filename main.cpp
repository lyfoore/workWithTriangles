#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

const double PRESIZE = 0.01;
double SIZE;

class Point {
public:
    double m_x;
    double m_y;
    double m_z;
    int number;
    double getDistFromInit() {
        return sqrt(m_x*m_x + m_y*m_y + m_z*m_z);
    };
};

bool comparison(Point a, Point b) {
    return (a.getDistFromInit() < b.getDistFromInit());
}

double difference(Point a, Point b) {
    return sqrt(pow(b.m_x - a.m_x, 2) + 
                pow(b.m_y - a.m_y, 2) +
                pow(b.m_z - a.m_z, 2));
}

int main() {
    Point a1 {1, 1, 1, 1};
    Point a2 {10, 4, 3, 2};
    Point a3 {9, 5, 2, 3};
    Point a4 {0, 4, 5, 4};
    Point a5 {10, 0, 3, 5};
    int n = 5;
    std::vector<Point> vec = {a1, a2, a3, a4, a5};
    for (auto now : vec) {
        std::cout << now.number << " ";
    }
    std::cout << std::endl;
    std::sort(vec.begin(), vec.end(), comparison);
    for (auto now : vec) {
        std::cout << now.number << " ";
    }
    std::cout << std::endl;
    for (auto now : vec) {
        std::cout << now.getDistFromInit() << " ";
    }
    SIZE = vec[vec.size() - 1].getDistFromInit() * 2;
    return 0;
}
