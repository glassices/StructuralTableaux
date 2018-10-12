/*
 * This project is to implement my idea about innex normal form and
 * structural tableaux. Hopefully it can beat traditional tableaux and
 * resolution algorithm
 */
#include <iostream>
#include "Tower.h"

using namespace std;

int main() {
    std::cout << "Hello, World!" << std::endl;

    std::vector<std::vector<Record>> res;
    Tower::search_raw(2, res);
    for (auto &e : res) cout << e.size() << ' ';
    cout << endl;

    return 0;
}