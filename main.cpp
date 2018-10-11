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

    auto sks = Tower::generate(2);
    for (auto &e : sks)
        cout << e.to_skeleton_string() << endl;

    return 0;
}