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

    auto towers = Tower::search_raw(2);
    for (auto &tower : towers) cout << tower.to_string() << endl;
    cout << towers.size() << endl;

    return 0;
}