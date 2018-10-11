//
// Created by Dieqiao Feng on 10/3/18.
//

#ifndef STRUCTURALTABLEAUX_TOWER_H
#define STRUCTURALTABLEAUX_TOWER_H

#include <vector>
#include <string>

const int undefined = -7777777;

struct Atom {
    bool is_positive;
    int arg1, arg2;     // the depth of the corresponding argument

    Atom(bool _isp, int _a1, int _a2) : is_positive(_isp), arg1(_a1), arg2(_a2) {}

    static int cmp(const Atom &a1, const Atom &a2)
    {
        if (a1.is_positive ^ a2.is_positive)
            return a2.is_positive ? -1 : 1;
        else if (a1.arg1 == -1 || a2.arg1 == -1)
            return undefined;
        else if (a1.arg1 != a2.arg1)
            return a1.arg1 < a2.arg1 ? -1 : 1;
        else if (a1.arg2 == -1 || a2.arg2 == -1)
            return undefined;
        else if (a1.arg2 != a2.arg2)
            return a1.arg2 < a2.arg2 ? -1 : 1;
        else
            return 0;
    }

    bool operator<(const Atom &other) const
    {
        return cmp(*this, other) == -1;
    }

    void flip()
    {
        is_positive = !is_positive;
    }
};

class Tower {

public:
    int natom;              // number of total atomic formulas in the tower
    int nvar;               // number of quantified variables
    std::vector<Atom> atoms;
    std::vector<Tower> towers;

    /*
     * Compare two Towers (each Tower must be self-normalized)
     * For Atom and Atom, just compare <is_positive, arg1, arg2>
     * For Tower and Tower, compare as following ordering
     *   1) natom
     *   2) nvar
     *   3) atoms.size()
     *   4) towers.size()
     *   5) pairwise compare atoms
     *   6) pairwise compare towers
     */
    static int cmp(const Tower &, const Tower &);
    bool operator<(const Tower &other) const
    {
        return cmp(*this, other) == -1;
    }

    void sort();

    void flip();

    std::string to_skeleton_string();

    static Tower make_tower(int, int, const std::vector<Atom> &, const std::vector<Tower> &);

    static std::vector<Tower> generate(int);

private:
    typedef std::pair<int, Tower> _Record;

    Tower() = default;
    Tower(int, int, std::vector<Atom>, std::vector<Tower>);
    static void _dfs(int, int, std::vector<std::vector<_Record>> &,
                        std::vector<_Record> &);
    static std::vector<Tower> generate_skeleton(int);
};


#endif //STRUCTURALTABLEAUX_TOWER_H
