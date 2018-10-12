//
// Created by Dieqiao Feng on 10/3/18.
//

#ifndef STRUCTURALTABLEAUX_TOWER_H
#define STRUCTURALTABLEAUX_TOWER_H

#include <vector>
#include <string>
#include <cassert>

const int undefined = -7777777;

struct Record {
    // Close bracket
    Record() : is_bracket(true), is_open_pos(false), int1(0), int2(0) {}
    // Open bracket
    Record(int _nq) : is_bracket(true), is_open_pos(true), int1(_nq), int2(0) {}
    // Atomic formula
    Record(bool _is_pos, int _arg1, int _arg2) : is_bracket(false), is_open_pos(_is_pos), int1(_arg1), int2(_arg2) {}

    bool is_close() { return is_bracket && !is_open_pos; }

    bool is_open() { return is_bracket && is_open_pos; }

    bool is_atom() { return !is_bracket; }

    int nq() { assert(is_open()); return int1; }

    int arg1() { assert(!is_bracket); return int1; }

    int arg2() { assert(!is_bracket); return int2; }

    bool is_positive() { assert(!is_bracket); return is_open_pos; }

private:
    bool is_bracket;
    bool is_open_pos;
    int int1, int2;
};

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

    static Tower make_tower(int, int, const std::vector<Atom> &, const std::vector<Tower> &);

    static void search_raw(int, std::vector<std::vector<Record>> &);

private:
    typedef std::pair<int, Tower> _Record;

    Tower(int, int, std::vector<Atom>, std::vector<Tower>);

    static void dfs(int, std::vector<Record> &, std::vector<std::vector<Record>> &);

};


#endif //STRUCTURALTABLEAUX_TOWER_H
