//
// Created by Dieqiao Feng on 10/3/18.
//

#include "Tower.h"
#include <cassert>
#include <algorithm>

struct Segment
{
    int lo, hi, dep;
    std::vector<Segment> ch;

    Segment(int _lo, int _hi, int _dep) : lo(_lo), hi(_hi), dep(_dep) {}
};

int Tower::cmp(const Tower &t1, const Tower &t2)
{
    if (t1.natom != t2.natom) return t1.natom < t2.natom ? -1 : 1;
    if (t1.nvar != t2.nvar) return t1.nvar < t2.nvar ? -1 : 1;
    if (t1.atoms.size() != t2.atoms.size())
        return t1.atoms.size() < t2.atoms.size() ? -1 : 1;
    if (t1.towers.size() != t2.towers.size())
        return t1.towers.size() < t2.towers.size() ? -1 : 1;

    for (int i = 0; i < t1.atoms.size(); i++) {
        int ret = Atom::cmp(t1.atoms[i], t2.atoms[i]);
        if (ret != 0) return ret; // include undefined
    }

    for (int i = 0; i < t1.towers.size(); i++) {
        int ret = Tower::cmp(t1.towers[i], t2.towers[i]);
        if (ret != 0) return ret; // include undefined
    }

    return 0;
}

void Tower::sort()
{
    std::sort(atoms.begin(), atoms.end());
    std::sort(towers.begin(), towers.end());
}

void Tower::flip()
{
    for (auto &e : atoms) e.flip();
    for (auto &e : towers) e.flip();
    sort();
}

std::string Tower::to_skeleton_string()
{
    std::vector<std::string> subs;
    for (auto &e : atoms) subs.emplace_back(e.is_positive ? "+" : "-");
    for (auto &e : towers) subs.push_back(e.to_skeleton_string());
    std::string res = std::to_string(nvar);

    if (subs.size() > 1) res.push_back('(');
    for (int i = 0; i < subs.size(); i++) {
        if (i > 0) res.push_back('|');
        res.append(subs[i]);
    }
    if (subs.size() > 1) res.push_back(')');
    return res;
}

Tower::Tower(int _natom, int _nvar,
             std::vector<Atom> _atoms,
             std::vector<Tower> _towers)
    : natom(_natom), nvar(_nvar), atoms(std::move(_atoms)), towers(std::move(_towers))
{}

Tower Tower::make_tower(int natom, int nvar,
                        const std::vector<Atom> &atoms,
                        const std::vector<Tower> &towers)
{
    Tower tower(natom, nvar, atoms, towers);
    tower.sort();
    return tower;
}

void Tower::_dfs(int n, int a, std::vector<std::vector<Tower::_Record>> &data,
                               std::vector<Tower::_Record> &dummy)
{
    int m = 0;
    for (auto &e : dummy) m += e.second.natom;

    if (a + m == n) {
        /*
         * maximum quantifers: 1 + sig_i(s_i - 1)
         * maximum slots: sig_i(s_i - 1) - nq + 1
         * dummy is in descending order according to tower
         */
        std::vector<Atom> atoms((unsigned)a, Atom(false, -1, -1));
        std::vector<Tower> towers;
        int max_nq = a + 1; // 1 + sig_i(s_i - 1)
        for (auto it = dummy.rbegin(); it != dummy.rend(); ++it) {
            towers.push_back(it->second);
            max_nq += it->first - 1;
        }

        // enumerate polarity of atomic formulas
        for (int nfalse = 0; nfalse <= a; nfalse++) {
            for (int i = 0; i < nfalse; i++) atoms[i].is_positive = false;
            for (int i = nfalse; i < a; i++) atoms[i].is_positive = true;

            // enumerate the number of quantifiers
            for (int nq = 1; nq <= max_nq; nq++) {
                Tower tower(n, nq, atoms, towers);
                data[n].emplace_back(max_nq - nq, tower);
            }
        }
        return;
    }

    int limit = dummy.empty() ? n - a - m : std::min(dummy.back().second.natom, n - a - m);
    for (int k = 1; k <= limit; k++)
        for (auto &record : data[k]) {
            if (!dummy.empty() && Tower::cmp(record.second, dummy.back().second) == 1)
                continue;
            dummy.push_back(record);
            _dfs(n, a, data, dummy);
            dummy.pop_back();
        }
}

std::vector<Tower> Tower::generate_skeleton(int n)
{
    std::vector<std::vector<_Record>> data;
    data.emplace_back(); // initiate data[0]
    std::vector<_Record> dummy;
    for (int k = 1; k <= n; k++) {
        printf("%d\n", k);
        /*
         * Deal with the case that at least one atom is under the quantifiers
         */
        data.emplace_back(); // initiate data[k]

        for (int a = 1; a <= k; a++) _dfs(k, a, data, dummy);

        for (int i = 0; i < data[n].size(); i++) {
            const Tower &tower = data[n][i].second;
            for (int nq = 1; nq <= data[n][i].first; nq++)
                data[n].emplace_back(data[n][i].first - nq,
                                     Tower(n, nq, std::vector<Atom>(),
                                           std::vector<Tower>(1, tower)));
        }
    }

    std::vector<Tower> towers;
    for (auto &e : data[n]) towers.push_back(e.second);
    return towers;
}

std::vector<Tower> Tower::generate(int n)
{
    /*
     * First step is to generate the skeletons of towers.
     * A skeleton is a structural tower with all arguments of Atom
     * set to -1. The whole skeleton should obey ordering requirement
     */
    std::vector<Tower> sks = generate_skeleton(n);
    return sks;
}
