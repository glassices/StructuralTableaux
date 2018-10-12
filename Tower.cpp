//
// Created by Dieqiao Feng on 10/3/18.
//

#include "Tower.h"
#include <cassert>
#include <algorithm>

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
    //assert(!subs.empty());
    assert(!atoms.empty() || !towers.empty());
    std::string res = std::to_string(nvar);
    res.push_back(':');

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

void Tower::_dfs1(int n, int a, std::vector<std::vector<Tower::_Record>> &data,
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

    int limit = std::min(n - a - m, n - 1);
    if (!dummy.empty()) limit = std::min(limit, dummy.back().second.natom);
    for (int k = 1; k <= limit; k++)
        for (auto &record : data[k]) if (record.first > 0) {
            if (!dummy.empty() && Tower::cmp(record.second, dummy.back().second) == 1)
                continue;
            dummy.push_back(record);
            _dfs1(n, a, data, dummy);
            dummy.pop_back();
        }
}

std::vector<Tower> Tower::generate_skeleton(int n)
{
    std::vector<std::vector<_Record>> data;
    data.emplace_back(); // initiate data[0]
    std::vector<_Record> dummy;
    for (int k = 1; k <= n; k++) {
        /*
         * Deal with the towers that contain at least two children
         */
        data.emplace_back(); // initiate data[k]

        for (int a = 0; a <= k; a++) _dfs1(k, a, data, dummy);

        for (int i = 0; i < data[k].size(); i++) {
            std::vector<Tower> towers(1, data[k][i].second);
            for (int nq = 1; nq <= data[k][i].first; nq++)
                data[k].emplace_back(data[k][i].first - nq,
                                     Tower(k, nq, std::vector<Atom>(), towers));
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

int _cmp(int l1, int r1, int l2, int r2, std::vector<Record> &sk)
{
    if (sk[l1].is_atom()) {
        if (sk[l2].is_atom()) {
            if (sk[l1].is_positive() ^ sk[l2].is_positive())
                return sk[l2].is_positive() ? -1 : 1;
            if (sk[l1].arg1() != sk[l2].arg1())
                return sk[l1].arg1() < sk[l2].arg1() ? -1 : 1;
            if (sk[l1].arg2() != sk[l2].arg2())
                return sk[l1].arg2() < sk[l2].arg2() ? -1 : 1;
            return 0;
        }
        else
            return -1;
    }
    else {
        if (sk[l2].is_atom())
            return 1;
        else {
            int n1 = 0, n2 = 0, prev1 = 0, prev2 = 0, last1 = -1, last2 = -1;
            std::vector<int> a1, a2;
            std::vector<std::pair<int,int>> t1, t2;
            for (int i = l1 + 1; i < r1; i++) {
                if (sk[i].is_atom()) {
                    n1++;
                    if (prev1 == 0) a1.push_back(i);
                }
                else if (sk[i].is_open()) {
                    if (prev1++ == 0) last1 = i;
                }
                else if (--prev1 == 0) {
                    t1.emplace_back(last1, i);
                }
            }
            for (int i = l2 + 1; i < r2; i++) {
                if (sk[i].is_atom()) {
                    n2++;
                    if (prev2 == 0) a2.push_back(i);
                }
                else if (sk[i].is_open()) {
                    if (prev2++ == 0) last2 = i;
                }
                else if (--prev2 == 0) {
                    t2.emplace_back(last2, i);
                }
            }
            if (n1 != n2) return n1 < n2 ? -1 : 1;
            if (sk[l1].nq() != sk[l2].nq())
                return sk[l1].nq() < sk[l2].nq() ? -1 : 1;
            if (a1.size() != a2.size())
                return a1.size() < a2.size() ? -1 : 1;
            if (t1.size() != t2.size())
                return t1.size() < t2.size() ? -1 : 1;
            for (int i = 0; i < a1.size(); i++) {
                int ret = _cmp(a1[i], a1[i], a2[i], a2[i], sk);
                if (ret != 0) return ret;
            }
            for (int i = 0; i < t1.size(); i++) {
                int ret = _cmp(t1[i].first, t1[i].second, t2[i].first, t2[i].second, sk);
                if (ret != 0) return ret;
            }
            return 0;
        }
    }
}

int _find(int k, std::vector<int> &p)
{
    return p[k] == k ? k : p[k] = _find(p[k], p);
}

void _union(int x, int y, std::vector<int> &p)
{
    p[_find(x, p)] = _find(y, p);
}

void Tower::dfs(int n, std::vector<Record> &sk)
{
    /*
     * Calculate some basic information. The reason we don't
     * put this information in arguments is that calculation
     * is cheap and we want to avoid complicated recursion
     * arguments
     * All closed braches should be ordered itself and ordered
     * with its preceding brothers
     */
    int natoms = 0;
    std::vector<int> nqs;
    std::vector<int> link(sk.size(), -1);
    std::vector<int> pos_left;
    for (int i = 0; i < sk.size(); i++) {
        if (sk[i].is_open()) {
            nqs.push_back(sk[i].nq());
            pos_left.push_back(i);
        }
        else if (sk[i].is_close()) {
            nqs.pop_back();
            link[pos_left.back()] = i;
            link[i] = pos_left.back();
            pos_left.pop_back();
        }
        else
            natoms++;
    }
    assert(!nqs.empty());
    int dep = 0, nq = nqs.back();
    for (int i = 0; i < nqs.size() - 1; i++) dep += nqs[i];

    /*
     * TODO: brach pruning here
     */

    /*
     * Case I: Open another bracket
     * Theoretically we can have maximum n - atoms + 1 new
     * quantifiers here. However we still need to leave a slot
     * for ancestor to come in, so this extra 1 should be removed
     */
    for (int k = 1; k <= n - natoms; k++) {
        sk.emplace_back(k);
        dfs(n, sk);
        sk.pop_back();
    }

    /*
     * Case II: Close one bracket. Need to check ordering now.
     * Note: There is no need to check ordering inside. Only need
     * to check ordering with its preceding brothers
     * Check quantifers are well distributed in subtrees
     */
    if ((nqs.size() > 1 || natoms == n) && !sk.back().is_open()) {
        auto k = (int)sk.size();
        sk.emplace_back();
        link[pos_left.back()] = (int)link.size();
        link.push_back(pos_left.back());

        if (link[k] == 0 || !sk[link[k]-1].is_close() ||
            _cmp(link[link[k]-1], link[k]-1, link[k], k, sk) <= 0) {
            // check distribution of quantifiers
            // [dep, dep + nq)
            std::vector<int> p(nq);
            for (int i = 0; i < nq; i++) p[i] = i;
            int prev = 0, last = -1;
            bool ok = true;
            for (int i = link[k] + 1; i < k; i++) {
                if (sk[i].is_atom()) {
                    if (prev == 0) {
                        /* one of arguments must be in [dep, dep + nq) is
                         * guaranteed when it was created
                         */
                        if (sk[i].arg1() >= dep && sk[i].arg2() >= dep)
                            _union(sk[i].arg1() - dep, sk[i].arg2() - dep, p);
                    }
                    else {
                        if (sk[i].arg1() >= dep) {
                            if (last != -1) _union(last, sk[i].arg1() - dep);
                            last = sk[i].arg1() - dep;
                        }
                        if (sk[i].arg2() >= dep) {
                            if (last != -1) _union(last, sk[i].arg2() - dep);
                            last = sk[i].arg2() - dep;
                        }
                    }
                }
                else if (sk[i].is_open()) {
                    if (prev++ == 0) last = -1;
                }
                else if (--prev == 0 && last == -1) {
                    // we have no quantifiers in this subtree
                    ok = false;
                    break;
                }
            }
            for (int i = 1; i < nq; i++) if (_find(0, p) != _find(i, p)) ok = false;
            if (ok) dfs(n, sk);
        }

        sk.pop_back();
    }

    /*
     * Case III: Create one new atomic formula
     */


}
