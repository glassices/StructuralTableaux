//
// Created by Dieqiao Feng on 10/3/18.
//

#include "Tower.h"
#include <cassert>
#include <algorithm>
#include <sstream>

std::string Tower::to_string(int dep)
{
    assert(!atoms.empty() || !towers.empty());

    std::stringstream ss;
    ss << '[';
    for (int i = 0; i < nvar; i++) {
        if (i > 0) ss << ',';
        ss << dep + i;
    }
    ss << "]:";

    if (atoms.size() + towers.size() > 1) ss << '(';

    bool is_first = true;
    for (auto atom : atoms) {
        if (is_first) is_first = false; else ss << '|';
        if (!atom.is_positive) ss << '~';
        ss << "E(" << atom.arg1 << ',' << atom.arg2 << ')';
    }
    for (auto tower : towers) {
        if (is_first) is_first = false; else ss << '|';
        ss << tower.to_string(dep + nvar);
    }

    if (atoms.size() + towers.size() > 1) ss << ')';

    return ss.str();
}

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
        if (ret != 0) return ret; // include utrefined
    }

    for (int i = 0; i < t1.towers.size(); i++) {
        int ret = Tower::cmp(t1.towers[i], t2.towers[i]);
        if (ret != 0) return ret; // include utrefined
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

Tower Tower::_to_tower(std::vector<Stack>::iterator lo,
                       std::vector<Stack>::iterator hi)
{
    assert(lo->is_open());
    std::vector<Atom> atoms;
    std::vector<Tower> towers;
    int n = 0, prev = 0;
    std::vector<Stack>::iterator last;
    for (auto it = lo + 1; it < hi; ++it) {
        if (it->is_atom()) {
            n++;
            if (prev == 0)
                atoms.emplace_back(it->is_positive(), it->arg1(), it->arg2());
        }
        else if (it->is_open()) {
            if (prev++ == 0) last = it;
        }
        else if (--prev == 0)
            towers.push_back(_to_tower(last, it));
    }
    assert(prev == 0);
    return Tower(n, lo->nq(), atoms, towers);
}

int _cmp(int l1, int r1, int l2, int r2, std::vector<Stack> &sk)
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
            assert(prev1 == 0);
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
            assert(prev2 == 0);
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

int _fitr(int k, std::vector<int> &p)
{
    return p[k] == k ? k : p[k] = _fitr(p[k], p);
}

void _union(int x, int y, std::vector<int> &p)
{
    assert(x < p.size() && y < p.size());
    p[_fitr(x, p)] = _fitr(y, p);
}

void _get_pmu(int k, int n, std::vector<int> &cnt, std::vector<std::vector<int>> &pmu)
{
    if (k == n) pmu.push_back(cnt);
    else {
        for (int i = 0; i < n; i++) {
            bool ok = true;
            for (auto &e : cnt) if (i == e) { ok = false; break; }
            if (ok) {
                cnt.push_back(i);
                _get_pmu(k + 1, n, cnt, pmu);
                cnt.pop_back();
            }
        }
    }
}

std::vector<Tower> Tower::search_raw(int n)
{
    std::vector<std::vector<std::vector<int>>> pmu;
    pmu.emplace_back();
    for (int i = 1; i <= n + 1; i++) {
        pmu.emplace_back();
        std::vector<int> cnt;
        _get_pmu(0, i, cnt, pmu.back());
    }

    std::vector<Tree> tr;
    std::vector<Stack> sk;
    std::vector<Tower> res;
    for (int i = 1; i <= n + 1; i++) {
        tr.emplace_back(true, i, 0);
        sk.emplace_back(i);
        dfs(n, 0, tr, sk, res, pmu);
        tr.pop_back();
        sk.pop_back();
    }
    return res;
}

void Tower::dfs(int n, int natoms, std::vector<Tree> &tr,
                std::vector<Stack> &sk, std::vector<Tower> &res,
                std::vector<std::vector<std::vector<int>>> &pmu)
{
    /*
     * is_unique itricates whether all permutation choosing of ancestors is unique
     * This is useful in pruning. If is_unique is true, then we can safely check
     * whether the current permutation of quantifiers satisfies the smallest setting
     * is_unique can be deduced from many sources. The most cheap atr convenient
     * is to deduce it from atomic children
     * If one node is unique then all its ancestors must be unique
     */
    /*
     * Calculate some basic information. The reason we don't
     * put this information in arguments is that calculation
     * is cheap atr we want to avoid complicated recursion
     * arguments
     * All closed braches should be ordered itself atr ordered
     * with its preceding brothers
     */

    // Check whether done
    if (tr.empty()) {
        res.push_back(_to_tower(sk.begin(), sk.end() - 1));
        return;
    }

    int dep = 0, nq = tr.back().nq;
    for (int i = 0; i < tr.size() - 1; i++) dep += tr[i].nq;

    /*
     * TODO: brach pruning here
     * Pruning one: the remaning slots in this branch are not enough for
     * non-satisfied ancestors
     */
    if (natoms < n) {
        // our greedy algorithm only works when there is a slot remaining
        // if the last element is close, we new a dummy open bracket with
        // nvar = 1
        bool is_new = sk.back().is_close();
        if (is_new) {
            sk.emplace_back(1);
            tr.emplace_back(false, 1, (int)sk.size());
        }

        // It's safe to set link[pos_left[k]] to dummy value
        for (auto &e : tr) sk[e.lpos].link = (int)sk.size();

        int _dep = 0;
        int need = 0;
        for (int k = 0; k < tr.size(); k++) {
            std::vector<int> p((unsigned)tr[k].nq);
            for (int i = 0; i < tr[k].nq; i++) p[i] = i;
            for (int i = tr[k].lpos + 1; i < sk.size(); ) {
                if (sk[i].is_atom()) {
                    if (sk[i].arg1() >= _dep && sk[i].arg2() >= _dep)
                        _union(sk[i].arg1() - _dep, sk[i].arg2() - _dep, p);
                    i++;
                }
                else {
                    for (int j = i + 1, last = -1; j < sk[i].link; j++)
                        if (sk[j].is_atom()) {
                            if (sk[j].arg1() >= _dep && sk[j].arg1() < _dep + tr[k].nq) {
                                if (last != -1) _union(last, sk[j].arg1() - _dep, p);
                                last = sk[j].arg1() - _dep;
                            }
                            if (sk[j].arg2() >= _dep && sk[j].arg2() < _dep + tr[k].nq) {
                                if (last != -1) _union(last, sk[j].arg2() - _dep, p);
                                last = sk[j].arg2() - _dep;
                            }
                        }
                    i = sk[i].link + 1;
                }
            }

            for (int i = 0; i < tr[k].nq; i++) if (_fitr(i, p) == i) need++;
            if (k + 1 < tr.size()) {
                // need to fitr whether current branch contains at least one quantifier,
                // if so then need += nparts - 1 else need += nparts
                for (int i = tr[k + 1].lpos + 1; i < sk.size(); i++)
                    if (sk[i].is_atom() && ((sk[i].arg1() >= _dep && sk[i].arg1() < _dep + tr[k].nq) ||
                                            (sk[i].arg2() >= _dep && sk[i].arg2() < _dep + tr[k].nq))) {
                        need--;
                        break;
                    }
            }
            else need--;

            _dep += tr[k].nq;
            if (need + natoms > n) break;
        }

        if (is_new) {
            sk.pop_back();
            tr.pop_back();
        }
        if (need + natoms > n) return;
    }



    /*
     * Case I: Open another bracket
     * Theoretically we can have maximum n - atoms + 1 new
     * quantifiers here. However we still need to leave a slot
     * for ancestor to come in, so this extra 1 should be removed
     * Check whether all atomic formulas in current subtree follows ordering
     * if unique is tree, and calculate unique value of its child
     */
    // only do unique checking the first time opening a bracket in current tree

    bool ok = true;
    if (sk.back().is_atom() || sk.back().is_open()) {
        // first time, need to calculate cuni
        if (sk.back().is_open()) tr.back().cuni = tr.back().uni && nq == 1;
        else if (nq == 1) tr.back().cuni = tr.back().uni;
        else if (!tr.back().uni) tr.back().cuni = false;
        else {
            // also need to check validity
            int count = 0;
            for (auto &p : pmu[nq]) {
                std::vector<Atom> atoms;
                for (int i = tr.back().lpos + 1; i < sk.size(); i++) {
                    assert(sk[i].is_atom());
                    atoms.emplace_back(sk[i].is_positive(),
                                       sk[i].arg1() < dep ? sk[i].arg1() : p[sk[i].arg1()-dep]+dep,
                                       sk[i].arg2() < dep ? sk[i].arg2() : p[sk[i].arg2()-dep]+dep);
                }
                std::sort(atoms.begin(), atoms.end());
                assert(!atoms.empty());
                int tmp = 0;
                for (int i = 0, j = tr.back().lpos + 1; i < atoms.size(); i++, j++) {
                    // we don't need to consider polarity
                    if (atoms[i].arg1 < sk[j].arg1() ||
                        (atoms[i].arg1 == sk[j].arg1() && atoms[i].arg2 < sk[j].arg2())) {
                        tmp = -1;
                        break;
                    }
                    if (atoms[i].arg1 > sk[j].arg1() ||
                        (atoms[i].arg1 == sk[j].arg1() && atoms[i].arg2 > sk[j].arg2())) {
                        tmp = 1;
                        break;
                    }
                }
                if (tmp == -1) { ok = false; break; }
                if (tmp == 0) count++;
            }
            if (ok) {
                assert(count >= 1);
                tr.back().cuni = (count == 1);
            }
        }
    }

    if (ok) {
        for (int k = 1; k <= n - natoms; k++) {
            tr.emplace_back(tr.back().cuni, k, sk.size());
            sk.emplace_back(k);
            dfs(n, natoms, tr, sk, res, pmu);
            tr.pop_back();
            sk.pop_back();
        }
    }


    /*
     * Case II: Close one bracket. Need to check ordering now.
     * Note: There is no need to check ordering inside. Only need
     * to check ordering with its preceding brothers
     * Check quantifers are well distributed in subtrees
     */
    if ((tr.size() > 1 || natoms == n) && !sk.back().is_open()) {
        auto k = (int)sk.size();
        sk.emplace_back();
        sk[tr.back().lpos].link = k;
        sk[k].link = tr.back().lpos;
        auto save = tr.back();
        tr.pop_back();


        if (sk[k].link == 0 || !sk[sk[k].link-1].is_close() ||
            _cmp(sk[sk[k].link-1].link, sk[k].link-1, sk[k].link, k, sk) < 0) {
            // check distribution of quantifiers
            // [dep, dep + nq)
            std::vector<int> p((unsigned)nq);
            for (int i = 0; i < nq; i++) p[i] = i;
            bool ok = true;

            if (!tr.empty()) {
                /* check whether contains direct parent
                 * [depp, dep), [dep, dep + nq)
                 */
                int depp = dep - tr.back().nq;
                ok = false;
                for (int i = sk[k].link + 1; i < k; i++)
                    if (sk[i].is_atom() && ((sk[i].arg1() >= depp && sk[i].arg1() < dep) ||
                                            (sk[i].arg2() >= depp && sk[i].arg2() < dep))) {
                        ok = true;
                        break;
                    }
            }

            if (ok) {
                for (int i = sk[k].link + 1; i < k; ) {
                    if (sk[i].is_atom()) {
                        // one of arguments must be in [dep, dep + nq) is
                        // guaranteed when it was created
                        if (sk[i].arg1() >= dep && sk[i].arg2() >= dep)
                            _union(sk[i].arg1() - dep, sk[i].arg2() - dep, p);
                        i++;
                    }
                    else {
                        for (int j = i + 1, last = -1; j < sk[i].link; j++)
                            if (sk[j].is_atom()) {
                                if (sk[j].arg1() >= dep && sk[j].arg1() < dep + nq) {
                                    if (last != -1) _union(last, sk[j].arg1() - dep, p);
                                    last = sk[j].arg1() - dep;
                                }
                                if (sk[j].arg2() >= dep && sk[j].arg2() < dep + nq) {
                                    if (last != -1) _union(last, sk[j].arg2() - dep, p);
                                    last = sk[j].arg2() - dep;
                                }
                            }
                        i = sk[i].link + 1;
                    }
                }
            }
            for (int i = 1; i < nq; i++) if (_fitr(0, p) != _fitr(i, p)) ok = false;
            if (ok) dfs(n, natoms, tr, sk, res, pmu);
        }

        sk.pop_back();
        tr.push_back(save);
    }

    /*
     * Case III: Create one new atomic formula
     */
    if (!sk.back().is_close() && natoms < n) {
        for (int i = 0; i < dep + nq; i++)
            for (int j = (i < dep ? dep : 0); j < dep + nq; j++)
                for (int k = 0; k < 2; k++) {
                    auto l = (int)sk.size();
                    sk.emplace_back((bool)k, i, j);
                    if (!sk[l-1].is_atom() || _cmp(l-1, l-1, l, l, sk) < 0) {
                        // check no complement
                        dfs(n, natoms + 1, tr, sk, res, pmu);
                    }
                    sk.pop_back();
                }
    }
}
