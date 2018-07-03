#include "stdafx.h"
#include "chain.h"
#include <algorithm>
#include "interaction.h"

int coord_3D::_r;
unsigned Cblock::rigidity;
vector<coord_3D> chains_act::chain_heads;
vector<coord_3D> chains_act::chain_tails;
vector<Cblock> chains_act::chains;
vector<Cpoint> chains_act::all;
vector<unsigned> chains_act::molCnt;

map<CInteract, double> Interaction::m;


const map<DIRECT, DIRECT> opposite = {
	{ f, b },{ b, f },{ l, r },{ r, l },{ u, d },{ d, u },{ farr, farr },{ st, st },
	{ ul, dr },{ dr, ul },{ ur, dl },{ dl, ur },{ uf, db },{ db, uf },
	{ lf, rb },{ rb, lf },{ lb, rf },{ rf, lb },{ df, ub },{ ub, df }
};

coord_3D adjacent(const coord_3D& _A, const coord_3D& _B) {
	coord_3D tmp;
	for (unsigned i = 1; i < 19; i++) {
		tmp = _A;
		tmp._move(DIRECT(i));
		if (tmp == _B)
			return dir_coord[i];
	}
	return coord_3D(2, 2, 2);
}

unsigned angle(const coord_3D& _A, const coord_3D& _B) {
	if (_A == coord_3D(0, 0, 0) || _B == coord_3D(0, 0, 0)) return 3;
	const double cosine_angle = _A * _B / (sqrt(_A.len()) * sqrt(_B.len()));
	for (unsigned i = 0; i < cosine_angles.size(); i++)
		if (abs(cosine_angle - cosine_angles[i]) < OFFSET) return i;
	return 3;
}

void p2P(vector<double>& _p) {
	double sum = _p[0];
	for (unsigned i = 1; i < _p.size(); i++) {
		_p[i] += sum;
		sum = _p[i];
	}
	for (unsigned i = 0; i < _p.size(); i++)
		_p[i] /= sum;
}

coord_3D Cpoint::movable(DIRECT d) const {
	coord_3D dest_coord = *this; dest_coord._move(d);
	const char dest_occup = chains_act::all[dest_coord.hash()].occup;
	if (dest_occup < 'C' || dest_occup == occup) return coord_3D(-1, -1, -1);
	//溶剂与链段的交换，仅从链段方进行；不与自身相同的分子交换
	return dest_coord;
}

double calcEnergy(char chain, const vector<double>& surroundings) {
	double res = 0.0;
	pair<char, char> p; p.first = chain;
	for (char c = 'A'; c < 'E'; c++) {
		p.second = c;
		res += surroundings[c - 'A'] * Interaction::m.at(p);
	}	
	return res;
}

double Cpoint::getDeltaEnergy(DIRECT d) const {
	coord_3D dest_coord = *this; dest_coord._move(d);
	const Cpoint& sol = chains_act::all[dest_coord.hash()];
	vector<double> me_surroundings, sol_surroundings;
	getSurroundings(d, me_surroundings);
	sol.getSurroundings(opposite.at(d), sol_surroundings);
	double sol_energy_delta, me_energy_delta;
	sol_energy_delta = calcEnergy(sol.occup, sol_surroundings) - calcEnergy(sol.occup, me_surroundings);
	me_energy_delta = calcEnergy(occup, me_surroundings) - calcEnergy(occup, sol_surroundings);
	return sol_energy_delta + me_energy_delta;
}

vector<double> Cpoint::getProb() {
	static vector<double> p; 
	p.clear(); p.resize(19);
	p[0] = 1.0;
	coord_3D dest_coord;
	for (unsigned i = 1; i < 19; i++) {
		dest_coord = movable(DIRECT(i));
		if (dest_coord.x != -1)
			p[i] = exp(getDeltaEnergy(DIRECT(i)));
	}
	p2P(p);
	return p;
}

Cblock::Cblock(int _x, int _y, int _z, int _c, char _o) noexcept : Cpoint(_x, _y, _z, _o) {
	out_border = false;
	before_dir = next_dir = st;	// 这里表示原位，首尾链段需要
	belong_chain = _c;
}

Cblock::Cblock(const coord_3D& c, int _c, char _o) noexcept : Cpoint(c, _o) {
	out_border = false;
	before_dir = next_dir = st; // 这里表示原位，首尾链段需要
	belong_chain = _c;
}

vector<Cblock>::iterator fast_find(vector<Cblock>& con, const coord_3D& val) {
	auto res = lower_bound(con.begin(), con.end(), val);
	if (res == con.end()) return res;
	if (*res == val) return res;
	return con.end();
}

Cblock Cblock::grow(int _c, char _o) {
	coord_3D new_moved;
	Cblock new_block;
	for (unsigned i = 1; i < 19; i++)
		if (angle(dir_coord[before_dir], dir_coord[i]) >= Cblock::rigidity) {
			new_moved = *this;
			out_border = new_moved._move(DIRECT(i));
			if (chains_act::all[new_moved.hash()].occup == 'D') {
				next_dir = DIRECT(i);
				new_block = Cblock(new_moved, _c, _o);
				new_block.before_dir = opposite.at(next_dir);
				return new_block;
			}
		}
	// exception
}

bool exist_bond_another_angles(const coord_3D& _first, const coord_3D& _second) {
	static pair<coord_3D, coord_3D> res;
	if ((_first - _second).len() == 2) {
		if (_first.x == _second.x) {
			res.first.x = res.second.x = _first.x;
			res.first.y = _second.y;
			res.first.z = _first.z;
			res.second.y = _first.y;
			res.second.z = _second.z;
		}
		else if (_first.y == _second.y) {
			res.first.y = res.second.y = _first.y;
			res.first.x = _second.x;
			res.first.z = _first.z;
			res.second.x = _first.x;
			res.second.z = _second.z;
		}
		else {
			res.first.z = res.second.z = _first.z;
			res.first.x = _second.x;
			res.first.y = _first.y;
			res.second.x = _first.x;
			res.second.y = _second.y;
		}
	}
	else return false;
	auto now = fast_find(chains_act::chains, res.first);
	if (now == chains_act::chains.end()) return false;
	coord_3D bef = *now; bef._move(now->before_dir);
	if (bef == res.second) return true;
	coord_3D nex = *now; nex._move(now->next_dir);
	return nex == res.second;
}

bool Cblock::cross(const coord_3D& _bef, const coord_3D& _nex, const coord_3D& _dest) const {
	if (exist_bond_another_angles(*this, _dest)) return true;	// 穿越
	if (exist_bond_another_angles(_dest, _bef)) return true;	// 交叉
	if (exist_bond_another_angles(_dest, _nex)) return true;	// 交叉
	return false;
}

Cblock& find_adjacent(Cblock i, bool is_before) {
	i._move(is_before ? i.before_dir : i.next_dir);
	return *fast_find(chains_act::chains, i);
}

coord_3D Cblock::movable(DIRECT d) const {
	const Cblock& before = find_adjacent(*this, true);
	const Cblock& next = find_adjacent(*this, false);
	coord_3D dest_coord = *this; dest_coord._move(d);
	if (chains_act::all[dest_coord.hash()].occup < 'C') return coord_3D(-1, -1, -1);	// 不允许链段之间的交换
	coord_3D new_before_dir, new_next_dir;
	new_before_dir = adjacent(dest_coord, before);
	new_next_dir = adjacent(dest_coord, next);
	if (new_before_dir.x == 2) return coord_3D(-1, -1, -1);		// 键长条件
	if (new_next_dir.x == 2) return coord_3D(-1, -1, -1);		// 键长条件
	if (angle(new_before_dir, new_next_dir) < rigidity) return coord_3D(-1, -1, -1);	// 键角（链刚度）条件
	if (angle(new_before_dir, dir_coord[opposite.at(before.before_dir)]) < rigidity) return coord_3D(-1, -1, -1);
	if (angle(new_next_dir, dir_coord[opposite.at(next.next_dir)]) < rigidity) return coord_3D(-1, -1, -1);
	if (cross(before, next, dest_coord)) return coord_3D(-1, -1, -1); // 不允许交叉与穿越
	return dest_coord;
}

unsigned short Cpoint::getSurroundings(DIRECT _d, bool) const {
	vector<unsigned> res; res.resize(3);
	for (unsigned i = 1; i < 7; i++)
		if (i != _d) {
			coord_3D dest_coord = *this; dest_coord._move(DIRECT(i));
			const char dest_occup = chains_act::all[dest_coord.hash()].occup;
			if (dest_occup < 'D') res[dest_occup - 'A']++;
		}
	return res[0] * 64 + res[1] * 8 + res[2];
}

void Cpoint::getSurroundings(DIRECT _d, vector<double>& res) const {
	res.clear(); res.resize(4);
	coord_3D dest_coord;
	for (unsigned i = 1; i < 19; i++)
		if (i != _d) {
			dest_coord = *this; dest_coord._move(DIRECT(i));
			const char dest_occup = chains_act::all[dest_coord.hash()].occup;
			if (i > 6) res[dest_occup - 'A'] += 1 / 8.0;
			else res[dest_occup - 'A']++;
		}
}

double Cblock::getDeltaEnergy(DIRECT d) const {
	coord_3D dest_coord = *this; dest_coord._move(d);
	const Cpoint& sol = chains_act::all[dest_coord.hash()];
	vector<double> chain_surroundings, sol_surroundings;
	getSurroundings(d, chain_surroundings);
	sol.getSurroundings(opposite.at(d), sol_surroundings);
	const double sol_energy_delta = calcEnergy(sol.occup, sol_surroundings) - calcEnergy(sol.occup, chain_surroundings);
	const Cblock& before = find_adjacent(*this, true);
	const Cblock& next = find_adjacent(*this, false);
	if (before_dir != st)
		if ((before - sol).len() == 1) sol_surroundings[before.occup - 'A']--;
		else sol_surroundings[before.occup - 'A'] -= 1 / 8.0;
	if (next_dir != st)
		if ((next - sol).len() == 1) sol_surroundings[next.occup - 'A']--;
		else sol_surroundings[next.occup - 'A'] -= 1 / 8.0;
	switch ((before - *this).len()) {
		case 1: chain_surroundings[before.occup - 'A']--; break;
		case 2: chain_surroundings[before.occup - 'A'] -= 1 / 8.0; break;
		default: break;
	}
	switch ((next - *this).len()) {
		case 1: chain_surroundings[next.occup - 'A']--; break;
		case 2: chain_surroundings[next.occup - 'A'] -= 1 / 8.0; break;
		default: break;
	}
	const double chain_energy_delta = calcEnergy(occup, chain_surroundings) - calcEnergy(occup, sol_surroundings);
	return chain_energy_delta + sol_energy_delta;
}

vector<double> Cblock::getProb() {
	static vector<double> p; 
	p.clear(); p.resize(19);
	p[0] = 1.0;
	coord_3D dest_coord;
	for (unsigned i = 1; i < 19; i++)
		if (i != next_dir && i != before_dir /*不重叠条件*/) {
			dest_coord = movable(DIRECT(i));
			if (dest_coord.x != -1)
				p[i] = exp(getDeltaEnergy(DIRECT(i)));
		}
	p2P(p);
	return p;
}

void Cpoint::move(DIRECT d) {
	coord_3D dest_coord = *this; dest_coord._move(d);
	swap(occup, chains_act::all[dest_coord.hash()].occup);
}

void my_insertion_sort(const Cblock b, bool smaller, const vector<Cblock>::iterator& p) {
	vector<Cblock>::iterator insert_posit;
	if (smaller) {
		insert_posit = lower_bound(chains_act::chains.begin(), p, b);
		copy_backward(insert_posit, p, p + 1);
		*insert_posit = b;
	}
	else {
		insert_posit = lower_bound(p + 1, chains_act::chains.end(), b);
		copy(p + 1, insert_posit, p);
		*(insert_posit - 1) = b;
	}
}

void Cblock::move(DIRECT d, vector<Cblock>::iterator posit) {
	coord_3D dest_coord = *this; dest_coord._move(d);
	const bool move_smaller = dest_coord < *this;
	chains_act::all[hash()].occup = chains_act::all[dest_coord.hash()].occup;
	chains_act::all[dest_coord.hash()].occup = occup;
	if (before_dir != st)
		before_dir = DIRECT(find(dir_coord.begin(), dir_coord.end(), (dir_coord[before_dir] - dir_coord[d])) - dir_coord.begin());
	else
		chains_act::chain_heads[belong_chain].real_move(d);
	if (next_dir != st)
		next_dir = DIRECT(find(dir_coord.begin(), dir_coord.end(), (dir_coord[next_dir] - dir_coord[d])) - dir_coord.begin());
	else
		chains_act::chain_tails[belong_chain].real_move(d);
	x = dest_coord.x; y = dest_coord.y; z = dest_coord.z;
	coord_3D temp1 = *this;
	out_border = temp1._move(next_dir);
	const Cblock temp2 = *this;
	my_insertion_sort(*this, move_smaller, posit);
	if (temp2.before_dir != st) {
		Cblock& before = find_adjacent(temp2, true);
		before.next_dir = opposite.at(temp2.before_dir);
		coord_3D temp3 = before;
		before.out_border = temp3._move(before.next_dir);
	}
	if (temp2.next_dir != st)
		find_adjacent(temp2, false).before_dir = opposite.at(temp2.next_dir);
}

void solInit(unsigned C_cnt, unsigned sol_cnt, unsigned solDis) {
	unsigned cnt = 0;
	switch (solDis) {
	case 0: {
		for (int i = coord_3D::_r - 1; i >= 0; i--)
			for (int j = coord_3D::_r - 1; j >= 0; j--)
				for (int k = coord_3D::_r - 1; k >= 0; k--)
					if (cnt++ < C_cnt)
						chains_act::all[coord_3D(k, j, i).hash()].occup = 'C';
					else return;
					break;
	}
	case 1: {
		for (unsigned i = 0; i < coord_3D::_r; i++)
			for (unsigned j = 0; j < coord_3D::_r; j++)
				for (unsigned k = 0; k < coord_3D::_r; k++)
					if (cnt < C_cnt) {
						auto& it = chains_act::all[coord_3D(k, j, i).hash()];
						if (it.occup == 'D') {
							it.occup = 'C'; cnt++;
						}
					}
					else return;
					break;
	}
	case 2: {
		const double thresh = double(C_cnt) / sol_cnt;
		for (int i = coord_3D::_r - 1; i >= 0; i--)
			for (int j = coord_3D::_r - 1; j >= 0; j--)
				for (int k = coord_3D::_r - 1; k >= 0; k--) {
					auto& it = chains_act::all[coord_3D(k, j, i).hash()];
					if (it.occup == 'D') {
						if (randomreal() < thresh) {
							it.occup = 'C'; cnt++;
						}
						if (cnt >= C_cnt) return;
					}
				}
		break;
	}
	}
}

void chains_act::init(const Input& _i) {
	Cblock::rigidity = _i.chain_rigid;
	unsigned As = 0, Bs = 0;
	for (auto& c : _i.chain_distr) {
		if (c == 'A') As++;
		if (c == 'B') Bs++;
	}
	coord_3D::_r = pow(double(_i.blockCnt) / _i.blockConcent * 100, 1.0 / 3) + 0.5;
	const unsigned chainCnt = _i.blockCnt / _i.chain_len;
	all.reserve(coord_3D::_r * coord_3D::_r * coord_3D::_r);
	molCnt.resize(4);
	molCnt[0] = _i.blockCnt * As / (As + Bs);
	molCnt[1] = _i.blockCnt - molCnt[0];
	molCnt[2] = double(all.capacity()) * _i.C_Concent / 100;
	molCnt[3] = all.capacity() - _i.blockCnt - molCnt[2];
	for (unsigned i = 0; i < coord_3D::_r; i++)
		for (unsigned j = 0; j < coord_3D::_r; j++)
			for (unsigned k = 0; k < coord_3D::_r; k++)
				all.emplace_back(k, j, i, 'D');
	Cblock startblock(0, 0, 0, 0, _i.chain_distr[0]);
	chain_heads.reserve(chainCnt);
	chain_tails.reserve(chainCnt);
	chains.reserve(_i.blockCnt);
	Cblock last_block;
	for (unsigned i = 0; i < chainCnt; i++) {
		chain_heads.push_back(startblock);
		chains.push_back(startblock);
		all[startblock.hash()].occup = startblock.occup;
		unsigned k = 1; if (k == _i.chain_distr.length()) k = 0;
		for (unsigned j = 1; j < _i.chain_len; j++) {
			chains.push_back(chains.back().grow(i, _i.chain_distr[k]));
			all[chains.back().hash()].occup = chains.back().occup;
			k++; if (k == _i.chain_distr.length()) k = 0;
			if (j == _i.chain_len - 1) {
				last_block = chains.back();
				chain_tails.push_back(last_block);
				startblock = last_block.grow(i + 1, _i.chain_distr[k]);
				startblock.before_dir = st;
			}
		}
	}
	if (_i.C_Concent != 0)
		solInit(molCnt[2], molCnt[2] + molCnt[3], _i.C_distr);
}

vector<unsigned> chains_act::move() {		// 1 MCS
	static vector<unsigned> res;
	res.clear(); res.resize(3);
	for (auto i = all.rbegin(); i != all.rend(); i++) {
		if (i->occup < 'C') {
			auto j = fast_find(chains, *i);
			const vector<double>& prob = j->getProb();
			DIRECT d = DIRECT(lower_bound(prob.begin(), prob.end(), randomreal()) - prob.begin());
			if (d != st) {
				res[i->occup - 'A']++; j->move(d, j);
			}
		}
		else if (i->occup == 'C') {
			const vector<double>& prob = i->getProb();
			DIRECT d = DIRECT(lower_bound(prob.begin(), prob.end(), randomreal()) - prob.begin());
			if (d != st) {
				res[2]++; i->move(d);
			}
		}
	}
	return res;
}