#pragma once
#include <map>
#include "point.h"
#include "input.h"


typedef pair<char, char> CInteract;

struct Cpoint : coord_3D {
	char occup;
	Cpoint() noexcept {}
	Cpoint(int _x, int _y, int _z, char _o) noexcept : coord_3D(_x, _y, _z), occup(_o) {}
	Cpoint(const coord_3D& c, char _o) noexcept : coord_3D(c), occup(_o) {}
	Cpoint(int _hash_val, char _o) noexcept : coord_3D(_hash_val), occup(_o) {}
	virtual ~Cpoint() noexcept {}
	virtual vector<double> getProb();
	virtual coord_3D movable(DIRECT d) const;
	virtual double getDeltaEnergy(DIRECT d) const;
	void getSurroundings(DIRECT _d, vector<double>& res) const;
	unsigned short getSurroundings(DIRECT _d, bool) const;
	virtual void move(DIRECT d);
};

struct Cblock : Cpoint {
	static unsigned rigidity;
	bool out_border;
	DIRECT next_dir, before_dir;
	int belong_chain;
	Cblock() noexcept {}
	Cblock(int _x, int _y, int _z, int _c, char _o) noexcept;
	Cblock(const coord_3D& c, int _c, char _o) noexcept;
	Cblock grow(int _c, char _o);
	vector<double> getProb() override;
	coord_3D movable(DIRECT d) const override;
	double getDeltaEnergy(DIRECT d) const override;
	bool cross(const coord_3D& _bef, const coord_3D& _nex, const coord_3D& _dest) const;
	void move(DIRECT d, vector<Cblock>::iterator posit);
};

vector<Cblock>::iterator fast_find(vector<Cblock>& con, const coord_3D& val);

struct chains_act {
	static vector<unsigned> molCnt;
	static vector<coord_3D> chain_heads;
	static vector<coord_3D> chain_tails;
	static vector<Cblock> chains;
	static vector<Cpoint> all;
	static map<CInteract, double> interactions;
	/* n聚合物链段的量，c1聚合物浓度，c2第二溶剂浓度，len平均聚合度
	distr 链段分布，rigid链刚度，溶剂用C、D表示，solDis溶剂初始分布情况 */
	static void init(const Input& _i);
	static vector<unsigned> move();
};

Cblock& find_adjacent(Cblock i, bool is_before);
double calcEnergy(char chain, const vector<double>& surroundings);