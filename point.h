// inline

#pragma once
#include <cmath>
#include <vector>
#define OFFSET 0.0001
#define MAXSPACE 1000000
using namespace std;

const double sq2 = sqrt(2);

typedef enum { st, f, r, b, l, u, d, ul, ur, uf, ub, dl, dr, df, db, lf, lb, rf, rb, farr } DIRECT;

const vector<double> cosine_angles = { sq2 / 2, 1.0 / 2, 0.0 };

struct coord_3D {
	static int _r; // 坐标上限，绝对不能是unsigned！
	int x, y, z;		
	coord_3D() noexcept {}
	coord_3D(int _x, int _y, int _z) noexcept {
		x = _x; y = _y; z = _z;
	}
	coord_3D(unsigned hash_val) {
		x = hash_val % _r;
		hash_val /= _r;
		y = hash_val % _r;
		hash_val /= _r;
		z = hash_val;
	}
	unsigned hash() const noexcept {
		return z * _r * _r + y * _r + x;
	}
	unsigned len() const noexcept {
		return x * x + y * y + z * z;
	}
	bool operator<(const coord_3D& another) const noexcept {
		if (z < another.z) return true;
		else if (z == another.z)
			if (y < another.y) return true;
			else if (y == another.y) return x < another.x;
		return false;
	}
	coord_3D operator-(const coord_3D& another) const noexcept {
		coord_3D res;
		res.x = x - another.x;
		res.y = y - another.y;
		res.z = z - another.z;
		return res;
	}
	bool operator==(const coord_3D& another) const noexcept {
		if (x == another.x && y == another.y && z == another.z) return true;
		else return false;
	}
	unsigned operator*(const coord_3D& another) const noexcept {
		return x * another.x + y * another.y + z * another.z;
	}
	// 输出值表示是否超出边界
	bool _move(DIRECT _d) noexcept {
		switch (_d) {
			case u: z++; break; 
			case d: z--; break;
			case l: y--; break;
			case r: y++; break;
			case f: x++; break;
			case b: x--; break;
			case ul: z++; y--; break;
			case ur: z++; y++; break;
			case uf: z++; x++; break;
			case ub: z++; x--; break;
			case dl: z--; y--; break;
			case dr: z--; y++; break;
			case df: z--; x++; break;
			case db: z--; x--; break;
			case lf: y--; x++; break;
			case lb: y--; x--; break;
			case rf: y++; x++; break;
			case rb: y++; x--; break;
			default: break;
		}
		bool res = false;
		if (x >= _r) { x -= _r; res = true; }
		else if (x < 0) { x += _r; res = true; }
		if (y >= _r) { y -= _r; res = true; }
		else if (y < 0) { y += _r; res = true; }
		if (z >= _r) { z -= _r; res = true; }
		else if (z < 0) { z += _r; res = true; }
		return res;
	}	
	void real_move(DIRECT _d) noexcept {
		switch (_d) {
			case u: z++; break;
			case d: z--; break;
			case l: y--; break;
			case r: y++; break;
			case f: x++; break;
			case b: x--; break;
			case ul: z++; y--; break;
			case ur: z++; y++; break;
			case uf: z++; x++; break;
			case ub: z++; x--; break;
			case dl: z--; y--; break;
			case dr: z--; y++; break;
			case df: z--; x++; break;
			case db: z--; x--; break;
			case lf: y--; x++; break;
			case lb: y--; x--; break;
			case rf: y++; x++; break;
			case rb: y++; x--; break;
			default: break;
		}
	}
};

// 定义DIRECT到coord_3D的转换
const vector<coord_3D> dir_coord = {
	coord_3D(0, 0, 0),		// st
	coord_3D(1, 0, 0),		// f
	coord_3D(0, 1, 0),		// r
	coord_3D(-1, 0, 0),		// b
	coord_3D(0, -1, 0),		// l
	coord_3D(0, 0, 1),		// u
	coord_3D(0, 0, -1),		// d
	coord_3D(0, -1, 1),		// ul
	coord_3D(0, 1, 1),		// ur
	coord_3D(1, 0, 1),		// uf
	coord_3D(-1, 0, 1),		// ub
	coord_3D(0, -1, -1),	// dl
	coord_3D(0, 1, -1),		// dr
	coord_3D(1, 0, -1),		// df
	coord_3D(-1, 0, -1),	// db
	coord_3D(1, -1, 0),		// lf
	coord_3D(-1, -1, 0),	// lb
	coord_3D(1, 1, 0),		// rf
	coord_3D(-1, 1, 0),		// rb
	coord_3D(2, 2, 2)		// farr
};

// 返回一个从_A到_B的方向向量，若another与this不相邻，返回(2, 2, 2)
coord_3D adjacent(const coord_3D& _A, const coord_3D& _B);

unsigned angle(const coord_3D& _A, const coord_3D& _B);