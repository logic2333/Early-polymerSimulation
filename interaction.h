#pragma once
#include "chain.h"
#include <interpolation.h>
#include <algorithm>
using namespace alglib;


struct Interaction {
	real_2d_array data;
	real_1d_array X;
	spline1dinterpolant* funcs;
	unsigned funcsCnt, colCnt;
	vector<double> res;
	bool is_continuous;
	static map<CInteract, double> m;
	Interaction(const string& fileName, bool c) : is_continuous(c) {
		read_csv(fileName.c_str(), ',', 0, data);
		funcsCnt = data.rows() - 1;
		colCnt = data.cols() - 1;
		if (is_continuous)
			funcs = new spline1dinterpolant[funcsCnt];
		res.resize(funcsCnt);
		X.attach_to_ptr(colCnt, &data[0][1]);
		if (is_continuous) {
			const ae_int_t natural = 2;
			for (unsigned i = 0; i < funcsCnt; i++) {
				real_1d_array Y;
				Y.attach_to_ptr(colCnt, &data[i + 1][1]);
				spline1dbuildcubic(X, Y, colCnt, natural, 0.0, natural, 0.0, funcs[i]);
			}
		}
	}
	~Interaction() {
		if (is_continuous)
			delete[] funcs;
	}
	void update_interactions(unsigned x) {
		double* Xbegin = X.getcontent();
		for (unsigned i = 0; i < funcsCnt; i++)
			if (is_continuous)
				res[i] = spline1dcalc(funcs[i], x);
			else {
				unsigned t = upper_bound(Xbegin, Xbegin + colCnt, x) - Xbegin;
				res[i] = data[i + 1][t];
			}
		for (unsigned i = 0; i < funcsCnt; i++)
			switch (i) {
				case 0: m[make_pair('A', 'A')] = res[i]; break;
				case 1: m[make_pair('A', 'B')] = m[make_pair('B', 'A')] = res[i]; break;
				case 2: m[make_pair('A', 'C')] = m[make_pair('C', 'A')] = res[i]; break;
				case 3: m[make_pair('A', 'D')] = m[make_pair('D', 'A')] = res[i]; break;
				case 4: m[make_pair('B', 'B')] = res[i]; break;
				case 5: m[make_pair('B', 'C')] = m[make_pair('C', 'B')] = res[i]; break;
				case 6: m[make_pair('B', 'D')] = m[make_pair('D', 'B')] = res[i]; break;
				case 7: m[make_pair('C', 'C')] = res[i]; break;
				case 8: m[make_pair('C', 'D')] = m[make_pair('D', 'C')] = res[i]; break;
				case 9: m[make_pair('D', 'D')] = res[i]; break;
				default: break;
			}
	}
};

