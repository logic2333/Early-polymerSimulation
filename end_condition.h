// inline

#pragma once
#include "drawing.h"
#include "interaction.h"

struct endCondition {
	unsigned MCSMax;
	unsigned which_one;		 // 0-A, 1-B, 2-C, 3-A+B, 4-A+B+C, 5-don't use
	unsigned mol_percent;
	unsigned completion;
	unsigned no_move_MCS;
	unsigned gif_shot_rate;
	endCondition() noexcept {}
	// MCSMax, which_one, mol_percent, no_move_MCS
	endCondition(unsigned _MCSMax, unsigned _w, unsigned _m, unsigned _n) noexcept {
		MCSMax = _MCSMax; which_one = _w; mol_percent = _m; no_move_MCS = _n;
		completion = 0;
		gif_shot_rate = (_MCSMax < 10000? 1 : _MCSMax / 10000);
	}
	bool check_reached(unsigned now_MCS, const vector<double>& movage) {
		if (now_MCS >= MCSMax) return true;
		if (which_one != 5)
			if (unsigned(movage[which_one] * 100) < mol_percent) completion++;
			else completion = 0;
		if (completion > no_move_MCS) return true;
		return false;
	}
};

bool one_mcs(Ccanvas& _canvas, endCondition& EC, Interaction& it, unsigned& _MCSCnt) {
	if (_canvas.size) {
		if (Input::shotRate && (_MCSCnt % Input::shotRate == 0)) _canvas.shotImage(_MCSCnt, Input::working_path);
		_canvas.init();
	}
	it.update_interactions(_MCSCnt);
	if (_MCSCnt == Input::add_time_stamp) {
		int max = Input::add_amount * chains_act::all.size() / 100;
#pragma omp parallel for
		for (int i = 0; i < max; i++) {
			unsigned t = randominteger(chains_act::all.size());
			while (chains_act::all[t].occup != 'D')
				t = randominteger(chains_act::all.size());
#pragma omp critical
			{
				chains_act::all[t].occup = 'C';
			}
		}
	}
	Cmonitor::movage = chains_act::move();
	if (_MCSCnt % EC.gif_shot_rate == 0) Cmonitor::getDistance();
	if (_canvas.size) {
		_canvas.drawChains(Input::display_elements);
		_canvas.drawSol(Input::display_elements);	
		if (Input::gif && (_MCSCnt % EC.gif_shot_rate == 0)) _canvas.writeGIF();
		if (Input::screen) _canvas.putScreen();
	}
	return EC.check_reached(++_MCSCnt, Cmonitor::getMovage());
}