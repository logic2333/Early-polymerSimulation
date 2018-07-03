// inline

#pragma once
#include "chain.h"
#include <omp.h>
#include <sstream>
#include "gif.h"
#define R_FILTER 0x00FF0000
#define G_FILTER 0x0000FF00
#define B_FILTER 0x000000FF

const double quarter_sq2 = sq2 / 4;

struct Ccanvas {
	unsigned size; int R;
	POINT origin;
	vector<POINT> dir_2D;
	double ratio;
	IMAGE* img;
	GifWriter w;
	BYTE* img_data;
	unsigned t;
	POINT from3D(const coord_3D& _c, bool relative = true) noexcept {
		POINT res;
		res.x = ratio * (_c.y - _c.x * quarter_sq2);
		res.y = ratio * (_c.x * quarter_sq2 - _c.z);
		if (!relative) {
			res.x += origin.x; res.y += origin.y;
		}
		return res;
	}
	Ccanvas() noexcept {}
	Ccanvas(unsigned _size) {
		size = (_size > 199 ? _size : 0);
		if (size) {
			R = (size - 100) * sq2 / (sq2 + 0.5);
			origin.x = 50 + R * quarter_sq2;
			origin.y = 50 + R;
			ratio = double(R) / (coord_3D::_r - 1);
			img = new IMAGE(size, size);
			dir_2D.resize(20);
			for (unsigned i = 0; i < 20; i++)
				dir_2D[i] = from3D(dir_coord[i]);
		}
	}
	void init() {
		SetWorkingImage(img);
		setbkcolor(WHITE);
		cleardevice();
		moveto(origin.x, origin.y);
		setlinestyle(PS_DASH);
		setlinecolor(DARKGRAY);
		linerel(-R * quarter_sq2, R * quarter_sq2);
		linerel(R, 0);
		linerel(R * quarter_sq2, -R * quarter_sq2);
		linerel(-R, 0);
		linerel(0, -R);
		linerel(-R * quarter_sq2, R * quarter_sq2);
		linerel(0, R);
		moveto(origin.x, 50);
		linerel(R, 0);
		linerel(0, R);
		moveto(origin.x - R * quarter_sq2, 50 + R * quarter_sq2);
		linerel(R, 0);
		linerel(0, R);
		line(50 + R, 50 + R * quarter_sq2, origin.x + R, 50);
	}
	void drawBlock(const Cpoint& p, bool is_erase) {
		const POINT coord_2D = from3D(p, false);
#pragma omp critical
		{
			COLORREF color;
			if (is_erase) color = WHITE;
			else switch (p.occup) {
			case 'A': color = RED; break;
			case 'B': color = GREEN; break;
			case 'C': color = BLUE; break;
			case 'D': color = WHITE; break;
			}
			setfillcolor(color);
			solidcircle(coord_2D.x, coord_2D.y, 3);
			moveto(coord_2D.x, coord_2D.y);
		}
	}
	void drawBond(const Cblock& block, bool is_forward) {
		setlinecolor(BLACK);
		setlinestyle(PS_SOLID | PS_ENDCAP_FLAT, 2);
		if (is_forward) {
			if (block.next_dir != st) {
				const POINT dir = dir_2D[block.next_dir];
				const POINT now_coord = from3D(block, false);
				moveto(now_coord.x, now_coord.y);
				if (block.out_border) {
					linerel(dir.x / 2, dir.y / 2);
					const POINT next_posit = from3D(find_adjacent(block, false), false);
					moveto(next_posit.x, next_posit.y);
					linerel(-dir.x / 2, -dir.y / 2);
				}
				else linerel(dir.x, dir.y);
			}
		}
		else {
			drawBond(find_adjacent(block, true), true);
		}
	}
	void drawChains(const bitset<5>& ipt) {
#pragma omp parallel for
		for (int i = 0; i < chains_act::chains.size(); i++) {
			switch (chains_act::chains[i].occup) {
				case 'A': if (ipt[0]) drawBlock(chains_act::chains[i], false); break;
				case 'B': if (ipt[1]) drawBlock(chains_act::chains[i], false); break;
				default: break;
			}
			if (ipt[4])
#pragma omp critical
			{
				drawBond(chains_act::chains[i], true);
			}		
		}
	}
	void drawSol(const bitset<5>& ipt) {
		if (ipt[2])
#pragma omp parallel for
			for (int i = 0; i < chains_act::all.size(); i++)
				if (chains_act::all[i].occup == 'C')
					drawBlock(chains_act::all[i], false);		
	}
	void drawPhase_init() {
		SetWorkingImage(img);
		cleardevice();
	}
	void shotImage(unsigned now_MCS, const string& path) {
		ostringstream os;
		os << path << now_MCS << ".bmp";
		saveimage(os.str().c_str(), img);
	}
	void writeGIF() {
		const auto data = GetImageBuffer(img);
#pragma omp parallel for
		for (int i = 0; i < size * size; i++) {
			img_data[4 * i] = (R_FILTER & data[i]) >> 16;
			img_data[4 * i + 1] = (G_FILTER & data[i]) >> 8;
			img_data[4 * i + 2] = B_FILTER & data[i];
			img_data[4 * i + 3] = 0x00;
		}
		GifWriteFrame(&w, img_data, size, size, t);
	}
	void initGIF(const char* gif_name, unsigned fps, unsigned stop) {
		t = 100 / fps;
		GifBegin(&w, gif_name, size, size, t);
		img_data = new BYTE[size * size * 4];
		for (unsigned i = 0; i < stop; i++)
			writeGIF();
	}
	void endGIF() {
		delete[] img_data;
		delete img;
		GifEnd(&w);
	}
	void openScreen() {
		initgraph(size, size, SHOWCONSOLE);
	}
	void putScreen() {
		SetWorkingImage();
		putimage(0, 0, img);
	}
	void closeScreen() {
		closegraph();
	}
};