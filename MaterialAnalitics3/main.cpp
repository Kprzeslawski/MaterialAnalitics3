#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>

#define DATA double

void euler_2(DATA y_at_beg, DATA beg, DATA end, int steps, DATA yprim(DATA x, DATA y, DATA y_t_tcr,
	DATA a1, DATA a2, DATA a3, DATA a8, DATA p_cr, DATA t_cr
), DATA* a);
DATA fun_ivm(DATA x, DATA y, DATA y_t_tcr,
	DATA a1, DATA a2, DATA a3, DATA a8, DATA p_cr, DATA t_cr
);

DATA e_dot = 1;
DATA t_maks = 1;
DATA tK = 575. + 273.15;

//stale 

DATA R = 8.314462;

//miedz 

DATA b = 0.25e-9;
DATA D = 30;
DATA u = 45000;
DATA Q = 238000;
DATA p0 = 1e4;

int main() {
	DATA a[] = {2.1e-3, 176., 19.5e3, 0.000148 * 3e10, 151.0e3, 0.973, 5.77, 1.0, 0.0, 0.262, 0.0e13, 0.000605e13, 0.167};
	euler_2(p0, 0, t_maks, 1000, fun_ivm, a);
}

DATA fun_ivm(DATA x, DATA y, DATA y_t_tcr, 
	DATA a1, DATA a2, DATA a3, DATA a8, DATA p_cr, DATA t_cr
	) 
	{
	if (x >= t_cr)
		return a1 * e_dot - a2 * y * e_dot - a3 * pow(y, a8) * y_t_tcr;
	else
		return a1 * e_dot - a2 * y * e_dot;
}

DATA calc_val(DATA t, DATA tcr, int steps, DATA beg, DATA end, std::vector<DATA> y) {
	DATA td = t - tcr;
	if (td <= beg)return y[0];
	td -= beg;
	DATA h = (end - beg) / steps;
	int ind = int(floor(td / h));
	DATA rest = td - ind * h;
	DATA proc = rest / h;
	return y[ind] * proc + (1 - proc) * y[ind + 1];
}

void euler_2(DATA y_at_beg, DATA beg, DATA end, int steps, DATA yprim(DATA x, DATA y, DATA y_t_tcr,
	DATA a1, DATA a2, DATA a3, DATA a8, DATA p_cr, DATA t_cr
), DATA* a) {

	DATA h = (end - beg) / steps;

	std::cout << "\n\nPrzedzial [" << beg << "," << end << "] warunek pocz y0=" << y_at_beg << " krok (h) " << h;
	std::cout << std::endl;
	std::vector<DATA> ys;
	ys.resize(steps + 1);
	ys[0] = y_at_beg;

	DATA Z = e_dot * exp(Q / (R *tK));
	DATA p_cr = -a[10] + a[11] * pow(Z,a[9]);
	std::cout << "PCR" << p_cr << std::endl;
	DATA l = a[0] / pow(Z, a[12]);
	DATA a1 = 1 / (b*l);
	DATA a2 = a[1] * pow(e_dot,-a[8]) * exp(-a[2]/(R * tK));
	DATA tau = 1e6 * u * b * b / 2;
	DATA a3 = a[3] * tau / D * exp(-a[4] / (R * tK));
	DATA t_cr = 100000;
	bool crit = true;

	for (int i = 1; i <= steps; i++) {
		ys[i] = ys[i - 1] + h * yprim(beg + (i - 1) * h, ys[i - 1], calc_val(beg + (i - 1) * h, t_cr, steps, beg, end, ys),
			a1,a2,a3,a[7], p_cr, t_cr);
		std::cout << ys[i] << " " << std::endl;
		if (ys[i] > p_cr && crit) {
			t_cr = beg + (i - 1) * h;
			crit = false;
		}
	}


	std::cout << "\nEuler Wartosc funkcji w punkcie " << ys[ys.size() - 1];

	std::ofstream plik("data.txt");
	for (int i = 0; i < ys.size(); i++)
		plik << ys[i] << std::endl;

	plik.close();
}

