#ifndef GRAPH_CPP
#define GRAPH_CPP

#include "graph.hpp"

graph::graph(uint ports, uint beamsplitters, uint directCouplers, uint waveplates)
{
	p = ports;
	bs = beamsplitters;
	dc = directCouplers;
	w = waveplates;
	comb.resize(bs + dc + w);
	
	//Инициализация comb
	{
		for (uint i = 0; i < bs; i++)					comb[i] = beamsplitter;
		for (uint i = bs; i < bs + dc; i++)				comb[i] = directCoupler;
		for (uint i = bs + dc; i < comb.size(); i++)	comb[i] = waveplate;
	}

	var.resize(bs + dc + 2 * w, 0.5);

	targetMatrix.resize(p, std::vector<std::complex<double> >(p, 0.0));

	traj.resize(p, std::vector<std::set<std::set<uint> > >(p));

	translate.resize(4, std::vector<std::vector<uint> >(4, std::vector<uint>(4)));

	translate[0][0] = {1,1, 3,3};
	translate[0][1] = {1,1, 2,3};
	translate[0][2] = {0,1, 3,3};
	translate[0][3] = {0,1, 2,3};

	translate[1][0] = {1,1, 3,2};
	translate[1][1] = {1,1, 2,2};
	translate[1][2] = {0,1, 3,2};
	translate[1][3] = {0,1, 2,2};
	
	translate[2][0] = {1,0, 3,3};
	translate[2][1] = {1,0, 2,3};
	translate[2][2] = {0,0, 3,3};
	translate[2][3] = {0,0, 2,3};
	
	translate[3][0] = {1,0, 3,2};
	translate[3][1] = {1,0, 2,2};
	translate[3][2] = {0,0, 3,2};
	translate[3][3] = {0,0, 2,2};
};

//! Возвращает строку с текстовым представлением графа для ввода-вывода
const std::string graph::print(void)
{
	std::stringstream tmp;
	for(auto i : edges)
	{
		int e = i;
		tmp << e << ' ';
	}

	tmp << std::endl;

	for(auto i : comb)
	switch(i)
	{
		case beamsplitter: tmp << "bs "; break;
		case directCoupler: tmp << "dc "; break;
		case waveplate: tmp << "wp "; break;
		default:; 
	}
	return tmp.str();
} 

void graph::set_edges(const std::vector<uint> _edges)
{
	edges = _edges;
}

void graph::set_variables(const std::vector<double> &_var)
{
	var= _var;
}

double graph::get_deviation()
{
	if(!sift()) return __DBL_MAX__;

	double deviation = 0.;

	for(size_t i = 0; i < p; ++i)
	for(size_t j = 0; j < p; ++j)
	{
		std::set<std::set<uint> > &t = traj[i][j];
		
		if(!t.empty())
		for(auto it = t.begin(); it != t.end(); ++it)
		deviation += abs(
			traj_to_ampl(*it);
		);
	}

	return deviation;
}

void graph::set_target_matrix(const cmatrix_t &_tM)
{
	targetMatrix = _tM;
}

void graph::make_matrix_traj()
{
	for(size_t i = 2*p; i < edges.size(); ++i)
	paths(i, *this, this->traj);
}

double graph::get_deviation(const std::vector<double> &x, std::vector<double> &grad, void* f_data)
{
	graph * gp = reinterpret_cast<graph*>(f_data);
	gp->var = x;

	return gp->get_deviation();
}

void graph::paths(uint start, graph &g, traj_t &matrix, std::vector<uint> way)
{
	way.push_back(start);
	way.push_back(g.edges[start]);

	if (g.edges[start] < g.q * 2) //Если наш порт смотрит в однокубитовый оператор
	{
		paths((g.edges[start] / 2) * 2, g, matrix, way);
		paths((g.edges[start] / 2) * 2 + 1, g, matrix, way);
	}
	else//Если наш порт смотрит в порт вывода
	{
		//Запишем получившуюся траекторию в массив
		matrix[way.front() - g.q * 2][way.back() - g.q * 2].insert(way);
	}
};

bool graph::sift(const smatrix_t &_sM)
{
	if(traj[0][0].empty())
	make_matrix_traj();

	bool sift_ok = true;

	for(size_t i = 0; i < p; ++i)
	for(size_t j = 0; j < p; ++j)
	{
		std::vector<uint> a = translate[i][j];
		if(
			traj[a[0]][a[1]].empty() || traj[a[2]][a[3]].empty() &&
			traj[a[0]][a[3]].empty() || traj[a[2]][a[1]].empty()
		) {sift_ok = false; break;}
	}

	return sift_ok;
}

std::complex<double> graph::get_func(uint oper_num, uint in, uint out)
{
	switch (g.comb[oper_num])
	{
	case beamsplitter:
	{
		double *tmp = var_num(oper_num);
		if (in == 0 && out == 0) return sqrt(*tmp);}; else
		if (in == 0 && out == 1) return sqrt(complex<double>(1,0) - *tmp);}; else
		if (in == 1 && out == 0) return sqrt(complex<double>(1,0) - *tmp);}; else
		if (in == 1 && out == 1) return -sqrt(*tmp);};
		break;
	}
	case directCoupler:
	{
		double *tmp = g.var_num(oper_num);

		if (in == 0 && out == 0) 
		return 
			sqrt(*tmp);}; else
		if (in == 0 && out == 1) 
		return 
			sqrt((complex<double>)1 - *tmp)*exp(complex<double>(0, M_PI / 2));}; else
		if (in == 1 && out == 0) 
		return
			sqrt((complex<double>)1 - *tmp)*exp(complex<double>(0, M_PI / 2));}; else
		if (in == 1 && out == 1)
		return
			sqrt(*tmp);};
		break;
	}
	case waveplate:
	{
		double *phi = g.var_num(oper_num);
		double *alpha = (phi + 1);

		if (in == 0 && out == 0) 
		return 
			exp(complex<double>(0, *phi * 2 * M_PI)) * 
			pow(cos(*alpha * 2 * M_PI), 2) + pow(sin(*alpha * 2 * M_PI),2); else
		if (in == 0 && out == 1) 
		return 
			(exp(complex<double>(0,*phi * 2 * M_PI)) - (complex<double>)1) * 
			cos(*alpha * 2 * M_PI)*sin(*alpha * 2 * M_PI); else
		if (in == 1 && out == 0)
		return 
			(exp(complex<double>(0, *phi * 2 * M_PI)) - (complex<double>)1) * 
			cos(*alpha * 2 * M_PI)*sin(*alpha * 2 * M_PI); else
		if (in == 1 && out == 1) 
		return 
			exp(complex<double>(0, *phi * 2 * M_PI)) * 
			pow(sin(*alpha * 2 * M_PI), 2) + pow(cos(*alpha * 2 * M_PI), 2);
		break;
	}
	default: return 0;
	}
}

graph::operators_types graph::oper_type(uint var_num)
{
	uint oper_num = 0;
	while (var_num > 0)
		switch (comb[oper_num])
		{
		case beamsplitter:
		case directCoupler:
			var_num--; break;
		case waveplate:
			var_num -= 2;
		}
	return comb[oper_num];
}

graph& operator= (const graph &other)
{
	if (this == &other) return *this;
	graph(other.p, other.bs, other.dc, other.w);
	edges = other.edges;
	var = other.var;
	busy = other.busy;
	comb = other.comb;
	
	return *this; 
}

std::complex<double> traj_to_ampl(const std::set<uint> &_traj)
{
	std::complex<double> ret = 1.0;

	for(size_t i = 1; i < _traj.size() - 1; i += 2)
	ret *= get_func(_traj[i] / 2, _traj[i] % 2, _traj[i+1] % 2);

	return ret;
}

#endif //! GRAPH_CPP