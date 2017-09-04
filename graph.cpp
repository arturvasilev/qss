#include <vector>
#include <string>
#include <sstream>

using namespace std;

//Перечисление типов однокубитовых операторов
enum operators_types
{
	beamsplitter,
	directCoupler,
	waveplate
};

//Содержит всю информацию об отдельном графе
struct graph
{
	/*
	p(orts) - число портов ввода-вывода
	cp (coupPlates) - число светоделительных пластинок
	dc (directionCoupler) - число направленных светоделителей
	w(aveplates) - число волновых пластинок (фазовращателей)
	*/
	graph(char ports = 0, char beamsplitters = 0, char directCouplers = 0, char waveplates = 0)
	{
		p = ports;
		bs = beamsplitters;
		dc = directCouplers;
		w = waveplates;
		q = bs + dc + w;
		size = p + 2 * (bs + dc + w);
		edges.resize(size);
		busy.resize(size, false);
		comb.resize(q);
		efficiency = 0;
		//Инициализация comb
		{
			for (char i = 0; i < bs; i++)			comb[i] = beamsplitter;
			for (char i = bs; i < bs + dc; i++)		comb[i] = directCoupler;
			for (char i = bs + dc; i < q; i++)		comb[i] = waveplate;
		}
		var.resize(bs + dc + 2 * w, 0);
	};

	//Направленный граф
	vector<char> edges;
	//Хранит в себе инфу о том, в какие из вершин уже смотрит ребро
	vector<bool> busy;
	//Хранит в себе комбинацию типов однокубитовых операторов
	vector<operators_types> comb;

	char size;//Размер направленного графа
	char p;//p(orts) - число портов ввода-вывода
	char q;//Число однокубитовых операторов
	char bs;//(beamsplitters) - число светоделительных пластинок
	char dc;//(direction couplers) - число направленных светоделителей
	char w;//(waveplates) - число волновых пластинок (фазовращателей)
	double efficiency;//Эффективность работы графа

	vector<double> var;//Внутренние параметры графа

	//Возвращает указатель на переменную из массива var[], соответствующей однокубитовому оператору comb->op[oper_num]
	//Если для данного однокубитового оператора нужно несколько переменных - то будет возвращёна ссылка на первую из них
	double* var_num(char oper_num)
	{
		//Вычислим номер необходимого оператора из массива var[]
		char var_num = 0;
		for (char i = 0; i < oper_num; i++)
			switch (comb[i])
			{
			case beamsplitter:
			case directCoupler: var_num++; break;
			case waveplate: var_num += 2; break;
			};
		return &(var[var_num]);
	};

	//Определяет какому типу однокубитового оператора принадлежит переменная с номером var_num
	operators_types oper_type(char var_num)
	{
		char oper_num = 0;
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
	};
	graph& operator= (const graph &other)
	{
		if (this == &other) return *this;
		graph(other.p, other.bs, other.dc, other.w);
		edges = other.edges;
		var = other.var;
		busy = other.busy;
		comb = other.comb;
		
		return *this; 
	};

	//Возвращает строку с текстовым представлением графа
	string print(void)
	{
		stringstream tmp;
		for(auto i : edges)
		{
			int e = i;
			tmp << e << ' ';
		}

		tmp << endl;

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
};