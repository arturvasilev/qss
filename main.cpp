// Optim_monte-karlo.cpp : Defines the entry point for the console application.
//Данное приложение делает полный перебор квантовых схем на основе направленных графов, фильтрует их по матрице амплитуд
//и далее для отфильтрованных схем находит глобальный оптимум внутренних параметров схемы (напр., соотношения светоделителей),
//при которых схема работает с наибольшей вероятностью.

#include <vector>
#include <set>
#include <complex>
#include <math.h>
#include <cmath>
#include <functional>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <algorithm>
#include <fstream>
#include <mpi.h>
#include <nlopt.hpp>

#include "graph.cpp"

using namespace std;

#define myfunc function<complex<double>()>

//Рекурсивно находит все возможные пути
void paths(char start, graph &g, vector<vector<set<vector<char> > > > &matrix, vector<char> way)
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

//Возвращает функцию, определяя тип оператора по номеру oper_num из графа g, а номера портов ввода-вывода по in и out
myfunc get_myfunc(graph &g, char oper_num, char in, char out)
{
	switch (g.comb[oper_num])
	{
	case beamsplitter:
	{
		double *tmp = g.var_num(oper_num);
		if (in == 0 && out == 0) return [tmp](void)->complex<double> {return sqrt(*tmp);}; else
		if (in == 0 && out == 1) return [tmp](void)->complex<double> {return sqrt(complex<double>(1,0) - *tmp);}; else
		if (in == 1 && out == 0) return [tmp](void)->complex<double> {return sqrt(complex<double>(1,0) - *tmp);}; else
		if (in == 1 && out == 1) return [tmp](void)->complex<double> {return -sqrt(*tmp);};
		break;
	}
	case directCoupler:
	{
		double *tmp = g.var_num(oper_num);

		if (in == 0 && out == 0) return [tmp](void)->complex<double> {return sqrt(*tmp);}; else
		if (in == 0 && out == 1) return [tmp](void)->complex<double> {return sqrt((complex<double>)1 - *tmp)*exp(complex<double>(0, M_PI / 2));}; else
		if (in == 1 && out == 0) return [tmp](void)->complex<double> {return sqrt((complex<double>)1 - *tmp)*exp(complex<double>(0, M_PI / 2));}; else
		if (in == 1 && out == 1) return [tmp](void)->complex<double> {return sqrt(*tmp);};
		break;
	}
	case waveplate:
	{
		double *phi = g.var_num(oper_num);
		double *alpha = (phi + 1);

		if (in == 0 && out == 0) return [phi, alpha](void)->complex<double>
			{return exp(complex<double>(0, *phi * 2 * M_PI))*pow(cos(*alpha * 2 * M_PI), 2) + pow(sin(*alpha * 2 * M_PI),2);}; else
		if (in == 0 && out == 1) return [phi, alpha](void)->complex<double>
			{return (exp(complex<double>(0,*phi * 2 * M_PI)) - (complex<double>)1)*cos(*alpha * 2 * M_PI)*sin(*alpha * 2 * M_PI);}; else
		if (in == 1 && out == 0) return [phi, alpha](void)->complex<double>
			{return (exp(complex<double>(0, *phi * 2 * M_PI)) - (complex<double>)1)*cos(*alpha * 2 * M_PI)*sin(*alpha * 2 * M_PI);}; else
		if (in == 1 && out == 1) return [phi, alpha](void)->complex<double>
			{return exp(complex<double>(0, *phi * 2 * M_PI))*pow(sin(*alpha * 2 * M_PI), 2) + pow(cos(*alpha * 2 * M_PI), 2);};
		break;
	}
	default: return []()->complex<double> {return 0;};
	}
}

//Возвращает матрицу амплитуд по матрице траекторий. Граф g нужен для определения типа однокубитового элемента и числа тех или иных элементов
vector<vector<myfunc> > make_matrix_amplitude(vector<vector<set<vector<char> > > > &matrix, graph &g)
{
	vector<vector<myfunc> > ans(g.p, vector<myfunc>(g.p));;//Функции будем формировать здесь

	for (char i = 0; i < g.p; i++)
		for (char j = 0; j < g.p; j++)
		{
			if (!matrix[i][j].empty())
			{
				while (!matrix[i][j].empty())
				{
					vector<char> one_path = *(matrix[i][j].begin());//Отдельный путь
					matrix[i][j].erase(matrix[i][j].begin());
					myfunc summand = nullptr;//Отдельное слагаемое
					for (char k = 1; k < (char)one_path.size() - 1; k += 2)
					{
						char oper_num = one_path[k] / 2;//Порядковый номер однокубитового оператора из g.comb
						char in = one_path[k] % 2;//Номер порта ввода
						char out = one_path[k + 1] % 2;//Номер порта вывода
						if (summand == nullptr)
						{
							//В этой точке мы знаем номер однокубитового порта и его порты ввода-вывода
							summand = get_myfunc(g, oper_num, in, out);
						}
						else
						{
							myfunc old_myfunc = summand;
							myfunc new_myfunc = get_myfunc(g, oper_num, in, out);
							summand = [old_myfunc, new_myfunc]()->complex<double> {return old_myfunc() * new_myfunc();};
						}
					}//end for(k)

					//В данной точке мы сформировали функцию, которая представляет собой отдельное слагаемое
					if (ans[i][j] == nullptr) ans[i][j] = summand;
					else
					{
						myfunc old_myfunc = ans[i][j];
						ans[i][j] = [old_myfunc, summand]()->complex<double> {return old_myfunc() + summand();};
					}
				}
			}
			else ans[i][j] = []()->complex<double> {return 0;};
		}
	return ans;
};

//Структура, необходимая для работы конвертера convert
struct convert_dummy {
	myfunc *function;
	int number;
	vector<double> *x;
};

//Конвертор функций в формат под Nlopt
double convert(const vector<double> x, vector<double> &grad, void *data)
{
	convert_dummy &d = reinterpret_cast<convert_dummy*>(data);
	d->x = x;
	return (d->function)();
}

//Реализация NLopt
//Возвращает максимальное значение целевой функции aim
//g - направленный граф
//L - логическая матрица из функций
//restrictions - массив условий равенства
//eps - точность удовлетворения условиям равенства
double NLopt (graph &g, myfunc obj_myfunc, vector<myfunc> &restrictions, double eps)
{
    //Зададим точность нахождения 
    //const double eps = 1e-3;
    int dim = g.var.size();
	g.var.assign(dim, 0.5);
    //Поиск глобальный оптимума, без производных
    nlopt::opt glob_problem(nlopt::AUGLAG, dim);

	vector<dummy> tmp(restrictions.size()+1);
	{
		for(size_t i = 0; i < tmp.size(); i++) 
			{tmp[i].x = &g.var; tmp[i].number = i;}
		for(size_t i = 0; i < tmp.size()-1; i++)
			tmp[i].function = &restrictions[i];
		tmp.back().function = obj_myfunc;
	}
    glob_problem.set_max_objective(convert, &tmp.back(), eps);
    //Устанавливаем границы изменения переменных
    vector<double> lb(dim, 0), ub(dim, 1);
    glob_problem.set_lower_bounds(lb);
    glob_problem.set_upper_bounds(ub);
	for (size_t i = 0; i < tmp.size()-1; i++)
		glob_problem.add_equality_constraint(convert, &tmp[i], eps);
			
	cout<<"I'm HERE**********";
    //Задаём конечную точность установления переменных
    glob_problem.set_xtol_abs(eps);

    {
        nlopt::opt loc_problem(nlopt::LN_COBYLA, dim);
        loc_problem.set_xtol_abs(eps);
        //ПРЕЖДЕ локальный оптимизатор надо конфигурировать ДО того как 
        //передать его для _копирования_ глобальному. 
        glob_problem.set_local_optimizer(loc_problem);
    }
}


//Данная рекурсия может получить на входе пустой граф, или его заготовку
//Возвращает наилучший граф из всех возможных для данной заготовки
//me - узел, с которого мы сейчас стартуем
void choose_best_graph(char me, graph g, graph &g_best)
{
	//Необходимо достроить направленный граф g
	if (me < g.q * 2)//Если мы сейчас стартуем из однокубитового оператора
	{
		for (char i = (me / 2) * 2 + 2; i < g.edges.size(); i++)//Узел, куда смотрит порт выхода однокубитового оператора
		{
			if (!g.busy[i])
			{
				g.busy[i] = true;
				g.edges[me] = i;
				choose_best_graph(me + 1, g, g_best);
				g.busy[i] = false;
			}
		}
	}
	else//В этой точке мы стартуем из порта ввода
	{
		if (me < g.edges.size())
			//Порт ввода может смотреть в любой узел графа
			for (int i = 0; i < g.edges.size(); i++)
			{
				if (!g.busy[i])
				{
					g.busy[i] = true;
					g.edges[me] = i;
					choose_best_graph(me + 1, g, g_best);
					g.busy[i] = false;
				}
			}
		else
		{
			//В этой ветке мы уже прошлись по всем исходящим вершинам графа.
			//Теперь в g полностью запоненный массив, представляющий собой направленный граф. Но он ещё абстрактен от светоделителей и фазовращателей

			//Накостыляем условия, которые скажут нам есть ли связь между нужным портом ввода и нужным(и) выводами,
			//и по этим условиям отсеять неподходящие графы
			{
				vector<vector<set<vector<char>>>> matrix(g.p, vector<set<vector<char>>>(g.p));//Матрица траекторий

				//Заполним матрицу траекторий
				//Стартанём алгоритм поиска траекторий из всех портов ввода
				//Результат будет в matrix
				for (int i = g.q * 2; i < g.size; i++) paths(i, g, matrix, vector<char>());

				//Теперь у нас есть матрица траекторий для сгенерированного графа. Теперь определим удовлетворяет ли этот граф условию.
				if (matrix[0][2].empty() && matrix[0][3].empty() && //Первая версия
					!matrix[1][2].empty() && !matrix[1][3].empty() &&
					!matrix[0][0].empty() && matrix[0][1].empty() && //Начало второй версии
					matrix[1][0].empty() && !matrix[1][1].empty() &&
					matrix[2][0].empty() && !matrix[2][1].empty() && //Начало четвёртой версии
					!matrix[2][2].empty() && !matrix[2][3].empty() &&
					matrix[3][0].empty() && !matrix[3][1].empty() &&
					!matrix[3][2].empty() && !matrix[3][3].empty())
				{
					//В этой точке у нас есть граф с матрицей траекторий, который удовлетворяет всем условиям
					//Теперь переберём все комбинации типов однокубитовых операторов

					set<vector<operators_types>> memory_oper_comb;//Все комбинации однокубитовых операторов

					do
					{
						if (memory_oper_comb.find(g.comb) == memory_oper_comb.end())
						{
							//В данной точке у нас уникальный заполненный граф g, с матрицей траекторий matrix, с уникальной комбинацией однокубитовых операторов
							memory_oper_comb.insert(g.comb);//Сохраним эту комбинацию, чтобы в будущем не повторяться

							//Теперь получим матрицу амплитуд из функций
							vector<vector<myfunc>> MAmpl = make_matrix_amplitude(matrix, g);

							//Теперь из этой матрицы амплитуд необходимо получить логическую матрицу, которая размером 4x4
							myfunc L[4][4];
							{
								L[0][0] = [&MAmpl]() {return MAmpl[1][1]() * MAmpl[3][3]() + MAmpl[1][3]() * MAmpl[3][1]();};
								L[0][1] = [&MAmpl]() {return MAmpl[1][1]() * MAmpl[2][3]() + MAmpl[1][3]() * MAmpl[2][1]();};
								L[0][2] = [&MAmpl]() {return MAmpl[0][1]() * MAmpl[3][3]() + MAmpl[0][3]() * MAmpl[3][1]();};
								L[0][3] = [&MAmpl]() {return MAmpl[0][1]() * MAmpl[2][3]() + MAmpl[0][3]() * MAmpl[2][1]();};

								L[1][0] = [&MAmpl]() {return MAmpl[1][1]() * MAmpl[3][2]() + MAmpl[1][2]() * MAmpl[3][1]();};
								L[1][1] = [&MAmpl]() {return MAmpl[1][1]() * MAmpl[2][2]() + MAmpl[1][3]() * MAmpl[2][1]();};
								L[1][2] = [&MAmpl]() {return MAmpl[0][1]() * MAmpl[3][2]() + MAmpl[0][2]() * MAmpl[3][1]();};
								L[1][3] = [&MAmpl]() {return MAmpl[0][1]() * MAmpl[2][2]() + MAmpl[0][2]() * MAmpl[2][1]();};

								L[2][0] = [&MAmpl]() {return MAmpl[1][0]() * MAmpl[3][3]() + MAmpl[1][3]() * MAmpl[3][0]();};
								L[2][1] = [&MAmpl]() {return MAmpl[1][0]() * MAmpl[2][3]() + MAmpl[1][3]() * MAmpl[2][0]();};
								L[2][2] = [&MAmpl]() {return MAmpl[0][0]() * MAmpl[3][3]() + MAmpl[0][3]() * MAmpl[3][0]();};
								L[2][3] = [&MAmpl]() {return MAmpl[0][0]() * MAmpl[2][3]() + MAmpl[0][3]() * MAmpl[2][0]();};

								L[3][0] = [&MAmpl]() {return MAmpl[1][0]() * MAmpl[3][2]() + MAmpl[1][2]() * MAmpl[3][0]();};
								L[3][1] = [&MAmpl]() {return MAmpl[1][0]() * MAmpl[2][2]() + MAmpl[1][2]() * MAmpl[2][0]();};
								L[3][2] = [&MAmpl]() {return MAmpl[0][0]() * MAmpl[3][2]() + MAmpl[1][2]() * MAmpl[3][0]();};
								L[3][3] = [&MAmpl]() {return MAmpl[0][0]() * MAmpl[2][2]() + MAmpl[0][2]() * MAmpl[2][0]();};
							}

							//Теперь необходимо составить список условий:
							vector<myfunc> restrictions;//Ограничения равенства " = 0 " - все они должны быть равны нулю
							{
								//унитарность
								{
									//Перебор всех пар столбцов.
									for (int k1 = 0; k1 < g.p - 1; k1++)//Номер первого столбца пары
										for (int k2 = k1 + 1; k2 < g.p; k2++)//Номер второго столбца пары
										{
											restrictions.push_back(
												[&MAmpl, k1, k2, g]()->complex<double> {
												complex<double> val;
												for (char i = 0; i < g.p; i++)
													val += MAmpl[i][k1]() * MAmpl[i][k2]();
												return val;
											});//end vector::push_back()
										}
								}//конец унитарность

								 //нулевые элементы в логической матрице
								{
									restrictions.push_back(L[0][0]);
									restrictions.push_back(L[0][2]);
									restrictions.push_back(L[0][3]);

									restrictions.push_back(L[1][1]);
									restrictions.push_back(L[1][2]);
									restrictions.push_back(L[1][3]);

									restrictions.push_back(L[2][0]);
									restrictions.push_back(L[2][1]);
									restrictions.push_back(L[2][3]);

									restrictions.push_back(L[3][0]);
									restrictions.push_back(L[3][1]);
									restrictions.push_back(L[3][2]);
								}//конец нулевые элементы логической матрицы

								 //Равенство действующих элементов логической матрицы
								{
									restrictions.push_back([&L]()->complex<double> {return L[0][1]() - L[1][0]();});
									restrictions.push_back([&L]()->complex<double> {return L[1][0]() - L[2][2]();});
									restrictions.push_back([&L]()->complex<double> {return L[2][2]() - L[3][3]();});
								}//конец равенство действующих элементов между собой
							}

							{

							}
							//В данной точке у нас есть функции для всех условий равенства и полностью собранный граф - осталось только найти глобальный оптимум
							if (g_best.efficiency < NLopt(g, L[3][3], restrictions, 1e-2))
							{
								cout << '.'; cout.flush();
								g_best = g;
							}
						}
					} while (next_permutation(g.comb.begin(), g.comb.end()));

					//В данной точке у нас есть граф g_best с оптимальными внутренними параметрами, при который он работает с вероятностью best
				}//Конец обработки графа, удовлетворяющего условиям матрицы траекторий
			}
		}
	}//Конец обработки сгенерированного направленного графа
}

//Функция возвращает заготовки графов, перебирая все комбинации того, куда смотрят первые count-портов ввода
void make_templates_graphs(char count, graph g, vector<graph> &answer, char deep = 0)
{
	if (deep < count)
	{
		for (char to = 0; to < g.edges.size(); to++)//Номер узла, куда смотрит порт ввода
		{
			if (!g.busy[to])
			{
				g.edges[g.q * 2 + deep] = to;
				g.busy[to] = true;
				make_templates_graphs(count, g, answer, deep + 1);
				g.busy[to] = false;
			}
		}
	}
	else answer.push_back(g);
}

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
    int MPI_size;	
	if (MPI_Comm_size(MPI_COMM_WORLD, &MPI_size) != MPI_SUCCESS) return EXIT_FAILURE;
    int MPI_rank;
	if (MPI_Comm_rank(MPI_COMM_WORLD, &MPI_rank) != MPI_SUCCESS) return EXIT_FAILURE;
	if (!MPI_rank) cout << "Processes: " << MPI_size << endl;
     	
    //Кусок кода, созданный для отладки написанного выше кода
	if (true)
	{
		graph g(6, 5, 0, 0);
		{
			char tmp_graph[] = { 10, 15, 4, 6, 8, 11, 9, 14, 12, 13, 0, 5, 2, 3, 7, 1 };
			for (char i = 0; i < g.size; i++)
			{
				g.edges[i] = tmp_graph[i];
				g.busy[i] = true;
			}
		}

		vector<vector<set<vector<char>>>> matrix_traj(g.p, vector<set<vector<char>>>(g.p));
		for (char i = 0; i < g.p; i++) paths(i + g.q * 2, g, matrix_traj, vector<char>());

		vector<vector<myfunc>> MAmpl = make_matrix_amplitude(matrix_traj, g);

		for (char i = 0; i < g.p; i++)
		{
			cout << endl;
			cout << setw(3) << real(MAmpl[i][0]());
			for (char j = 1; j < g.p; j++) cout << '\t' << setw(3) << real(MAmpl[i][j]());
		}
		myfunc L[4][4];
		{
			L[0][0] = [&MAmpl]() {return MAmpl[1][1]() * MAmpl[3][3]() + MAmpl[1][3]() * MAmpl[3][1]();};
			L[0][1] = [&MAmpl]() {return MAmpl[1][1]() * MAmpl[2][3]() + MAmpl[1][3]() * MAmpl[2][1]();};
			L[0][2] = [&MAmpl]() {return MAmpl[0][1]() * MAmpl[3][3]() + MAmpl[0][3]() * MAmpl[3][1]();};
			L[0][3] = [&MAmpl]() {return MAmpl[0][1]() * MAmpl[2][3]() + MAmpl[0][3]() * MAmpl[2][1]();};

			L[1][0] = [&MAmpl]() {return MAmpl[1][1]() * MAmpl[3][2]() + MAmpl[1][2]() * MAmpl[3][1]();};
			L[1][1] = [&MAmpl]() {return MAmpl[1][1]() * MAmpl[2][2]() + MAmpl[1][3]() * MAmpl[2][1]();};
			L[1][2] = [&MAmpl]() {return MAmpl[0][1]() * MAmpl[3][2]() + MAmpl[0][2]() * MAmpl[3][1]();};
			L[1][3] = [&MAmpl]() {return MAmpl[0][1]() * MAmpl[2][2]() + MAmpl[0][2]() * MAmpl[2][1]();};

			L[2][0] = [&MAmpl]() {return MAmpl[1][0]() * MAmpl[3][3]() + MAmpl[1][3]() * MAmpl[3][0]();};
			L[2][1] = [&MAmpl]() {return MAmpl[1][0]() * MAmpl[2][3]() + MAmpl[1][3]() * MAmpl[2][0]();};
			L[2][2] = [&MAmpl]() {return MAmpl[0][0]() * MAmpl[3][3]() + MAmpl[0][3]() * MAmpl[3][0]();};
			L[2][3] = [&MAmpl]() {return MAmpl[0][0]() * MAmpl[2][3]() + MAmpl[0][3]() * MAmpl[2][0]();};

			L[3][0] = [&MAmpl]() {return MAmpl[1][0]() * MAmpl[3][2]() + MAmpl[1][2]() * MAmpl[3][0]();};
			L[3][1] = [&MAmpl]() {return MAmpl[1][0]() * MAmpl[2][2]() + MAmpl[1][2]() * MAmpl[2][0]();};
			L[3][2] = [&MAmpl]() {return MAmpl[0][0]() * MAmpl[3][2]() + MAmpl[1][2]() * MAmpl[3][0]();};
			L[3][3] = [&MAmpl]() {return MAmpl[0][0]() * MAmpl[2][2]() + MAmpl[0][2]() * MAmpl[2][0]();};
		}

		cout << endl;
		for (char i = 0; i < 4; i++)
		{
			cout << endl << setw(3) << real(L[i][0]());
			for (char j = 1; j < 4; j++) cout << '\t' << setw(3) << real(L[i][j]());
		}
		cout << endl << endl;
		graph g_best(g);
		choose_best_graph(g.size, g, g_best);
		cout << endl << "Variables:" << endl;
		for (int i = 0; i < g_best.var.size(); i++) cout << "var[" << i << "]=" << g_best.var[i] << endl;

		cout << "\nLogic matrix:";
		for (char i = 0; i < 4; i++)
		{
			cout << endl << setw(3) << abs(L[i][0]());
			for (char j = 1; j < 4; j++) cout << '\t' << setw(3) << abs(L[i][j]());
		}
		cout << endl;
	}
	else
	//Основной боевой алгоритм
	{
		for (char p = 6; p <= 6; p++)
			for (char bs = 5; bs <= 5; bs++)
				for (char dc = 0; dc <= 0; dc++)
					for (char w = 0; w <= 0; w++)
					{
						//Создадим массив заготовок - переберём все варианты куда могут смотреть первые три порта ввода
						//это 16!/13! = 3360 для p = 6 и q = 5
						{
							if (!MPI_rank) cout << "Making templates..." << endl;
							graph g(p, bs, dc, w);
							vector<graph> templates_graphs;
							make_templates_graphs(3, g, templates_graphs);
							if (!MPI_rank) cout << "Made " << templates_graphs.size() << " graphs" << endl;
							MPI_Barrier(MPI_COMM_WORLD);

                            //start master
                            if (MPI_rank == 0)
                            {
                                vector<double> inbuf(MPI_size, MPI_PROC_NULL);
								double best_eff = 0;
                                
								MPI_Request reqs[MPI_size];
                                for(int i = 0; i < MPI_size; i++)
                                MPI_Irecv(&inbuf[i], 1, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, &reqs[i]);

                                for(int i = 0; i < templates_graphs.size(); i++)
                                {
									cout << i << ' ';
									cout.flush();
									int num_task;
									MPI_Status status;
									MPI_Waitany(MPI_size, reqs, &num_task, &status);
									best_eff = max(inbuf[num_task], best_eff);

									MPI_Send(&i, 1, MPI_INT, num_task, 0, MPI_COMM_WORLD);
                                	MPI_Irecv(&inbuf[num_task], 1, MPI_DOUBLE, num_task, MPI_ANY_TAG, MPI_COMM_WORLD, &reqs[num_task]);
									
                                    if (i % MPI_size == 0 && i != 0)//Каждый раз, когда мы раскидаем на все потоки новые заготовки
                                    {
                                        //Мы записываем текущий прогресс в текстовый файл
                                        cout << "Persentage: " << (double)i / templates_graphs.size() * 100 << "%" << endl;
                                        cout << "Best: " << best_eff << endl;
										cout.flush();
                                    }
                                }//end master

								int i = -1;
								MPI_Bcast(&i, 1, MPI_INT, 0, MPI_COMM_WORLD);
								MPI_Barrier(MPI_COMM_WORLD);
                            }
                            else
                            //start slaves 
                            {
                                graph loc_best_graph;
                                //Номер заготовки графа для обработки
                                do{
                                    MPI_Send(&loc_best_graph.efficiency, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
                                    int templ_number;
									MPI_Status status;
                                    MPI_Recv(&templ_number, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

									cout << templates_graphs[templ_number].print() << endl;
									
                                    if (templ_number < 0)
                                    {
										MPI_Barrier(MPI_COMM_WORLD);
                                        MPI_Finalize();
										cout << loc_best_graph.print() << endl;
                                        return EXIT_SUCCESS;
                                    } else choose_best_graph(0, templates_graphs[templ_number], loc_best_graph);
                                } while(true);
                            }//end slaves
						}
					}
	}//Конец реализации основного боевого алгоритма

    MPI_Finalize();
	return EXIT_SUCCESS;
}
