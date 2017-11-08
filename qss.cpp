/*
 * Смысл жизни данного приложения в том, чтобы
 * по стандартному потоку принять предельные размеры графов и
 * булеву матрицу PxP (P - число портов ввода-вывода),
 * по которой будет выполняться просеивание.
 * 
 * В процессе своей работы графы, которые пройдут все проверки будут выведены 
 * на стандартный вывод для дальнейшей оптимизации на хостовой системе.
 * 
 * В стандартный поток ошибок выводится текущий прогресс.
 * 
 */

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include <string>
#include <vector>
#include <complex>
#include <algorithm>

#include <nlopt.hpp>
#include "graph.hpp"

double target_function(const std::vector<double> &x, std::vector<double> &grad, void * data);

/*
 * @brief Реализация NLopt. 
 *
 * @param g         направленный граф
 * @param L         логическая матрица из функций
 * @param restrictions массив условий равенства
 * @param eps       точность удовлетворения условиям равенства
 */
void NLopt(graph &_g, double _eps)
{
    const uint v = _g.get_variables().size();

    // std::cout << "Setting glob_problem" << std::endl;
    // //! Поиск глобального оптимума, без производных
    nlopt::opt glob_problem(nlopt::AUGLAG, v);
    
    // std::cout << "Setting min_objective" << std::endl;
    glob_problem.set_min_objective(&target_function, (void*)&_g);

    //! Устанавливаем границы изменения переменных
    std::vector<double> lb(v, 0), ub(v, 1);
    glob_problem.set_lower_bounds(lb);
    glob_problem.set_upper_bounds(ub);

    //Задаём конечную точность установления переменных
    glob_problem.set_xtol_abs(_eps);

    {
        // std::cout << "Setting loc_problem" << std::endl;
        //! Локальный оптимизатор
        nlopt::opt loc_problem(nlopt::LN_COBYLA, v);
        loc_problem.set_xtol_abs(_eps);
        //ПРЕЖДЕ локальный оптимизатор надо конфигурировать ДО того как 
        //передать его для _копирования_ глобальному. 
        glob_problem.set_local_optimizer(loc_problem);
    }

    double result;
    std::vector<double> grad(v);
    std::vector<double> x(v, 0.5);

    glob_problem.set_maxtime(1e-2);
    // std::cout << "Optimizing..." << std::endl;
    nlopt::result res = glob_problem.optimize(x, result);
}

bool compare_graph (graph &_a, graph &_b) { return (_a.get_deviation() < _b.get_deviation()); }

double target_function(const std::vector<double> &x, std::vector<double> &grad, void * data)
{
    graph *g = reinterpret_cast<graph *>(data);

    g->set_variables(x);
    
    double ret = g->get_deviation();
    // std::cout << "\tDeviation = " << ret << std::endl;
    return ret;
}