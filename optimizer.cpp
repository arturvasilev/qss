#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <unistd.h>
#include <omp.h>
#include <string>
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

int main(int argc, char ** argv)
{
    using namespace std;

    if(argc == 1)
    {
        cerr << "Enter file name with graphs" << endl;
        return 1;
    }
    
    ifstream gfile(argv[1]);
    if(!gfile.is_open())
    {
        cerr << "Cannot open file" << endl;
        return 2;
    }
    
    size_t numGraphs = 0;
    {
        std::string line;

        while (std::getline(gfile, line))
            ++numGraphs;

        gfile.clear(ios_base::eofbit);
        gfile.seekg(0);
    }

    uint p, bs, dc, w;
    {
        string str;
        gfile >> str;
        p = stoi(str);
        gfile >> str;
        bs = stoi(str);
        gfile >> str;
        dc = stoi(str);
        gfile >> str;
        w = stoi(str);

        cout << p << '\t' << bs << '\t' << dc << '\t' << w << endl;
        
        --numGraphs;
    }
    const uint gSize = p + 2*(bs+dc+w);

    graph::cmatrix_t targetMatrix(p, std::vector<std::complex<double> >(p));
    {
        for(size_t i = 0; i < p; ++i)
        {
            for(size_t j = 0; j < p; ++j)
            {
                gfile >> targetMatrix[i][j];
                cout << targetMatrix[i][j] << '\t';
            }
            cout << endl;
        }

        numGraphs -= p;
    }

    vector<graph> best;

    double best_dev = __DBL_MAX__;

    #pragma omp parallel for schedule(guided)
    for(size_t i = 0; i < numGraphs - 1; ++i)
    {
        // cout << "Graph #" << i << endl;
        graph g(p, bs, dc, w);
        vector<uint> edges(p+2*(bs+dc+w));
        g.set_target_matrix(targetMatrix);

        #pragma omp critical(graphs)
        {
            for(uint i = 0; i < gSize; ++i)
            gfile >> edges[i];
        }

        g.set_edges(edges);
        NLopt(g, 1e-2);
        const double dev = g.get_deviation();

        #pragma omp critical(best)
        {
            if(dev < best_dev) 
            {
                best_dev = dev;
                cout << endl << "best.size() = " << best.size() << endl;
                cout << "Graph #" << i << " deviation: " << dev << endl;
                best.push_back(g);

                cout << "Edges:" << endl;
                for(auto j : best.back().get_edges())
                cout << j << '\t';
                cout << endl;
    
                cout << "Variables:" << endl << '\t';
                for(auto j : best.back().get_variables())
                cout << j << '\t';
                cout << endl;

                cout << "\tMatrix amplitude:" << endl;
                graph::cmatrix_t MAmp = best.back().get_matrix_amplitude();
                for(size_t i = 0; i < p; ++i)
                {
                    cout << "\t\t";
                    for(size_t j = 0; j < p; ++j)
                    cout << MAmp[i][j] << '\t';
    
                    cout << endl;
                }
    
                cout << "\tMatrix truth:" << endl;
                graph::cmatrix_t MTruth = best.back().get_matrix_truth();
                for(size_t i = 0; i < p; ++i)
                {
                    cout << "\t\t";
                    for(size_t j = 0; j < p; ++j)
                    cout << MTruth[i][j] << '\t';
    
                    cout << endl;
                }
            }
        }

        #pragma omp critical(stderr)
        {
            static uint toShow = 50;
            static uint processed = 0;
            if(++processed % (numGraphs/toShow) == 0)
            cerr << round(100 * float(processed) / numGraphs) << "% graphs" << endl;
        }
    }

    // // Графы упорядочены по убыванию deviation
    // // Необходимо обратить этот порядок
    // {
    //     vector<graph> best_rev(best.size());
    //     for(size_t i = 0; i < best.size(); ++i)
    //     best_rev[i] = best[best.size() - 1 - i];

    //     swap(best, best_rev);
    // }

    // // Вывод на экран лучших десяти графов
    // {
    //     for(size_t i = 0; i < 10; ++i)
    //     {
    //         cout << best[i].get_deviation() << '\t';
    //         for(auto j : best[i].get_edges())
    //         cout << j << '\t';

    //         cout << endl;

    //         cout << '\t';
    //         for(auto j : best[i].get_variables())
    //         cout << j << '\t';
    //         cout << endl;

    //         graph::cmatrix_t MTruth = best[i].get_matrix_truth();
    //         for(size_t i = 0; i < p; ++i)
    //         {
    //             for(size_t j = 0; j < p; ++j)
    //             cout << MTruth[i][j] << '\t';

    //             cout << endl;
    //         }
    //     }
    // }

    return 0;
}