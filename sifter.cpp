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
#include <vector>
#include <string>

#include <omp.h>

#include "graph.hpp"

#define DEBUG_LOG   2

//! Находит факториал _top!. Если указан _bot, то вычисляется факториал _top!/_bot!.
uint fact(const uint _top, const uint _bot = 1);

/*
 * @brief Генерация заготовок графов
 *  Функция рекусивно создаёт заготовки графов, перебирая все комбинации того,
 *  куда смотрят последние (_i - _p) портов ввода.
 * 
 * @param _i        Порт ввода с которого начинаем перебор
 * @param _p        Число портов ввода
 * @param _g        Пустой граф (все узлы "0") с отметкой о занятости входного узла
 * @param answer    Сюда пишется результат
 */
void make_templates_graphs(
    uint _i,
    uint _p, 
    std::vector<std::pair<uint, bool> > _g, 
    std::vector<std::vector<uint> > &answer);

int main(int argc, char ** argv)
{
    using namespace std;

    if(argc < 5)
    {
        cerr << "Недостаточно аргументов" << endl;
        return 1; 
    }

    uint p = stoi(string(argv[1]));
    uint bs = stoi(string(argv[2]));
    uint dc = stoi(string(argv[3]));
    uint w = stoi(string(argv[4]));

    #if DEBUG_LOG >= 1
    cout << "p = " << p << endl
        << "bs = " << bs << endl
        << "dc = " << dc << endl
        << "w = " << w << endl;    
    #endif

    if(argc == 5)
    {
        cerr << "Матрица не введена" << endl;
        return 2; 
    }

    if(argc < 5 + p*p)
    {
        cerr << "Матрица введена неполностью" << endl;
        return 3; 
    }

    //! Матрица для просеивания
    graph::smatrix_t sM(p, vector<bool>(p));
    for(size_t row = 0; row < p; row++)
    for(size_t col = 0; col < p; col++)
        sM[col][row] = bool(stoi(string(argv[5 + row*p + col])));

    #if DEBUG_LOG >= 1
    cout << endl << "--Sift matrix--" << endl;
    for(size_t i = 0; i < p; ++i)
    {
        for(size_t j = 0; j < p; ++j)
        cout << sM[i][j] << '\t';

        cout << endl;
    }
    #endif

    int ompThreads;
    #pragma omp parallel
    if(omp_get_thread_num() == 0) ompThreads = omp_get_num_threads();

    #if DEBUG_LOG >= 1
    cout << "Number of threads: " << ompThreads << endl;
    #endif

    set<vector<uint> > templates;
    // Создание заготовок
    {
        uint T; //!< Число заготовок
        uint i = p + 1; //!< Порт ввода
        do
        {
            --i;
            T = fact(p + 2*(bs+dc+w), 2*(bs+dc+w) + i);
        } while(T < 10 * ompThreads && i != 0);

        #if DEBUG_LOG >= 1
        cout << "Templates to create (" << i << " to " << p << "): " << T << endl;
        #endif

        vector<vector<uint> > templates;
        vector<pair<uint, bool> > g(p + 2*(bs+dc+w), pair<uint, bool>(0, false));
        make_templates_graphs(i, p, g, templates);

        #if DEBUG_LOG >= 1
        cout << "Templates quantity: " << templates.size() << endl;
        #endif

        #if DEBUG_LOG >= 2
        cout << "Ten templates: " << endl;
        const size_t templShow = 10;
        for(size_t i = 0; i < templShow; ++i)
        {
            const size_t num = i * templates.size() / templShow;
            cout << i << ": ";
            for(size_t i = 0; i < templates[num].size(); ++i)
            cout << templates[num][i] << '\t';
            cout << endl;
        }
        #endif
    }

    #pragma omp parallel for schedule(guided)
    for(uint templ = 0; templ < templates.size(); ++templ)
    {
        graph g(p, bs, dc, w);
        g.set_sift_matrix(sM);


        g.set_edges(templNum_to_edges(templ));
        
    //     if(g.sifting())
    //     {
    //         vector<uint> we = g.get_edges();
    //         #pragma omp atomic
    //         cout.write((char*)we, we.size() * sizeof(we[0]));
    //     }
    }

    return 0;
}

uint fact(const uint _top, const uint _bot/* = 1*/)
{
    if(_top == _bot)
    return 1;
    else return _top * fact(_top - 1, _bot);
}

void make_templates_graphs(
    uint _i,
    uint _p, 
    std::vector<std::pair<uint, bool> > _g, 
    std::vector<std::vector<uint> > &answer)
{
    const size_t &s = _g.size();
	if (_i < _p)
	{
		for (uint to = 0; to < _g.size(); to++)//Номер узла, куда смотрит порт ввода
		{
			if (!_g[to].second)
			{
				_g[s - _p + _i].first = to;
				_g[to].second = true;
				make_templates_graphs(_i + 1, _p, _g, answer);
				_g[to].second = false;
			}
		}
	}
    else 
    {
        answer.push_back(std::vector<uint>(_g.size(), 0));
        for(size_t i = 0; i < _g.size(); ++i)
        answer.back()[i] = _g[i].first;
    }
}


