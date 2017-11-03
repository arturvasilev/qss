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
#include <unistd.h>
#include <stdlib.h>

#include <omp.h>

#include "graph.hpp"

#define SIFTER_DEBUG_LOG  0
u_int64_t   graphs_generated = 0;

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

//Данная рекурсия может получить на входе пустой граф, или его заготовку
//Возвращает наилучший граф из всех возможных для данной заготовки
//me - узел, с которого мы сейчас стартуем
//_e    Заготовка графа
//_sP   Порт ввода, от которого все последующие порты ввода (включая _sP) инициализированы
//_fP   Число портов ввода-вывода
void print_sifted_graphs(
    uint me, 
    std::vector<uint> _e, 
    std::vector<bool> _busy,
    const uint _sP, 
    const uint _fP,
    graph &_g,
    const graph::smatrix_t &_sM
);

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

    cerr << "There must be" << endl;
    const u_int64_t graphs_to_generate = fact(p) * pow(p*(p-1), bs+dc+w);
    cerr << '\t' << float(graphs_to_generate) << " graphs" << endl;

    #if SIFTER_DEBUG_LOG >= 1
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

    #if SIFTER_DEBUG_LOG >= 1
    cout << endl << "--Sift matrix--" << endl;
    for(size_t i = 0; i < p; ++i)
    {
        for(size_t j = 0; j < p; ++j)
        cout << sM[i][j] << '\t';

        cout << endl;
    }
    #endif

    #if SIFTER_DEBUG_LOG >= 3
    omp_set_num_threads(1);
    #endif

    int ompThreads;
    #pragma omp parallel
    if(omp_get_thread_num() == 0) ompThreads = omp_get_num_threads();

    #if SIFTER_DEBUG_LOG >= 1
    cout << "Number of threads: " << ompThreads << endl;
    #endif

    vector<vector<uint> > templates;
    // Создание заготовок
    uint startPort = p + 1; //!< Порт ввода
    {
        uint T; //!< Число заготовок
        do
        {
            --startPort;
            T = fact(p + 2*(bs+dc+w), 2*(bs+dc+w) + startPort);
        } while(T < 10 * ompThreads && startPort != 0);

        #if SIFTER_DEBUG_LOG >= 1
        cout << "Templates to create (" << startPort << " to " << p << "): " << T << endl;
        #endif

        templates.reserve(T);
        vector<pair<uint, bool> > g(p + 2*(bs+dc+w), pair<uint, bool>(0, false));
        make_templates_graphs(startPort, p, g, templates);

        #if SIFTER_DEBUG_LOG >= 1
        cout << "Templates quantity: " << templates.size() << endl;
        #endif

        #if SIFTER_DEBUG_LOG >= 2
        const size_t templShow = 10;
        cout << templShow << " templates: " << endl;
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

    cout << p << '\t' << bs << '\t' << dc << '\t' << w << endl;
    #pragma omp parallel for schedule(guided)
    for(size_t templ = 0; templ < templates.size(); ++templ)
    {
        //! Заготовка
        vector<uint> &T = templates[templ];
        graph g(p, bs, dc, w);
        
        vector<uint> e(T);

        vector<bool> busy(T.size(), false);
        for(size_t i = 2*(bs+dc+w)+startPort; i < T.size(); ++i)
        busy[e[i]] = true;

        print_sifted_graphs(0, e, busy, startPort, p, g, sM);

        #pragma omp critical(stderr)
        {
            static uint toShow = 50;
            static uint processed = 0;
            if(++processed % (templates.size()/toShow) == 0)
            cerr << round(100 * float(processed) / templates.size()) << "% templates (" 
                << round(100 * float(graphs_generated)/graphs_to_generate) << "% graphs generated)" << endl;
        }
    }

    cerr << "Generated graphs: " << graphs_generated << endl;

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

void print_sifted_graphs(
    uint me, 
    std::vector<uint> _e, 
    std::vector<bool> _busy,
    const uint _sP, 
    const uint _fP,
    graph &_g,
    const graph::smatrix_t &_sM)
{
    #if SIFTER_DEBUG_LOG >= 3
    #pragma omp critical(stdout)
    {
        std::cout << "templ #" << omp_get_thread_num() << " recursive: ";
        for(auto i : _e)
        std::cout << i << '\t';
        std::cout << std::endl;
    }
    #endif

    const std::vector<bool> BUSY_TRUE(_busy.size(), true);
    if(_busy == BUSY_TRUE)
    {
        // В этой ветке мы уже прошлись по всем исходящим вершинам графа.
        // Теперь в _e полностью запоненный массив, представляющий собой направленный граф.
        #pragma omp atomic
        ++graphs_generated;        

        _g.set_edges(_e);
        
        #if SIFTER_DEBUG_LOG >= 3
        #pragma omp critical(stdout)
        {
            std::cout << "thread #" << omp_get_thread_num() << " try: ";
            for(auto i : _e)
            std::cout << i << '\t';
            std::cout << std::endl;
        }
        #endif

        if(_g.sift(_sM))
        {
            #pragma omp critical(stdout)
            {
                for(auto i : _e)
                std::cout << i << '\t';
                std::cout << std::endl;
            }

        }
    }
    else
	//Необходимо достроить направленный граф g
	if (me < (_e.size() - _fP))//Если мы сейчас стартуем из однокубитового оператора
	{
        //! Узел, куда смотрит порт выхода однокубитового оператора
		for (uint i = (me / 2) * 2 + 2; i < _e.size(); i++)
		{
			if (!_busy[i])
			{
				_busy[i] = true;
				_e[me] = i;
				print_sifted_graphs(me + 1, _e, _busy, _sP, _fP, _g, _sM);
				_busy[i] = false;
			}
		}
	}
	else//В этой точке мы стартуем из порта ввода
    {
        if (me < (_e.size() - _sP))
        //Порт ввода может смотреть в любой узел графа
        for (uint i = 0; i < _e.size(); i++)
        {
            if (!_busy[i])
            {
                _busy[i] = true;
                _e[me] = i;
                print_sifted_graphs(me + 1, _e, _busy, _sP, _fP, _g, _sM);
                _busy[i] = false;
            }
        }
        else
        {
            
        }//Конец обработки сгенерированного направленного графа}
    }
}
