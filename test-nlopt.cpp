#include <nlopt.hpp>
#include <vector>

#define DIM 2

/*
 * @brief Реализация NLopt. 
 *
 * @param g         направленный граф
 * @param L         логическая матрица из функций
 * @param restrictions массив условий равенства
 * @param eps       точность удовлетворения условиям равенства
 *
 * @return максимальное значение целевой функции
 */
float NLopt (std::vector<float> &x, float eps)
{
    //! Поиск глобального оптимума, без производных
    nlopt::opt glob_problem(nlopt::AUGLAG, DIM);
    
    glob_problem.set_min_objective(target_function);

    //! Устанавливаем границы изменения переменных
    vector<float> lb(DIM, -1), ub(DIM, 1);
    glob_problem.set_lower_bounds(lb);
    glob_problem.set_upper_bounds(ub);

    //Задаём конечную точность установления переменных
    glob_problem.set_xtol_abs(eps);

    {
        //! Локальный оптимизатор
        nlopt::opt loc_problem(nlopt::LN_COBYLA, DIM);
        loc_problem.set_xtol_abs(eps);
        //ПРЕЖДЕ локальный оптимизатор надо конфигурировать ДО того как 
        //передать его для _копирования_ глобальному. 
        glob_problem.set_local_optimizer(loc_problem);
    }

    float result;
    std::vector<float> grad;

    result = glob_problem.optimize(x, grad, NULL);

    return result;
}

float target_function(const std::vector<float> &x, std::vector<float> &grad, void * data)
{
    float sum = 0;

    for(auto i : x)    
        sum += i*i;
    
    return sum;
}

int main(void)
{


    return 0;
}