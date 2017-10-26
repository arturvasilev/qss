#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>
#include <set>
#include <complex>
#include <string>
#include <stdlib.h>

class graph {
public:
	
	//! Тип комплексной 2D-матрицы
	typedef std::vector<std::vector<std::complex<double> > > cmatrix_t;

	//! Тип матрицы траекторий
	typedef std::vector<std::vector<std::set<std::set<uint> > > > traj_t;

	//! Тип матрицы для просеивания
	typedef std::vector<std::vector<bool> > smatrix_t;
	
	//! Поддерживаемые типы однокубитовых элементов
	enum operators_types
    {
        beamsplitter,
        directCoupler,
        waveplate
	};
	/*
     * @brief Конструктор по умолчанию
     *
	 * @param p     число портов ввода-вывода
	 * @param cp    число светоделительных пластинок
	 * @param dc    число направленных светоделителей
	 * @param w     число волновых пластинок (фазовращателей)
	 */
	graph(uint ports = 0, uint beamsplitters = 0, uint directCouplers = 0, uint waveplates = 0);

	//! Установить рёбра графа
	void set_edges(const std::vector<uint> _edges);

	/* 
	 * Вернуть эффективность текущего графа относительно целевой матрицы
	 * 
	 * @param _tM		Целевая матрица
	 * 
	 * @return Текущее значение эффективности
	 */
	double get_deviation();
	static double get_deviation(
		const std::vector<double> &x, 
		std::vector<double> &grad, 
		void* f_data);
	
	/*
	 * @brief Устанавливает целевую матрицу истинности
	 * 
	 * @param _tM		Целевая матрица
	 */
	void set_target_matrix(const cmatrix_t &_tM);
	
	//! Возвращает строку с текстовым представлением графа для ввода-вывода
	const std::string print(void);

	/*
	 * @brief Устанавливает внутренние параметры графа
	 * 
	 * @param _var		Новые внутренние параметры
	 */
	void set_variables(const std::vector<double> &_var);

	/*
	 * @brief Функция для просеивания графа на основе матрицы траекторий.
	 * 	Все ненулевые элементы матрицы истинности должны иметь хоть одну траекторию.
	 * 
	 * @param _sM		Булевая целевая матрица истинности
	 * 
	 * @return true в случае просеянного графа. Иначе false.
	 */
	bool sift(const smatrix_t &_sM);

protected:

	//! Направленный граф
	std::vector<uint> edges;
	
	// //! Хранит в себе инфу о том, в какие из вершин ещё не смотрит ребро
	// std::set<uint> busy;
	
	//! Хранит в себе комбинацию типов однокубитовых операторов
	std::vector<operators_types> comb;

	//! Матрица траекторий
	traj_t traj;

	uint p;     //!< p(orts) - число портов ввода-вывода
	uint bs;    //!< (beamsplitters) - число светоделительных пластинок
	uint dc;    //!< (direction couplers) - число направленных светоделителей
	uint w;     //!< (waveplates) - число волновых пластинок (фазовращателей)
	

	std::vector<double> var;//Внутренние параметры графа

	//! Целевая матрица
	cmatrix_t targetMatrix;

	//! Матрица конвертации матрицы амплитуд (или траекторий) в матрицу истинности
	std::vector<std::vector<std::vector<uint> > > translate;

	//! Создаёт матрицу траекторий
	void make_matrix_traj();

    /*
     * @brief Возвращает указатель на переменную из массива var[], соответствующей
     *  однокубитовому оператору comb->op[oper_num]. Если для данного
     *  однокубитового оператора нужно несколько переменных - то будет
     *  возвращёна ссылка на первую из них.
     *
     * @param oper_num Номер оператора
     *
     * @return Указатель на переменную из массива var[]
	 */ 
	double* var_num(uint oper_num)
	{
		//Вычислим номер необходимого оператора из массива var[]
		uint var_num = 0;
		for (uint i = 0; i < oper_num; i++)
			switch (comb[i])
			{
			case beamsplitter:
			case directCoupler: var_num++; break;
			case waveplate: var_num += 2; break;
			};
		return &(var[var_num]);
	};

	/*
	 * @param oper_num		Порядковый номер оператора (находится из рёбер графа)
	 * @param in			Порт ввода
	 * @param out			Порт вывода
	 * 
	 * @return Комплексная амплитуда однокубитового оператора
	 */
	std::complex<double> get_func(uint oper_num, uint in, uint out);
	
	//Определяет какому типу однокубитового оператора принадлежит переменная с номером var_num
	operators_types oper_type(uint var_num);

	graph& operator= (const graph &other);

	//! Рекурсивно находит все возможные пути
	static void paths(
		uint start, 
		graph &g, 
		traj_t &matrix, 
		std::vector<uint> way = std::vector<uint>());	

	/*
	 * @brief Конвертирует тректорию в комплексную амплитуду
	 */
	std::complex<double> traj_to_ampl(const std::set<uint> &_traj);
};

#endif //! GRAPH_HPP