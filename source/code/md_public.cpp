#include "../headers/Vector2D.hpp"
#include "../headers/MD_model.hpp"
#include <random>
#include <iostream>

molecular_dinamics::molecular_dinamics(molecule &molecul, int number_of_particles, int bound_x, int bound_y)
{
    // Инициализация переменных модели
    Number_of_particles = number_of_particles;
    mol = molecul;
    Bound_x = mol.radius * bound_x;
    Bound_y = mol.radius * bound_y;

    // Выделение памяти для массивов частиц
    particles = new particle[number_of_particles];
    start_posistions = new Vector2D[number_of_particles];
    displacement = new Vector2D[number_of_particles];
    previous_step_acceleration = new Vector2D[number_of_particles];

    // Параметры для вычисления корреляционной функции
    correlation_max_r = mol.radius * 5;                      // Максимальная дистанция, для которой будет считаться кор.функция
    correlation_points = 100;                                // Количество точек на графике корреляционной функции
    correlation_dr = correlation_max_r / correlation_points; // Шаг корреляционной функции
    Correlation = new double[correlation_points];

    // Метод, распределяющий атомы на участке (подробное описание внутри функции)
    set_start_position();

    std::cout << "Density is " << double(Number_of_particles) / (bound_x * bound_y) << std::endl;
}

molecular_dinamics::~molecular_dinamics()
{
    delete[] particles;
    delete[] start_posistions;
    delete[] displacement;
    delete[] previous_step_acceleration;
}

void molecular_dinamics::set_temperature(double temperatur)
{
    Temperature = temperatur;
    Target_temperature = temperatur;
    // Расчет средней скорости молекулы.
    const double average_velosity = sqrt(2 * K_b * Target_temperature / mol.mass);
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::normal_distribution<double> normal_dist{average_velosity, average_velosity / 10}; // Нормальное распределение
    std::uniform_real_distribution<double> random(-1.0, 1.0);                              // Равномерное распределение

    // Метод для каждой частицы выбирает направление вектора скорости.
    // Затем устанавливает абсолютное значение вектора скорости, полученное из нормального распределения.
    Vector2D sum_velosity = null_vec;
    for (int i = 0; i < Number_of_particles; i++)
    {
        particles[i].velosity.x = random(rd);
        particles[i].velosity.y = random(rd);
        particles[i].velosity.abs(normal_dist(gen));
        sum_velosity = sum_velosity + particles[i].velosity;
    }

    sum_velosity = sum_velosity / Number_of_particles;
    for (int i = 0; i < Number_of_particles; i++)
        particles[i].velosity = particles[i].velosity - sum_velosity;
    
}

void molecular_dinamics::simulate(int steps, double delta_t, int frame_rate, std::string filename)
{
    File.open(filename);
    number_of_steps = steps;
    time_step = delta_t;
    double lyamda = 1;

    // Метод записывает в файл заголовок с параметрами модели.
    // Частота сохранения кадров передаётся напрямую, т.к. не схраняется в классе.
    write_header(frame_rate);

    // Основной цикл моделирования.
    for (int step = 0; step < steps; step++)
    {

        // Индикатор выполнения программы (progress bar)
        if (step % (steps / 100) == 0)
            std::cout << "-" << std::flush;

        // Первый шаг интегратора.
        // 1.Вычисление координаты частицы, проверка граничных условий и обнуление сил/ускорений.
        Kin_Energy = 0;
        for (int i = 0; i < Number_of_particles; i++)
        {
            particles[i].position = particles[i].position + particles[i].velosity * delta_t + (particles[i].acceleration * pow(delta_t, 2)) / 2;
            // Проверка граничных условий и перенос частицы.
            boundaries_check(i);
            // Зануление ускорений частицы.
            reset_acceleration(i);
            Kin_Energy += mol.mass * pow(particles[i].velosity.abs(), 2) / 2;
        }

        // Macштабирование скоростей
        Temperature = (Kin_Energy) / (K_b * Number_of_particles);
        lyamda = std::sqrt(Target_temperature / Temperature);

        // Второй и третий шаг интегратора.
        // 2.Расчет сил, действующих на частицы из потенциала Леннарда-Джонса.
        // 3.Расчет скорости частиц.
        for (int i = 0; i <= Number_of_particles - 1; i++)
        {
            // Все пары частиц.
            for (int j = i + 1; j < Number_of_particles; j++)
                LJ_interact(i, j);
            // Метод принимает номер двух частиц.
            // Рассчитывает потенциал Леннарда-Джонса и силы, действующие на частицы. Применяет их к обоим частицам.

            // 3.Расчет скорости частиц.
            Vector2D Half_step_acceleration = ((particles[i].acceleration + previous_step_acceleration[i]) / 2) * delta_t;
            particles[i].velosity = lyamda * particles[i].velosity + Half_step_acceleration;
            // Расчет наблюдаемых: отклонение частицы от начального положения.
            displacement[i] = displacement[i] + particles[i].velosity * delta_t + Half_step_acceleration * delta_t / 2;
        }

        // Метод учитывает кадр при постройке гистрограммы количества частиц от расстояния между ними
        if (step > steps / 2)
            save_in_correlation();

        // Расчет наблюдаемых и вывод кадра в файл.
        if (step % frame_rate == 0)
        {
            Potential = 0;
            Square_displacment = 0;
            for (int i = 0; i < Number_of_particles; i++)
            {
                // Расчет наблюдаемых.
                Potential += particles[i].potential;
                Square_displacment += pow((displacement[i] - start_posistions[i]).abs(), 2);
                // Вывод информации о частице в файл.
                write_mol_info(i);
            }
            // Расчет температуры и вывод наблюдаемых в файл.
            Temperature = (Kin_Energy) / (K_b * Number_of_particles);
            write_observables(step);
        }
    }

    // Метод рассчитывает кореляционную функцию и выводит её в файл.
    write_correlation_data();

    File.close();
    std::cout << "\nDone! \n";
}