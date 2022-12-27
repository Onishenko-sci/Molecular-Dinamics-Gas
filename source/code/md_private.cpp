#include "../headers/Vector2D.hpp"
#include "../headers/MD_model.hpp"
#include <random>
#include <iostream>

void molecular_dinamics::reset_acceleration(int i)
{
    previous_step_acceleration[i] = particles[i].acceleration;
    particles[i].acceleration = null_vec;
    particles[i].potential = 0;
}

void molecular_dinamics::set_start_position()
{
    // Генератор псевдослучайных величин. Равноменое распределение.
    std::random_device rd;
    std::uniform_real_distribution<double> random_x(0.0, Bound_x - mol.interaction_radius);
    std::uniform_real_distribution<double> random_y(0.0, Bound_y - mol.interaction_radius);

    // Метод генерирует случайное положение для частицы, и проверяет его. Если все остальные частицы находятся
    // находится достаточно далеко (т.е. дальше радиуса взаимодействия), случайное положение принимается.
    // Если новая частица находится слишком близко к другой частице, выбирается новое случайное пложение.
    // Так происходит пока у всех частиц будет ненулевая координата.
    Vector2D random_position;
    for (int i = 0; i < Number_of_particles; i++)
    {
        while (particles[i].position == null_vec)
        {
            bool far_enough = true;
            random_position.x = random_x(rd);
            random_position.y = random_y(rd);
            for (int j = 0; j < i; j++)
            {
                if ((particles[j].position - random_position).abs() <= mol.radius*1.5)
                    far_enough = false;
            }
            if (far_enough)
                particles[i].position = random_position;
        }
        start_posistions[i] = particles[i].position;
        displacement[i] = start_posistions[i];
    }
        std::cout << "Molecules positioned" << std::endl;
}

Vector2D molecular_dinamics::radius_vector(int &i, int &j) const
{

    Vector2D nearest_copy = particles[j].position;

    // Вектор между частицами без учета периодических условий
    Vector2D delta_position = particles[j].position - particles[i].position;

    // Учет граничных условий по оси Х
    if (delta_position.x >= Bound_x / 2)
        nearest_copy.x = nearest_copy.x - Bound_x;
    else if (delta_position.x <= -Bound_x / 2)
        nearest_copy.x = nearest_copy.x + Bound_x;

    // Учет граничных условий по оси Y
    if (delta_position.y >= Bound_y / 2)
        nearest_copy.y = nearest_copy.y - Bound_y;
    else if (delta_position.y <= -Bound_y / 2)
        nearest_copy.y = nearest_copy.y + Bound_y;

    // Итоговый вектор
    delta_position = nearest_copy - particles[i].position;
    return delta_position;
}

void molecular_dinamics::save_in_correlation()
{
    // По всем парам частиц
    for (int i = 0; i <= Number_of_particles - 1; i++)
        for (int j = i + 1; j < Number_of_particles; j++)
        {
            // Находим расстояние между частицами.
            double distance = radius_vector(i, j).abs();
            // Подсчитываем только частицы, попавшие в диапазон [0,correlation_max_r]
            if (distance < correlation_max_r)
                Correlation[int(distance / correlation_dr)] += 2;
        }
}

void molecular_dinamics::LJ_interact(int &i, int &j)
{
    // Расчет констант
    static const double e4sig6 = 4 * mol.pit_debt * pow(mol.radius, 6);
    static const double e4sig12 = 4 * mol.pit_debt * pow(mol.radius, 12);
    static const double shift = (e4sig12 / pow(mol.interaction_radius, 12)) - (e4sig6 / pow(mol.interaction_radius, 6));

    // Расчет радиус вектора между i-й частицой и ближайшей копией (из граничных условий) j-й частицы.
    Vector2D radius_vec = radius_vector(i, j);
    double distance = radius_vec.abs();

    // Расчет потенциалов и сил.
    Vector2D force;
    double PotentialLJ;
    if (distance < mol.interaction_radius)
    {
        PotentialLJ = (e4sig12 / pow(distance, 12)) - (e4sig6 / pow(distance, 6)) - shift;
        // Учет взаимодействия между частицами в суммарном потенциале чистицы.
        particles[i].potential = particles[i].potential + PotentialLJ;
        particles[j].potential = particles[j].potential + PotentialLJ;

        force = -12 * ((e4sig12 / pow(distance, 14)) - (0.5 * e4sig6 / pow(distance, 8))) * radius_vec;
        // Учет взаимодействия между частицами в суммарной силе, действующей на частицу.
        particles[i].acceleration = particles[i].acceleration + force / mol.mass;
        particles[j].acceleration = particles[j].acceleration - force / mol.mass;
    }
}

void molecular_dinamics::boundaries_check(int i)
{
    // При выходе частицы за пределы участка c одной стороны, частица переносится на другую сторону.

    if (particles[i].position.x > Bound_x)
        particles[i].position.x = particles[i].position.x - Bound_x;
    else if (particles[i].position.X() < 0)
        particles[i].position.x = particles[i].position.x + Bound_x;

    if (particles[i].position.y > Bound_y)
        particles[i].position.y = particles[i].position.y - Bound_y;
    else if (particles[i].position.y < 0)
        particles[i].position.y = particles[i].position.y + Bound_y;
}

void molecular_dinamics::write_header(int frame_rate)
{
    // Используется формат данных csv.
    File << Bound_x << ";"
         << Bound_y << ";"
         << Number_of_particles << ";"
         << mol.radius << ";"
         << number_of_steps << ";"
         << time_step << ';'
         << frame_rate << ';'
         << correlation_points
         << "\n";
}

void molecular_dinamics::write_mol_info(int i)
{
    File << particles[i].position.x << ";"
         << particles[i].position.y << ";"
         << particles[i].potential << ";"
         << particles[i].velosity.abs() << ";"
         << particles[i].acceleration.x << "\n";
}

void molecular_dinamics::write_observables(int current_step)
{
    File << "///----End----///;"
         << current_step << ";"
         << Kin_Energy << ';'
         << Potential << ';'
         << Square_displacment / Number_of_particles << ';'
         << Temperature << ';' << std::endl;
}

void molecular_dinamics::write_correlation_data()
{
    File << "101; Correlation data section" << std::endl;
    // Расчет констант
    const double volume = Bound_x * Bound_y;
    const double normalization = (number_of_steps / 2) * (Number_of_particles - 1) * Number_of_particles;
    // Расчет корреляционный функции
    for (int i = 1; i < correlation_points; i++)
        File << i * correlation_dr << ';'
             << ((volume / normalization) * Correlation[i]) / (3.1415926535 * (pow(i * correlation_dr + correlation_dr, 2) - pow(i * correlation_dr, 2))) << ';'
             << std::endl;
}
