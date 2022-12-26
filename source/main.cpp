#include "./headers/MD_model.hpp"

//Практическое задание по молекулярной динамке.
//Задача №1 Двумерный газ 1
//Периодические граничные условия.
//Исследовать зависимость среднеквадратичного смещения и парной корреляционный фукнции от плотности газа и температуры.

int main()
{
    //Инициализация и задание параметров исследуемой молекулы.
    molecule Argon;
    Argon.mass = 6.63 * pow(10, -23);
    Argon.radius = 0.34 * pow(10, -9);
    Argon.pit_debt = 1.66 * pow(10, -21);
    Argon.interaction_radius = 2.5 * Argon.radius;

    
    molecular_dinamics gas(Argon, 100, 80, 40); // gas, number of particles, bound_x,bound_y in atom_radius
    gas.set_temperature(1.5);
    gas.simulate(1e+5, 1e-13);
    return 0;

}