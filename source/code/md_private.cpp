#include "../headers/Vector2D.hpp"
#include "../headers/MD_model.hpp"
#include <random>
#include <iostream>
#include <omp.h>

void molecular_dinamics::reset_acceleration(int i)
{
    previous_step_acceleration[i] = particles[i].acceleration;
    particles[i].acceleration = null_vec;
    particles[i].potential = 0;
}

void molecular_dinamics::set_start_position()
{
    Vector2D random_position;
    std::random_device rd;
    std::uniform_real_distribution<double> random_x(0.0, Bound_x);
    std::uniform_real_distribution<double> random_y(0.0, Bound_y);
    for (int i = 0; i < Number_of_particles; i++) // init positions and speeds
    {
        while (particles[i].position == null_vec)
        {
            bool far_enough = true;
            random_position.x = random_x(rd);
            random_position.y = random_y(rd);
            for (int j = 0; j < i; j++)
            {
                if ((particles[j].position - random_position).abs() < mol.interaction_radius)
                    far_enough = false;
            }
            if (far_enough)
                particles[i].position = random_position;
        }
        start_posistions[i] = particles[i].position;
        displacement[i] = start_posistions[i];
    }
}

Vector2D molecular_dinamics::radius_vec(int &i, int &j) const
{
    Vector2D nearest_copy = particles[j].position;
    Vector2D delta_position = particles[j].position - particles[i].position;
    if (delta_position.x >= Bound_x / 2)
        nearest_copy.X(nearest_copy.x - Bound_x);
    else if (delta_position.x <= -Bound_x / 2)
        nearest_copy.X(nearest_copy.x + Bound_x);

    if (delta_position.y >= Bound_y / 2)
        nearest_copy.Y(nearest_copy.y - Bound_y);
    else if (delta_position.y <= -Bound_y / 2)
        nearest_copy.Y(nearest_copy.y + Bound_y);

    delta_position = nearest_copy - particles[i].position;
    return delta_position;
}

void molecular_dinamics::save_in_correlation()
{
    for (int i = 0; i < Number_of_particles - 1; i++)
        for (int j = i + 1; j < Number_of_particles; j++)
        {
            double distance = radius_vec(i, j).abs();
            if (distance < Bound_x / 2)
                Correlation[int(distance / correlation_dr)] += 2;
        }
}

void molecular_dinamics::LJ_interact(int &i, int &j, const Vector2D &radius_vec)
{
    static const double e4sig6 = 4 * mol.pit_debt * pow(mol.radius, 6);
    static const double e4sig12 = 4 * mol.pit_debt * pow(mol.radius, 12);
    static const double shift = (e4sig12 / pow(mol.interaction_radius, 12)) - (e4sig6 / pow(mol.interaction_radius, 6));
    double PotentialLJ;
    Vector2D force;
    double distance = radius_vec.abs();
    if (distance < mol.interaction_radius)
    {
        PotentialLJ = (e4sig12 / pow(distance, 12)) - (e4sig6 / pow(distance, 6)) - shift;
        particles[i].potential = particles[i].potential + PotentialLJ;
        particles[j].potential = particles[j].potential + PotentialLJ;

        force = 12 * ((e4sig12 / pow(distance, 14)) - (0.5 * e4sig6 / pow(distance, 8))) * radius_vec;
        particles[i].acceleration = particles[i].acceleration + force / mol.mass;
        particles[j].acceleration = particles[j].acceleration - force / mol.mass;
    }
}

void molecular_dinamics::boundaries_check(int i)
{
    if (particles[i].position.X() > Bound_x)
        particles[i].position.X(0);
    else if (particles[i].position.X() < 0)
        particles[i].position.X(Bound_x);

    if (particles[i].position.Y() > Bound_y)
        particles[i].position.Y(0);
    else if (particles[i].position.Y() < 0)
        particles[i].position.Y(Bound_y);
}

void molecular_dinamics::write_header(std::ostream &file, int steps, double delta_t, int save_every_frame)
{
    file << Bound_x << ";"
         << Bound_y << ";"
         << Number_of_particles << ";"
         << mol.radius << ";"
         << steps << ";"
         << delta_t << ';'
         << save_every_frame << ';' // save every frame
         << correlation_points
         << "\n"; // correlation points
}

void molecular_dinamics::write_mol_info(std::ostream &file, int i)
{
    file << particles[i].position.x << ";"
         << particles[i].position.y << ";"
         << particles[i].potential << ";"
         << particles[i].velosity.abs() << ";"
         << particles[i].acceleration.abs() << "\n";
}

void molecular_dinamics::write_step_info(std::ostream &file, int current_step)
{
    file << "///----End----///;"
         << current_step << ";"
         << Kin_Energy << ';'
         << Potential << ';'
         << Square_displacment / Number_of_particles << ';'
         << Temperature << ';' << std::endl;
}

void molecular_dinamics::write_correlation_data(std::ostream &file, int steps)
{
    file << "1488;" << std::endl;
    const double p_dencity = mol.mass * Number_of_particles / (Bound_x * Bound_y);
    const double correlation_constanta = (4 * 3.1415 * p_dencity * Number_of_particles);
    const double normalization = (steps / 2) * (Number_of_particles - 1) * Number_of_particles / 2;
    for (int i = 1; i < correlation_points; i++)
        file
            << i * correlation_dr << ';'
            //   << Correlation[i] / normalization << ';'
            << Correlation[i] / (normalization * correlation_constanta * pow(i * correlation_dr, 2)) << ';'
            << std::endl;
}

