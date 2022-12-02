#include "../headers/Vector2D.hpp"
#include "../headers/MD_model.hpp"
#include <random>
#include <iostream>
#include <omp.h>

molecular_dinamics::molecular_dinamics(molecule &molecul, int number_of_particles, int bound_x, int bound_y)
{
    Number_of_particles = number_of_particles;
    particles = new particle[number_of_particles];
    start_posistions = new Vector2D[number_of_particles];
    displacement = new Vector2D[number_of_particles];
    previous_step_acceleration = new Vector2D[number_of_particles];
    correlation_points = 100;
    Correlation = new double[correlation_points];
    mol = molecul;
    Bound_x = mol.radius * bound_x;
    Bound_y = mol.radius * bound_y;
    correlation_dr = (Bound_x / 2) / correlation_points;
    set_start_position();
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
    const double K_b = 1.38 * pow(10, -23);
    const double average_velosity = sqrt(3 * K_b * Temperature / mol.mass);

    std::random_device rd;
    std::mt19937 gen{rd()};
    std::normal_distribution<double> normal_dist{average_velosity, average_velosity / 5};
    std::uniform_real_distribution<double> random(-1.0, 1.0);

    for (int i = 0; i < Number_of_particles; i++)
    {
        particles[i].velosity.x = random(rd);
        particles[i].velosity.y = random(rd);
        particles[i].velosity.abs(normal_dist(gen));
    }
}

void molecular_dinamics::simulate(int steps, double delta_t, int save_every_frame, std::string filename)
{
    std::ofstream File(filename);

    write_header(File, steps, delta_t, save_every_frame);

    for (int step = 0; step < steps; step++)
    {
        if (step % (steps / 100) == 0)
            std::cout << "-" << std::flush;

        for (int i = 0; i < Number_of_particles; i++)
        {
            particles[i].position = particles[i].position + particles[i].velosity * delta_t + (particles[i].acceleration * pow(delta_t, 2)) / 2;
            boundaries_check(i);
            reset_acceleration(i);
        }

        for (int i = 0; i < Number_of_particles - 1; i++)
        {
            for (int j = i + 1; j < Number_of_particles; j++)
                LJ_interact(i, j, radius_vec(i, j));
            particles[i].velosity = particles[i].velosity +
                                    ((particles[i].acceleration + previous_step_acceleration[i]) / 2) * delta_t;
            displacement[i] = displacement[i] + particles[i].velosity * delta_t;
        }

        if (step >= steps / 2)
            save_in_correlation();

        if (step % save_every_frame == 0)
        {
            Kin_Energy = 0;
            Potential = 0;
            Square_displacment = 0;
            for (int i = 0; i < Number_of_particles; i++)
            {
                Kin_Energy += mol.mass * pow(particles[i].velosity.abs(), 2) / 2;
                Potential += particles[i].potential;
                Square_displacment += pow((displacement[i] - start_posistions[i]).abs(), 2);
                write_mol_info(File, i);
            }
            Temperature = (Kin_Energy * 2) / (3 * K_b * Number_of_particles);
            write_step_info(File, step);
        }
    }

    write_correlation_data(File, steps);

    File.close();
    std::cout << "\n Done! \n";
}