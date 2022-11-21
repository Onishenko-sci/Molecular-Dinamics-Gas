#include "../headers/Vector2D.hpp"
#include "../headers/MD_model.hpp"
#include <random>

molecular_dinamics::molecular_dinamics(molecule &molecul, int number_of_particles, int bound_x, int bound_y)
{
    particles = new particle[number_of_particles];
    start_posistions = new Vector2D[number_of_particles];
    displacement = new Vector2D[number_of_particles];
    Number_of_particles = number_of_particles;
    Bound_x = mol.radius * bound_x;
    Bound_y = mol.radius * bound_y;

    set_start_position();
}

molecular_dinamics::~molecular_dinamics()
{
    delete[] particles;
    delete[] start_posistions;
    delete[] displacement;
}

void molecular_dinamics::reset_acceleration(Vector2D *previous, int i)
{
    previous[i] = particles[i].acceleration;
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
                if ((particles[j].position - random_position).abs() < mol.radius)
                    far_enough = false;
            }
            if (far_enough)
                particles[i].position = random_position;
        }
        start_posistions[i] = particles[i].position;
    }
}

void molecular_dinamics::set_temperature(double temperatur)
{
    Temperature = temperatur;
    const double average_velosity = sqrt(3 * K_b * Temperature / mol.mass);

    std::random_device rd;
    std::mt19937 gen{rd};
    std::normal_distribution<double> normal_dist{average_velosity, average_velosity / 5};
    std::uniform_real_distribution<double> random(-1.0, 1.0);

    for (int i = 0; i < Number_of_particles; i++)
    {
        particles[i].velosity.x = random(rd);
        particles[i].velosity.y = random(rd);
        particles[i].velosity.abs(normal_dist(gen));
    }
}

Vector2D molecular_dinamics::radius_vec(int i, int j)
{
    Vector2D nearest_copy = particles[j].position;
    Vector2D delta_position = particles[j].position - particles[i].position;
    if (delta_position.x > Bound_x / 2)
        nearest_copy.X(nearest_copy.x - Bound_x);
    else if (delta_position.x <= -Bound_x / 2)
        nearest_copy.X(nearest_copy.x + Bound_x);

    if (delta_position.y > Bound_y / 2)
        nearest_copy.Y(nearest_copy.y - Bound_y);
    else if (delta_position.y <= -Bound_y / 2)
        nearest_copy.Y(nearest_copy.y + Bound_y);

    delta_position = nearest_copy - particles[i].position;
    return delta_position;
}

void molecular_dinamics::LJ_interact(int i, int j, Vector2D radius_vec)
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

void molecular_dinamics::simulate(int steps, double delta_t)
{

    Vector2D *previous_step_acceleration = new Vector2D[Number_of_particles];
    for (int step = 0; step < steps; step++)
    {
        for (int i = 0; i < Number_of_particles; i++)
        {
            particles[i].position = particles[i].position +
                                    particles[i].velosity * delta_t +
                                    (particles[i].acceleration * pow(delta_t, 2)) / 2;
            boundaries_check(i);
            reset_acceleration(previous_step_acceleration, i);
        }

        Vector2D vector_between;
        for (int i = 0; i < Number_of_particles; i++)
        {
            for (int j = i + 1; j < Number_of_particles; j++)
            {
                vector_between = radius_vec(i, j);
                LJ_interact(i, j, vector_between);
            }
            particles[i].velosity = particles[i].velosity +
                                    ((particles[i].acceleration + previous_step_acceleration[i]) / 2) * delta_t;
        }
        
        find();
        step_out();
        step_info();
    }
}