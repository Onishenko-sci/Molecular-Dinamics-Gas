#include "Vector2D.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <iomanip>
#include <omp.h>

#include <ctime>
#include <cstdlib>

struct molecula 
{
    Vector2D position{0, 0};
    Vector2D velosity{0, 0};
    Vector2D acceleration{0, 0};
    double potential;
};

const double K_b = 1.38 * pow(10, -23);

const double mol_radius = 0.34 * pow(10, -9);
const double mol_mass = 6.63 * pow(10, -23);
const double pit_depth = 1.66 * pow(10, -21);
const double interaction_radius = mol_radius * 2.5;

const double bound_x = 0.0000001;
const double bound_y = 0.0000001;

const double model_duration = 0.000000005;
const double delta_t = 0.00000000000001;
const int steps = model_duration / delta_t;
const int Save_every_frame = 100;

const int Number_of_particles = 50;
const double temperature =100.1;
const double average_velosity = sqrt(3 * K_b * temperature / mol_mass);
const double p_dencity = mol_mass*Number_of_particles/(bound_x*bound_y);

const int correlation_points = 50;
const double correlation_dr = bound_x*sqrt(2) / correlation_points;
const double correlation_radius = correlation_dr / 2;
const double correlation_constanta = (4*3.1415*p_dencity*Number_of_particles);

const double e4sig6 = 4 * pit_depth * pow(mol_radius, 6);
const double e4sig12 = 4 * pit_depth * pow(mol_radius, 12);
const double shift = (e4sig12 / pow(interaction_radius, 12)) - (e4sig6 / pow(interaction_radius, 6));


int binary_search(double a)
{
    int first = 0;
    int last = correlation_points;
    int b = 10;
    while(b<20)
    {
        b++;
        
        int mid = first + (last - first) / 2;
        double mid_value = mid* correlation_dr;
        if (a <= mid_value + correlation_radius && a >= mid_value - correlation_radius)
            return mid;
        if (a < mid_value)
            last = mid;
        else
            first = mid;
      //  std::cout << first << " " << last << " " << mid_value << " " << a << std::endl;
    }
    return first + (last - first) / 2;
}

int main()
{
    using namespace std;
    // srand(time(0));  //real random
    unsigned int start_time = clock();
    cout << "Average velosity: " << average_velosity << endl;
    cout << "Shift: " << shift << endl;
    cout << "Dencity: " << p_dencity << endl;

    ofstream File("render_out.txt");
    File << bound_x << ";"
         << bound_y << ";"
         << Number_of_particles << ";"
         << mol_radius << ";"
         << model_duration << ";"
         << delta_t << ';'
         << Save_every_frame << ';'
         << correlation_points << "\n";
    File << setprecision(4);
    // inicialization

    molecula Molecules[Number_of_particles];
    Vector2D forces_before[Number_of_particles];

    random_device rd{};
    mt19937 gen{rd()};
    normal_distribution<> d{average_velosity, average_velosity / 10};

    Vector2D random_position;
    bool far_enough;
    Vector2D start_positions[Number_of_particles];
    Vector2D displacment[Number_of_particles];
    double square_displacment;

    for (int i = 0; i < Number_of_particles; i++) // init positions and speeds
    {
        while (Molecules[i].position.X() == 0.0 || Molecules[i].position.Y() == 0.0)
        {
            random_position.X((double)bound_x * rand() / RAND_MAX);
            random_position.Y((double)bound_y * rand() / RAND_MAX);
            far_enough = true;
            for (int j = 0; j < i; j++)
            {
                if ((Molecules[j].position - random_position).abs() < 2 * mol_radius)

                    far_enough = false;
            }
            if (far_enough)
                Molecules[i].position = random_position;
        }
        start_positions[i] = Molecules[i].position;
        displacment[i] = start_positions[i];

        Molecules[i].velosity.X(rand() - RAND_MAX / 2);
        Molecules[i].velosity.Y(rand() - RAND_MAX / 2);
        Molecules[i].velosity.abs(d(gen));
    }

    double correlation[correlation_points];
    int current_step = 0;
    cout << "Calculating:" << endl
         << "1/10: [";
    int counter = 1;
    for (double t = 0; t < model_duration; t += delta_t) // Main loop
    {

        current_step += 1;
        if (current_step % (steps / 1000) == 0)
        {
            cout << /*current_step / (steps / 100) <<*/ "-" << flush;
            if ((current_step / (steps / 1000)) % 100 == 0)
                cout << "]" << endl
                     << ++counter
                     << "/10: [";
        }

        for (int i = 0; i < Number_of_particles; i++) // Calculating Position for next step
        {
            forces_before[i] = Molecules[i].acceleration;
            Molecules[i].acceleration = null_vec;
            Molecules[i].position = Molecules[i].position + Molecules[i].velosity * delta_t + (forces_before[i] * pow(delta_t, 2)) / 2;
            Molecules[i].potential = 0;
        }

        Vector2D delta_position;
        Vector2D nearest_copy;
        Vector2D force;
        for (int i = 0; i < Number_of_particles; i++) // Calculating Force for next step
        {
            for (int j = i + 1; j < Number_of_particles; j++)
            {
                delta_position = Molecules[j].position - Molecules[i].position;
                correlation[binary_search(delta_position.abs())]+=2;
                { // Finding nearest copy in tor
                    nearest_copy = Molecules[j].position;
                    if (delta_position.x > bound_x / 2)
                        nearest_copy.X(nearest_copy.x - bound_x);
                    else if (delta_position.x <= -bound_x / 2)
                        nearest_copy.X(nearest_copy.x + bound_x);

                    if (delta_position.y > bound_y / 2)
                        nearest_copy.Y(nearest_copy.y - bound_y);
                    else if (delta_position.y <= -bound_y / 2)
                        nearest_copy.Y(nearest_copy.y + bound_y);
                }

                double PotentialLJ = 0;
                force = null_vec;
                delta_position = nearest_copy - Molecules[i].position;
                double distance = delta_position.abs();

                if (distance < interaction_radius)
                {
                    PotentialLJ = (e4sig12 / pow(distance, 12)) - (e4sig6 / pow(distance, 6)) - shift;
                    Molecules[i].potential = Molecules[i].potential + PotentialLJ;
                    Molecules[j].potential = Molecules[j].potential + PotentialLJ;

                    force = 12 * ((e4sig12 / pow(distance, 14)) - (0.5 * e4sig6 / pow(distance, 8))) * delta_position;
                    Molecules[i].acceleration = Molecules[i].acceleration + force / mol_mass;
                    Molecules[j].acceleration = Molecules[j].acceleration - force / mol_mass;
                }
            }
        }

        for (int i = 0; i < Number_of_particles; i++) // Calculating Velosity for next step
        {
            Molecules[i].velosity = Molecules[i].velosity + ((Molecules[i].acceleration + forces_before[i]) / 2) * delta_t;
        }

        for (int i = 0; i < Number_of_particles; i++) // teleport on other bound
        {
            if (Molecules[i].position.X() > bound_x)
            {
                Molecules[i].position.X(0);
            }
            else if (Molecules[i].position.X() < 0)
            {
                Molecules[i].position.X(bound_x);
            }

            if (Molecules[i].position.Y() > bound_y)
            {
                Molecules[i].position.Y(0);
            }
            else if (Molecules[i].position.Y() < 0)
            {
                Molecules[i].position.Y(bound_y);
            }
        }

        if (current_step % Save_every_frame == 0) // Out in File
        {
            double cinetic_energy = 0;
            double potential_energy = 0;
            square_displacment =0;
            for (int i = 0; i < Number_of_particles; i++)
            {
                cinetic_energy += mol_mass * pow(Molecules[i].velosity.abs(), 2) / 2;
                potential_energy += Molecules[i].potential;
                displacment[i] = displacment[i] + Molecules[i].velosity*delta_t;
                square_displacment+= pow((displacment[i]-start_positions[i]).abs(),2);
                File << Molecules[i].position.X() << ";"
                     << Molecules[i].position.Y() << ";"
                     << Molecules[i].potential << ";"
                     << Molecules[i].velosity.abs() << ";"
                     << Molecules[i].acceleration.abs() << "\n";
            }
            double current_temperature = (cinetic_energy * 2) / (3 * K_b * Number_of_particles);
            File << "///----End----///;"
                 << current_step << ";"
                 << cinetic_energy << ';'
                 << potential_energy << ';'
                 << square_displacment / Number_of_particles << ';'
                 << current_temperature << ';' << endl;
        }
    }
    File << "1488;" << endl;

    for (int i = 0; i < correlation_points; i++)
    {
        File << i*correlation_dr << ';' << (correlation[i]/steps)/(correlation_constanta*pow(i*correlation_dr,2)) << ';' << endl;
    }
    

    cout << "Время выполнения программы: " << (float)(clock() - start_time) / CLOCKS_PER_SEC << endl;
    File.close();
    return 0;
}