#include "./headers/MD_model.hpp"

int main()
{
    molecule Argon;
    Argon.mass = 6.63 * pow(10, -23);
    Argon.radius =  0.34 * pow(10, -9);
    Argon.pit_debt = 1.66 * pow(10, -21);
    Argon.interaction_radius = 2.5*Argon.radius;

    molecular_dinamics gas(Argon,100,100,100);
    gas.set_temperature(10.5);
    gas.simulate(1e+5,1e-13,500);
    
    return 0;
}