#include"Vector2D.hpp"
#include<string>
#include<fstream>

struct molecule 
{
    double mass;
    double radius;
    double pit_debt;
    double interaction_radius = 2.5*radius;
};

struct particle
{
    Vector2D position{0, 0};
    Vector2D velosity{0, 0};
    Vector2D acceleration{0, 0};
    double potential;
};

const double K_b = 1.38 * pow(10, -23);

class molecular_dinamics
{
private:
    molecule mol;
    particle* particles;
    Vector2D *previous_step_acceleration;
    double Number_of_particles;
    double Bound_x;
    double Bound_y;

    void set_start_position();
    void reset_acceleration(int i);
    Vector2D radius_vec(int& i,int& j) const;
    void LJ_interact(int& i,int& j,const Vector2D& vector_between);
    void boundaries_check(int i);
    
    void write_header(std::ostream& file,int steps,double delta_t, int save_every_frame);
    void write_mol_info(std::ostream& file, int i);
    void write_step_info(std::ostream& file,int step);
    void write_correlation_data(std::ostream& file,int steps);

    double Temperature;
    double Potential;
    double Kin_Energy;
    double Square_displacment;
    Vector2D* start_posistions;
    Vector2D* displacement;

    double* Correlation;
    int correlation_points;
    double correlation_dr;
    void save_in_correlation();

public:
    molecular_dinamics(molecule& molecul, int number_of_particles, int bound_x = 100,int bound_y = 100);
    void set_temperature(double temperatur);
    void simulate(int steps,double delta_t, int save_every_frame = 50, std::string filename = "../render/md_render.txt");
    ~molecular_dinamics();
};

