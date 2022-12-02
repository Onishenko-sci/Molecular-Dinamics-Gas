#include "../headers/Vector2D.hpp"

void Vector2D::abs(double len)
{
    double angle = this->fi();
    x = len*cos(angle);
    y = len*sin(angle);
}

Vector2D::Vector2D(double x_in, double y_in)
{
    x = x_in;
    y = y_in;
};
Vector2D::Vector2D()
{
    x = 0;
    y = 0;
};

std::ostream &operator<<(std::ostream &os, const Vector2D Vec)
{
    os << "(" << Vec.x << "," << Vec.y << ")";
    return os;
};

Vector2D Vector2D::operator*(double mult)
{
    return Vector2D(mult * x, mult * y);
};

Vector2D Vector2D::operator/(double del)
{
    return Vector2D(x / del, y / del);
};