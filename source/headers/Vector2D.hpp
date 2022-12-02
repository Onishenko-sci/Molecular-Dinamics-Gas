#ifndef THIS
#define THIS

#include <iostream>
#include <cmath>

class Vector2D
{
private:
public:
    double x;
    double y;
    double X() {return x;};
    void X(double X_in) {x= X_in;};
    double Y() {return y;};
    void Y(double Y_in) {y= Y_in;};

    Vector2D(double x_in,double y_in);
    Vector2D();
    friend std::ostream& operator<<(std::ostream& os, const Vector2D Vec);

    Vector2D operator*(double mult);
    friend Vector2D operator*( double mult, Vector2D Vec) { return Vec*mult; };

    Vector2D operator/(double del);
    Vector2D operator+(Vector2D AddVec) { return Vector2D(x+AddVec.x,y+AddVec.y);}
    Vector2D operator-(Vector2D AddVec) { return Vector2D(x-AddVec.x,y-AddVec.y);}
    bool operator==(Vector2D Vec) { if (x == Vec.x && y == Vec.y) {return 1;} return 0;}
    double abs() const { return sqrt(x*x+y*y);}
    void abs(double len);
    double fi() { return acos(x/(this->abs()));}

};

const Vector2D null_vec(0,0);

#endif