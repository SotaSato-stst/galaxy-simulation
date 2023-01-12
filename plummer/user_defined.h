typedef struct _vec
{
    double x, y, z;
} vec;

typedef struct full_particle
{
    double epot;       // potential enegy
    double ekin;       // kinetic enegy
    double etot;       // total enegy
    double vel_escape; // escape velocity
    vec pos;           // position
    vec vel;           // velocity
    vec acc;           // acceleration
} Full_particle;