double t_max = 50;
double x_max = 60;
double t_step = 0.01;
double x_step = 0.01;

double func(double t, double x){
    return x*t;
}

double fi(double x){
    return x*x*x / 12;
}

double ksi(double t){
    return t*t*t / 12;
}