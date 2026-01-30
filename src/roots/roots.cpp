#include "roots.hpp"

bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root) {
    const double tol = 1e-7;
    const int max_iter = 1000;            
    
    if (f(a) * f(b) >= 0) {
        return false; // No sign change
    }      

    for (int i = 0; i < max_iter; ++i) {
        double c = (a + b) / 2;
        if (std::abs(f(c)) < tol) {
            *root = c;
            return true;
        }
        if (f(a) * f(c) < 0) {
            b = c;
        } else {
            a = c;
        }
    }

    return false;
}

bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root) {
    const double tol = 1e-7;
    const int max_iter = 1000000;
    
    if (f(a) * f(b) >= 0) {
        return false; // No sign change
    }
    for (int i = 0; i < max_iter; ++i) {
        double c = a - (f(a) * (b - a)) / (f(b) - f(a));
        
        if (std::abs(f(c)) < tol) {
            *root = c;
            return true;
        }
        
        if (f(a) * f(c) < 0) {
            b = c;
        } 
        if (f(b) * f(c) < 0) {
            a = c;
        }
        
        if (std::abs(b-a) < tol) {
            *root = c;
            return true;
        }
        
    }

    return false;
}

bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root) {
    const double tol = 1e-7;
    const int max_iter = 1000;
    
    if (f(a) * f(b) >= 0) {
        return false; // No sign change
    }

    for (int i = 0; i < max_iter; ++i) {
        double fc = f(c);
        double gc = g(c);
        
        if (std::abs(fc) < tol) {
            *root = c;
            return true;
        }
        
        if (gc == 0) {
            return false; // Denominator is zero
        }
        
        double c_new = c - fc / gc;
        
        if (c_new < a || c_new > b) {
            c_new = (a + b) / 2; // Fallback to midpoint
        }
        if (std::abs(c_new - c) < tol) {
            *root = c_new;
            return true;
        }
        if (f(a) * f(c_new) < 0) {
            b = c_new;
        } else {
            a = c_new;
        }
        
        c = c_new;
    }

    return false;
}

bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root) {
    const double tol = 1e-7;
    const int max_iter = 1000;

    double x_prev = a;
    double x_curr = b;

    for (int i = 0; i < max_iter; ++i) {
        double f_prev = f(x_prev);
        double f_curr = f(x_curr);

        if (std::abs(f_curr) < tol) {
            *root = x_curr;
            return true;
        }

        if (f_curr - f_prev == 0) {
            return false; // Denominator is zero
        }

        double x_new = x_curr - f_curr * (x_curr - x_prev) / (f_curr - f_prev);

        if (x_new < a || x_new > b) {
            x_new = (a + b) / 2; // Fallback to midpoint
        }

        if (std::abs(x_new - x_curr) < tol) {
            *root = x_new;
            return true;
        }

        x_prev = x_curr;
        x_curr = x_new;
    }
    
    return false;
}

