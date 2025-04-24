#include <iostream>
#include <cmath>
#include <vector>
#include <limits>
#include <iomanip>
#include <fstream>

long double airy_ai_numeric(long double x) {
    const long double PI = acosl(-1.0L);
    const long double T  = 100.0L;   
    const int    N  = 5000000;    
    const long double h = T / N;

    long double sum = 0.0L;
    for (int i = 0; i <= N; ++i) {
        long double t = i * h;
        long double f = cosl(t*t*t/3.0L + x*t);
        if (i == 0 || i == N) {
            sum += f;
        } 
        else if (i & 1) {
            sum += 4.0L * f;
        } 
        else {
            sum += 2.0L * f;
        }
    }
    long double integral = sum * h / 3.0L;
    return integral / PI;
}

long double secant_method(long double (*f)(long double), long double x0, long double x1, long double tol = 1e-18L,int max_iter = 100) {
    long double f0 = f(x0);
    long double f1 = f(x1);
    for (int iter = 0; iter < max_iter; ++iter) {
        long double denom = f1 - f0;
        if (fabsl(denom) < std::numeric_limits<long double>::epsilon()) {
            return std::numeric_limits<long double>::quiet_NaN();
        }
        long double x2 = x1 - f1 * (x1 - x0) / denom;
        if (fabsl(x2 - x1) < tol) {
            std::cout << "Quantity iterations for root " << x2 << " is " << iter << "\n";
            return x2;
        }
        x0 = x1; 
        f0 = f1;
        x1 = x2; 
        f1 = f(x1);
    }
    return std::numeric_limits<long double>::quiet_NaN();
}

std::vector<long double> find_airy_zeros(long double a, long double b, long double scan_step = 0.1L, long double tol = 1e-18L, int max_iter = 100) {
    std::vector<long double> roots;
    long double x = a;
    while (x + scan_step <= b) {
        long double f1 = airy_ai_numeric(x);
        long double f2 = airy_ai_numeric(x + scan_step);
        if (f1 * f2 < 0.0L) {
            long double root = secant_method(airy_ai_numeric, x, x + scan_step, tol, max_iter);
            if (!std::isnan(root)) {
                bool unique = true;
                for (auto r : roots) {
                    if (fabsl(r - root) < tol * 10) {
                        unique = false;
                        break;
                    }
                }
                if (unique) {
                    roots.push_back(root);
                    std::cout << "Found root: x = " << std::setprecision(18) << root << ", Ai(x) = " << airy_ai_numeric(root) << "\n";
                }
            }
        }
        x += scan_step;
    }
    return roots;
}

int main() {
    long double a = -15.0L;
    long double b = 5.0L;
    long double dx = 0.1L;
    std::ofstream fout("airy_ai_values.txt");
    fout << std::fixed << std::setprecision(18);

    for (long double x = a; x <= b; x += dx) {
        fout << x << " " << airy_ai_numeric(x) << "\n";
    }
    fout.close();

    a = -15.0L;
    b = -1.0L;
    long double scan_step = 0.1L;
    long double tol = 1e-16L;
    int max_iter = 100;
    std::cout << std::fixed << std::setprecision(18);
    std::cout << "\nSearching for roots of Ai(x)=0 in [" << a << ", " << b << "]\n";
    auto roots = find_airy_zeros(a, b, scan_step, tol, max_iter);
    std::cout << "\nTotal roots found: " << roots.size() << "\n";

    return 0;
}
