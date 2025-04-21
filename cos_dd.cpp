#include <iostream>
#include <cmath>
#include <iomanip>

// Структура для double-double числа
struct DoubleDouble {
    double hi, lo;
    DoubleDouble(double h = 0.0, double l = 0.0) : hi(h), lo(l) {}
};

// TwoSum: Полное сложение двух double
void twoSum(double a, double b, double& s, double& e) {
    s = a + b;
    double v = s - a;
    e = (a - (s - v)) + (b - v);
}

// TwoDiff: Полное вычитание двух double
void twoDiff(double a, double b, double& s, double& e) {
    s = a - b;
    double v = s - a;
    e = (a - (s - v)) - (b + v);
}

// QuickTwoSum: Для нормализации
void quickTwoSum(double a, double b, double& s, double& e) {
    s = a + b;
    e = b - (s - a);
}

// TwoProd: Умножение двух double с учётом ошибки
void twoProd(double a, double b, double& p, double& e) {
    p = a * b;
    e = std::fma(a, b, -p);
}

// Сложение double-double
DoubleDouble operator+(const DoubleDouble& a, const DoubleDouble& b) {
    double s1, s2, t1, t2, v1, v2;
    twoSum(a.hi, b.hi, s1, t1);
    twoSum(a.lo, b.lo, s2, t2);
    s2 += t1;
    quickTwoSum(s1, s2, s1, v1);
    v1 += t2;
    quickTwoSum(s1, v1, s1, v2);
    return DoubleDouble(s1, v2);
}

// Вычитание double-double
DoubleDouble operator-(const DoubleDouble& a, const DoubleDouble& b) {
    double s1, s2, t1, t2, v1, v2;
    twoDiff(a.hi, b.hi, s1, t1);
    twoDiff(a.lo, b.lo, s2, t2);
    s2 += t1;
    quickTwoSum(s1, s2, s1, v1);
    v1 += t2;
    quickTwoSum(s1, v1, s1, v2);
    return DoubleDouble(s1, v2);
}

// Умножение double-double с TwoProd
DoubleDouble operator*(const DoubleDouble& a, const DoubleDouble& b) {
    double p1, e1, p2;
    twoProd(a.hi, b.hi, p1, e1);
    p2 = a.hi * b.lo + a.lo * b.hi + e1 + a.lo * b.lo;
    return DoubleDouble(p1, p2);
}

// Деление double-double на double
DoubleDouble operator/(const DoubleDouble& a, double b) {
    double q1 = a.hi / b;
    DoubleDouble prod = DoubleDouble(q1, 0.0) * DoubleDouble(b, 0.0);
    DoubleDouble diff = a - prod;
    double q2 = diff.hi / b;
    quickTwoSum(q1, q2, q1, q2);
    return DoubleDouble(q1, q2);
}

// Вычисление косинуса через ряд Тейлора
DoubleDouble cos_dd(DoubleDouble x) {
    // Точные константы
    const DoubleDouble pi_dd(3.14159265358979323846, 1.2246467991473532072e-16);
    const DoubleDouble two_pi = pi_dd + pi_dd;
    const DoubleDouble pi_over_2 = pi_dd / 2.0;
    const DoubleDouble three_pi_over_2 = pi_over_2 + pi_dd;

    // Приведение x к [0, 2π]
    while (x.hi >= 2.0 * M_PI) x = x - two_pi;
    while (x.hi < 0.0) x = x + two_pi;

    // Проверка на x ≈ π/2
    DoubleDouble diff = x - pi_over_2;
    if (std::abs(diff.hi) < 1e-15) {
        return DoubleDouble(0.0, 0.0); // cos(π/2) = 0
    }

    // Приведение к [0, π/2]
    bool negate = false;
    if (x.hi > M_PI_2 && x.hi <= M_PI) {
        x = pi_dd - x; // cos(x) = -cos(π - x)
        negate = true; // Исправлено: знак меняется
    } else if (x.hi > M_PI && x.hi <= three_pi_over_2.hi) {
        x = x - pi_dd; // cos(x) = -cos(x - π)
        negate = true;
    } else if (x.hi > three_pi_over_2.hi && x.hi < 2.0 * M_PI) {
        x = two_pi - x; // cos(x) = -cos(2π - x)
        negate = true;
    }

    // Ряд Тейлора для cos(x)
    DoubleDouble sum(1.0, 0.0);
    DoubleDouble term(1.0, 0.0);
    DoubleDouble x2 = x * x;
    double n = 1.0;

    for (int i = 1; i < 100; ++i) {
        term = term * x2;
        term = term / (2.0 * n * (2.0 * n - 1.0));
        term = term * DoubleDouble(-1.0, 0.0);
        sum = sum + term;
        n += 1.0;
        if (std::abs(term.hi) < 1e-40) break;
    }

    return negate ? (DoubleDouble(-1.0, 0.0) * sum) : sum;
}

int main() {
    double x;
    std::cout << "Enter x (in radians): ";
    if (!(std::cin >> x)) {
        std::cerr << "Invalid input" << std::endl;
        return 1;
    }

    DoubleDouble xd(x, 0.0);
    DoubleDouble result = cos_dd(xd);

    std::cout << std::setprecision(40);
    std::cout << "cos(" << x << ") = " << result.hi << " (hi) + " << result.lo << " (lo)" << std::endl;

    return 0;
}