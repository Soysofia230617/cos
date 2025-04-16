#include <iostream>
#include <cmath>
#include <cassert>
#include <limits>

using namespace std;

struct DoubleDouble {
    double hi; // старшая часть
    double lo; // младшая часть
};

// Константа для числа π
const double m_PI = 3.141592653589793238462643383279502884197169399375105820974944;

// Функции для работы с Double-Double числами

// Преобразование int в DoubleDouble
DoubleDouble dd_from_int(int x) {
    return { static_cast<double>(x), 0.0 };
}

// Преобразование double в DoubleDouble
DoubleDouble dd_from_double(double x) {
    return { x, 0.0 };
}

// Сложение DoubleDouble чисел
DoubleDouble dd_add(const DoubleDouble& a, const DoubleDouble& b) {
    double hi = a.hi + b.hi;
    double lo = a.lo + b.lo + ((a.hi - hi) + (b.hi - b.hi));
    return { hi, lo };
}

// Вычитание DoubleDouble чисел
DoubleDouble dd_sub(const DoubleDouble& a, const DoubleDouble& b) {
    double hi = a.hi - b.hi;
    double lo = a.lo - b.lo + ((a.hi - hi) - (b.hi - b.hi));
    return { hi, lo };
}

// Умножение DoubleDouble чисел
DoubleDouble dd_mul(const DoubleDouble& a, const DoubleDouble& b) {
    double hi = a.hi * b.hi;
    double lo = (a.hi * b.lo) + (a.lo * b.hi);
    return { hi, lo };
}

// Деление DoubleDouble чисел
DoubleDouble dd_div(const DoubleDouble& a, const DoubleDouble& b) {
    double hi = a.hi / b.hi;
    double lo = ((a.hi - hi * b.hi) + a.lo) / b.hi;
    return { hi, lo };
}

// Абсолютное значение DoubleDouble числа
DoubleDouble dd_abs(const DoubleDouble& a) {
    return { std::abs(a.hi), std::abs(a.lo) };
}

// Функция модуля для DoubleDouble чисел
DoubleDouble dd_mod(const DoubleDouble& a, const DoubleDouble& b) {
    DoubleDouble result = a;
    while (result.hi >= b.hi) {
        result = dd_sub(result, b);
    }
    return result;
}

// Основная функция вычисления косинуса с помощью DoubleDouble
DoubleDouble dd_cos(const DoubleDouble& x) {
    if (x.lo == 0.0) {
        return dd_from_double(std::cos(x.hi));
    }

    DoubleDouble r = dd_mod(x, dd_from_double(2.0 * M_PI)); // Ожидаем значение в пределах 0..2π

    if (r.hi > M_PI) {
        r = dd_sub(r, dd_from_double(2.0 * M_PI)); // Сворачиваем значение в диапазон [-π, π]
    }
    else if (r.hi < -M_PI) {
        r = dd_add(r, dd_from_double(2.0 * M_PI));
    }

    DoubleDouble sum = dd_from_double(1.0);
    DoubleDouble term = dd_from_double(1.0);
    DoubleDouble r2 = dd_mul(r, r);
    int n = 0;
    const double tol = 1e-40; // Уменьшаем погрешность

    // Суммирование членов ряда Тейлора до нужной точности
    while (true) {
        n++;
        DoubleDouble factor = dd_div(dd_from_double(-1.0), dd_from_int((2 * n - 1) * (2 * n)));
        term = dd_mul(term, dd_mul(r2, factor));
        sum = dd_add(sum, term);

        if (dd_abs(term).hi < tol) {
            break;
        }

        if (n > 500) {
            std::cerr << "Предупреждение: достигнут лимит итераций (" << n << ")\n";
            break;
        }
    }
    return sum;
}

int main() {
    cout.precision(40);
    double x;

    while (true) {
        cout << "Enter x (Ctrl+D/Ctrl+Z to exit): ";
        if (cin >> x) {
            DoubleDouble input = dd_from_double(x);
            DoubleDouble my_result = dd_cos(input);
            double std_result = std::cos(x);
            double dd_full = my_result.hi + my_result.lo;
            double abs_error = std::abs(dd_full - std_result);

            cout << "dd_cos( x )  = " << dd_full << endl;
            cout << "std::cos( x ) = " << std_result << endl;
            cout << "abs error             = " << abs_error << endl;
            cout << "------------------------------------------" << endl;
        } else if (cin.eof()) {
            // Выход при Ctrl+D/Ctrl+Z
            break;
        } else {
            cout << "Ошибка: введите число\n";
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
        }
    }

    return 0;
}
