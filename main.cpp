#include <iostream>
#include <cmath>
#include <cassert>
#include <limits>
#include <iomanip>
#include <cstdint>
using namespace std;

struct DoubleDouble {
    double hi; // старшая часть
    double lo; // младшая часть
};

// Константа для числа π
const double m_PI = 3.141592653589793238462643383279502884197169399375105820974944;

// Функция для вывода двоичного представления числа IEEE 754
void printBinaryIEEE754(double x) {
    // Преобразуем double в 64-битное представление
    union {
        double d;
        uint64_t bits;
    } u;
    u.d = x;

    // Извлекаем знак, экспоненту и мантиссу
    uint64_t sign = (u.bits >> 63) & 1;
    uint64_t exponent = (u.bits >> 52) & 0x7FF; // 11 бит
    uint64_t mantissa = u.bits & 0xFFFFFFFFFFFFF; // 52 бита

    cout << "Binary IEEE 754 representation:\n";
    cout << "Sign: " << sign << "\n";
    cout << "Exponent (biased, 11 bits): ";
    for (int i = 10; i >= 0; --i) {
        cout << ((exponent >> i) & 1);
    }
    cout << " (decimal: " << exponent << ", actual: " << (exponent - 1023) << ")\n";
    cout << "Mantissa (52 bits + implicit 1): 1.";
    for (int i = 51; i >= 0; --i) {
        cout << ((mantissa >> i) & 1);
        if (i % 4 == 0) cout << " "; // Для читаемости
    }
    cout << "\n";
}

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
    double v = hi - a.hi;
    double lo = (a.hi - (hi - v)) + (b.hi - v);
    lo += a.lo + b.lo;
    double t = hi + lo;
    lo = lo - (t - hi);
    hi = t;
    return { hi, lo };
}

// Вычитание DoubleDouble чисел
DoubleDouble dd_sub(const DoubleDouble& a, const DoubleDouble& b) {
    double hi = a.hi - b.hi;
    double v = hi - a.hi;
    double lo = (a.hi - (hi - v)) - (b.hi + v);
    lo += a.lo - b.lo;
    double t = hi + lo;
    lo = lo - (t - hi);
    hi = t;
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
    DoubleDouble r = dd_mod(x, dd_from_double(2.0 * m_PI)); // Ожидаем значение в пределах 0..2π

    if (r.hi > m_PI) {
        r = dd_sub(r, dd_from_double(2.0 * m_PI)); // Сворачиваем значение в диапазон [-π, π]
    }
    else if (r.hi < -m_PI) {
        r = dd_add(r, dd_from_double(2.0 * m_PI));
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

// Функция для вывода значений DoubleDouble с заданной точностью
void printDoubleDouble(const DoubleDouble& dd, int precision) {
    cout  << setprecision(precision) << dd.hi << endl;
    cout  << setprecision(precision) << dd.lo << endl;
    cout  << setprecision(precision) << fmod(dd.hi * pow(10, 16), 10.) + dd.lo * pow(10, 16) << endl;
}

// Функция тестирования операций DoubleDouble
void testDoubleDoubleOperations(int precision) {
    cout << "Testing DoubleDouble Arithmetic Operations\n";
    cout << "=======================================\n";
    cout.precision(40);

    // Test 1: Addition
    {
        cout << "\n=== Addition Test ===\n";
        DoubleDouble a = dd_from_double(1.234567890123456789);
        DoubleDouble b = dd_from_double(2.345678901234567890);
        DoubleDouble result = dd_add(a, b);
        double expected = 1.234567890123456789 + 2.345678901234567890;
        double actual = result.hi + result.lo;

        cout << "Input a:\n";
        printBinaryIEEE754(a.hi);
        cout << "Input b:\n";
        printBinaryIEEE754(b.hi);
        cout << "DoubleDouble result:\n";
        printDoubleDouble(result, precision);
        cout << "Binary representation of result.hi:\n";
        printBinaryIEEE754(result.hi);
        cout << "Binary representation of result.lo:\n";
        printBinaryIEEE754(result.lo);
        cout << "Standard double result: " << expected << "\n";
        cout << "Absolute error: " << std::abs(actual - expected) << "\n";
    }
    // Test 2: Subtraction
    {
        cout << "\n=== Subtraction Test ===\n";
        DoubleDouble a = dd_from_double(3.14159265358979323846);
        DoubleDouble b = dd_from_double(1.41421356237309504880);
        DoubleDouble result = dd_sub(a, b);
        double expected = 3.14159265358979323846 - 1.41421356237309504880;
        double actual = result.hi + result.lo;

        cout << "Input a:\n";
        printBinaryIEEE754(a.hi);
        cout << "Input b:\n";
        printBinaryIEEE754(b.hi);
        cout << "DoubleDouble result:\n";
        printDoubleDouble(result, precision);
        cout << "Binary representation of result.hi:\n";
        printBinaryIEEE754(result.hi);
        cout << "Binary representation of result.lo:\n";
        printBinaryIEEE754(result.lo);
        cout << "Standard double result: " << expected << "\n";
        cout << "Absolute error: " << std::abs(actual - expected) << "\n";
    }
    // Test 3: Multiplication
    {
        cout << "\n=== Multiplication Test ===\n";
        DoubleDouble a = dd_from_double(2.71828182845904523536);
        DoubleDouble b = dd_from_double(1.41421356237309504880);
        DoubleDouble result = dd_mul(a, b);
        double expected = 2.71828182845904523536 * 1.41421356237309504880;
        double actual = result.hi + result.lo;

        cout << "Input a:\n";
        printBinaryIEEE754(a.hi);
        cout << "Input b:\n";
        printBinaryIEEE754(b.hi);
        cout << "DoubleDouble result:\n";
        printDoubleDouble(result, precision);
        cout << "Binary representation of result.hi:\n";
        printBinaryIEEE754(result.hi);
        cout << "Binary representation of result.lo:\n";
        printBinaryIEEE754(result.lo);
        cout << "Standard double result: " << expected << "\n";
        cout << "Absolute error: " << std::abs(actual - expected) << "\n";
    }
    // Test 4: Division
    {
        cout << "\n=== Division Test ===\n";
        DoubleDouble a = dd_from_double(3.14159265358979323846);
        DoubleDouble b = dd_from_double(2.0);
        DoubleDouble result = dd_div(a, b);
        double expected = 3.14159265358979323846 / 2.0;
        double actual = result.hi + result.lo;

        cout << "Input a:\n";
        printBinaryIEEE754(a.hi);
        cout << "Input b:\n";
        printBinaryIEEE754(b.hi);
        cout << "DoubleDouble result:\n";
        printDoubleDouble(result, precision);
        cout << "DoubleDouble result:\n";
        printDoubleDouble(result, precision);
        cout << "Binary representation of result.hi:\n";
        printBinaryIEEE754(result.hi);
        cout << "Binary representation of result.lo:\n";
        printBinaryIEEE754(result.lo);
        cout << "Standard double result: " << expected << "\n";
        cout << "Absolute error: " << std::abs(actual - expected) << "\n";
    }
    // Test 5: Absolute Value
    {
        cout << "\n=== Absolute Value Test ===\n";
        DoubleDouble a = dd_from_double(-2.71828182845904523536);
        DoubleDouble result = dd_abs(a);
        double expected = std::abs(-2.71828182845904523536);
        double actual = result.hi + result.lo;

        cout << "Input a:\n";
        printBinaryIEEE754(a.hi);
        cout << "DoubleDouble result:\n";
        printDoubleDouble(result, precision);
        cout << "Binary representation of result.hi:\n";
        printBinaryIEEE754(result.hi);
        cout << "Binary representation of result.lo:\n";
        printBinaryIEEE754(result.lo);
        cout << "Standard double result: " << expected << "\n";
        cout << "Absolute error: " << std::abs(actual - expected) << "\n";
    }
    // Test 6: Modulo Operation
    {
        cout << "\n=== Modulo Operation Test ===\n";
        DoubleDouble a = dd_from_double(10.5);
        DoubleDouble b = dd_from_double(3.0);
        DoubleDouble result = dd_mod(a, b);
        double expected = fmod(10.5, 3.0);
        double actual = result.hi + result.lo;

        cout << "Input a:\n";
        printBinaryIEEE754(a.hi);
        cout << "Input b:\n";
        printBinaryIEEE754(b.hi);
        cout << "DoubleDouble result:\n";
        printDoubleDouble(result, precision);
        cout << "Binary representation of result.hi:\n";
        printBinaryIEEE754(result.hi);
        cout << "Binary representation of result.lo:\n";
        printBinaryIEEE754(result.lo);
        cout << "Standard double result: " << expected << "\n";
        cout << "Absolute error: " << std::abs(actual - expected) << "\n";
    }
    // Test 7: Binary Representation of a Special Case (small number)
    {
        cout << "\n=== Small Number Binary Representation Test ===\n";
        double small = 1.0e-300; // Очень малое число
        DoubleDouble a = dd_from_double(small);
        cout << "Input small number: " << setprecision(20) << small << "\n";
        printBinaryIEEE754(a.hi);
        cout << "DoubleDouble representation:\n";
        printDoubleDouble(a, precision);
    }
}
int computeErrorOrder(const DoubleDouble& error) {
    double abs_error = std::abs(error.hi + error.lo);
    if (abs_error == 0.0) {
        return 0; // Порядок не определён для нулевой ошибки
    }
    return static_cast<int>(std::floor(std::log10(abs_error)));
}
int main() {
    int precision=19;


    // Запускаем тесты
    testDoubleDoubleOperations(precision);
    return 0;
    double x;
    cout << "Enter x (non-numeric input to exit): ";
    while (cin >> x) {
        DoubleDouble input = dd_from_double(x);
        DoubleDouble my_result = dd_cos(input);
        DoubleDouble std_result = dd_from_double(std::cos(x));
        DoubleDouble abs_error = dd_abs(dd_sub(my_result, std_result));
        int error_order = computeErrorOrder(abs_error);
        cout << "x = " << setprecision(precision) << x << endl;
        cout << "dd_cos(x):\n";
        printDoubleDouble(my_result, precision);
        cout << "std::cos(x):\n";
        printDoubleDouble(std_result, precision);
        cout << "Absolute error:\n";
        printDoubleDouble(abs_error, precision);
        if (error_order == 0 && (abs_error.hi + abs_error.lo) == 0.0) {
            cout << "Error order: undefined (error is zero)\n";
        } else {
            cout << "Error order: " << error_order << "\n";
        }
        cout << "------------------------------------------" << endl;
        cout << "Enter x (non-numeric input to exit): ";
    }

    return 0;
}
