#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <sstream>

const double PI = 3.14159265358979323846;
const double PI_2 = 1.57079632679489661923;
const int TABLE_SIZE = 524288;

struct DoubleDouble {
    double hi;
    double lo;
};

struct TableEntry {
    double x;    // Координата узла (в радианах)
    double cos;  // Значение cos(x)
    double sin;  // Значение sin(x)
};

std::vector<TableEntry> trig_table;

// Точные базовые операции
DoubleDouble two_sum(double a, double b) {
    double s = a + b;
    double v = s - a;
    double e = (a - (s - v)) + (b - v);
    return {s, e};
}

DoubleDouble two_prod(double a, double b) {
    double p = a * b;
    double e = std::fma(a, b, -p);
    return {p, e};
}

DoubleDouble add_dd(DoubleDouble a, DoubleDouble b) {
    DoubleDouble s = two_sum(a.hi, b.hi);
    s.lo += a.lo + b.lo;
    return two_sum(s.hi, s.lo);
}

// Параметры CORDIC
const int CORDIC_ITERATIONS = 50;
const double CORDIC_K = 0.6072529350088812561694;
const double CORDIC_ANGLES[] = {
        0.7853981633974483, 0.4636476090008061, 0.24497866312686414,
        0.12435499454676144, 0.06241880999595735, 0.031239833430268277,
        0.015623728620476831, 0.007812341060101111, 0.0039062301319669718,
        0.0019531225164788188, 0.0009765621895593195, 0.0004882812111948983,
        0.0002441406201493618, 0.0001220703118936702, 0.00006103515617420877,
        0.0000305175781155261, 0.00001525878906131576, 0.00000762939453110197,
        0.00000381469726560650, 0.00000190734863281019
};

void cordic(double theta, double* cos_val, double* sin_val) {
    double x = CORDIC_K;
    double y = 0.0;
    double z = theta;

    if (z > PI_2) {
        z -= PI;
        x = -x;
    } else if (z < -PI_2) {
        z += PI;
        x = -x;
    }

    for (int i = 0; i < CORDIC_ITERATIONS; i++) {
        double d = (z >= 0) ? 1.0 : -1.0;
        double pow2 = 1.0 / (1LL << i);

        double x_new = x - y * d * pow2;
        double y_new = y + x * d * pow2;
        z -= d * (i < sizeof(CORDIC_ANGLES) / sizeof(CORDIC_ANGLES[0]) ?
                  CORDIC_ANGLES[i] : atan(pow2));

        x = x_new;
        y = y_new;
    }

    *cos_val = x;
    *sin_val = y;
}

void init_tables() {
    trig_table.resize(TABLE_SIZE);

    std::cout << "Инициализация таблицы значений...\n";
    std::cout << "0% [";
    for (int i = 0; i < 50; ++i) std::cout << " ";
    std::cout << "] 100%\n0% [";
    std::cout.flush();

    for (int i = 0; i < TABLE_SIZE; ++i) {
        double theta = (2 * i + 1) * PI / (2 * TABLE_SIZE);
        trig_table[i].x = PI_2 / 2 * (1 + cos(theta));

        cordic(trig_table[i].x, &trig_table[i].cos, &trig_table[i].sin);

        if (i % (TABLE_SIZE / 50) == 0) {
            std::cout << "=";
            std::cout.flush();
        }
    }

    std::sort(trig_table.begin(), trig_table.end(),
              [](const TableEntry& a, const TableEntry& b) { return a.x < b.x; });

    std::cout << "] 100%\nТаблица успешно инициализирована!\n\n";
}

DoubleDouble lagrange_interpolate(double x, bool use_cos = true, int degree = 3) {
    auto it = std::lower_bound(trig_table.begin(), trig_table.end(), x,
                               [](const TableEntry& entry, double val) { return entry.x < val; });

    int center_idx = std::distance(trig_table.begin(), it);
    int start_idx = std::max(0, center_idx - degree / 2);
    int end_idx = std::min(TABLE_SIZE - 1, start_idx + degree);
    start_idx = end_idx - degree;

    DoubleDouble result = {0.0, 0.0};

    std::cout << "Используемые узлы интерполяции:\n";
    for (int i = start_idx; i <= end_idx; ++i) {
        std::cout << "  x[" << i << "] = " << std::setprecision(16) << trig_table[i].x
                  << ", cos(x) = " << trig_table[i].cos << "\n";
    }
    std::cout << "\n";

    for (int i = start_idx; i <= end_idx; ++i) {
        double y = use_cos ? trig_table[i].cos : trig_table[i].sin;
        double weight = 1.0;

        for (int j = start_idx; j <= end_idx; ++j) {
            if (j != i) {
                weight *= (x - trig_table[j].x) / (trig_table[i].x - trig_table[j].x);
            }
        }

        std::cout << "Узел " << i << ": вес = " << std::setprecision(16) << weight << "\n";

        DoubleDouble term = two_prod(y, weight);
        result = add_dd(result, term);
    }

    return result;
}

void range_reduce(double x_deg, int* quadrant, double* reduced) {
    double x_rad = x_deg * PI / 180.0; // Преобразование в радианы
    double n = std::floor(x_rad / (2 * PI));
    double remainder = x_rad - n * 2 * PI;

    *quadrant = static_cast<int>(remainder / PI_2) % 4;
    if (*quadrant < 0) *quadrant += 4;
    *reduced = remainder - (*quadrant) * PI_2;
}

DoubleDouble cos_dd(double x_deg, int degree) {
    static bool initialized = false;
    if (!initialized) {
        init_tables();
        initialized = true;
    }

    if (x_deg == 0.0) return {1.0, 0.0};

    int quadrant;
    double reduced;
    range_reduce(x_deg, &quadrant, &reduced);

    std::cout << "\nВычисление cos(" << x_deg << "°):\n";
    std::cout << "1. Редукция аргумента:\n";
    std::cout << "   Исходный x (градусы): " << std::setprecision(16) << x_deg << "\n";
    std::cout << "   Исходный x (радианы): " << x_deg * PI / 180.0 << "\n";
    std::cout << "   Приведенный x (радианы): " << reduced << " (квадрант " << quadrant << ")\n";

    DoubleDouble result;
    switch (quadrant) {
        case 0:
            std::cout << "2. Используем интерполяцию cos\n";
            result = lagrange_interpolate(reduced, true, degree);
            break;
        case 1:
            std::cout << "2. Используем интерполяцию sin с коррекцией знака\n";
            result = lagrange_interpolate(reduced, false, degree);
            result.hi = -result.hi;
            result.lo = -result.lo;
            break;
        case 2:
            std::cout << "2. Используем интерполяцию cos с коррекцией знака\n";
            result = lagrange_interpolate(reduced, true, degree);
            result.hi = -result.hi;
            result.lo = -result.lo;
            break;
        case 3:
            std::cout << "2. Используем интерполяцию sin\n";
            result = lagrange_interpolate(reduced, false, degree);
            break;
    }

    std::cout << "3. Итоговый результат: " << std::setprecision(16) << result.hi << " + " << result.lo << "\n\n";
    return result;
}

void interactive_mode() {
    std::cout << "\n=== Режим интерактивного ввода ===\n";
    std::cout << "Вводите значения x в градусах для вычисления cos(x)\n";
    std::cout << "Для выхода введите 'q'\n\n";

    while (true) {
        std::cout << "Введите x (градусы): ";
        std::string input;
        std::getline(std::cin, input);

        if (input == "q" || input == "Q") {
            break;
        }

        std::istringstream iss(input);
        double x_deg;
        int degree;
        if (!(iss >> x_deg)) {
            std::cout << "Ошибка: введите число или 'q' для выхода\n";
            continue;
        }

        std::cout << "Введите степень интерполяции (например, 3, 7, 15): ";
        std::getline(std::cin, input);
        iss.clear();
        iss.str(input);
        if (!(iss >> degree) || degree < 1 || degree > TABLE_SIZE) {
            std::cout << "Ошибка: степень должна быть числом от 1 до " << TABLE_SIZE << "\n";
            continue;
        }

        DoubleDouble result = cos_dd(x_deg, degree);
        double x_rad = x_deg * PI / 180.0;
        double exact = std::cos(x_rad);
        double error = std::abs(result.hi - exact);

        std::cout << std::setprecision(16);
        std::cout << "=== Результат ===\n";
        std::cout << "Наш алгоритм: " << result.hi << " + " << result.lo << "\n";
        std::cout << "Библиотечный cos: " << exact << "\n";
        std::cout << "Абсолютная ошибка: " << error << "\n";
        std::cout << "Относительная ошибка: " << (exact != 0 ? error / std::abs(exact) : 0) << "\n\n";
    }
}

int main() {
    std::cout << "Программа вычисления cos(x) с повышенной точностью\n";
    std::cout << "Используется:\n";
    std::cout << "- Чебышёвские узлы интерполяции\n";
    std::cout << "- Алгоритм CORDIC\n";
    std::cout << "- Арифметика двойной точности\n\n";

    interactive_mode();

    std::cout << "\nПрограмма завершена.\n";
    return 0;
}