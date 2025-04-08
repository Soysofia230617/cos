#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <sstream>

// Структура для double-double арифметики
struct DoubleDouble {
    double hi;
    double lo;

    DoubleDouble(double h = 0.0, double l = 0.0) : hi(h), lo(l) {}
};

// Точные базовые операции
DoubleDouble two_sum(double a, double b) {
    double s = a + b;
    double v = s - a;
    double e = (a - (s - v)) + (b - v);
    return {s, e};
}

DoubleDouble two_diff(double a, double b) {
    double s = a - b;
    double v = s - a;
    double e = (a - (s - v)) - (b + v);
    return {s, e};
}

DoubleDouble two_prod(double a, double b) {
    double p = a * b;
    double e = std::fma(a, b, -p);
    return {p, e};
}

DoubleDouble add_dd(DoubleDouble a, DoubleDouble b) {
    DoubleDouble s = two_sum(a.hi, b.hi);
    DoubleDouble t = two_sum(a.lo, b.lo);
    s.lo += t.hi;
    s = two_sum(s.hi, s.lo);
    s.lo += t.lo;
    return two_sum(s.hi, s.lo);
}

DoubleDouble sub_dd(DoubleDouble a, DoubleDouble b) {
    DoubleDouble s = two_diff(a.hi, b.hi);
    DoubleDouble t = two_diff(a.lo, b.lo);
    s.lo += t.hi;
    s = two_sum(s.hi, s.lo);
    s.lo += t.lo;
    return two_sum(s.hi, s.lo);
}

DoubleDouble mul_dd(DoubleDouble a, DoubleDouble b) {
    DoubleDouble p1 = two_prod(a.hi, b.hi);
    double p2 = a.hi * b.lo;
    double p3 = a.lo * b.hi;
    DoubleDouble s = two_sum(p1.lo, p2);
    s.lo += p3;
    return two_sum(p1.hi, s.lo);
}

DoubleDouble div_dd(DoubleDouble a, DoubleDouble b) {
    double q1 = a.hi / b.hi;
    DoubleDouble r = sub_dd(a, mul_dd(b, {q1, 0.0}));
    double q2 = r.hi / b.hi;
    DoubleDouble s = two_sum(q1, q2);
    r = sub_dd(r, mul_dd(b, {q2, 0.0}));
    double q3 = r.hi / b.hi;
    return two_sum(s.hi, q3);
}

DoubleDouble neg_dd(DoubleDouble a) {
    return {-a.hi, -a.lo};
}

bool gt_dd(DoubleDouble a, DoubleDouble b) {
    if (a.hi > b.hi) return true;
    if (a.hi < b.hi) return false;
    return a.lo > b.lo;
}

bool lt_dd(DoubleDouble a, DoubleDouble b) {
    if (a.hi < b.hi) return true;
    if (a.hi > b.hi) return false;
    return a.lo < b.lo;
}

bool geq_dd(DoubleDouble a, DoubleDouble b) {
    if (a.hi > b.hi) return true;
    if (a.hi < b.hi) return false;
    return a.lo >= b.lo;
}

// Константы в формате DoubleDouble
const DoubleDouble PI = {3.14159265358979323846, 1.2246467991473532072e-16};
const DoubleDouble PI_2 = {1.57079632679489661923, 6.1232339957367660359e-17};
const DoubleDouble TWO_PI = {6.28318530717958647692, 2.4492935982947064143e-16};
const int TABLE_SIZE = 524288;

// Структура для таблицы
struct TableEntry {
    DoubleDouble x; // Координата узла (в радианах)
    double cos;     // Значение cos(x)
    double sin;     // Значение sin(x)
};

std::vector<TableEntry> trig_table;

// Параметры CORDIC
const int CORDIC_ITERATIONS = 50;
const DoubleDouble CORDIC_K = {0.6072529350088812561694, 0.0};
const double CORDIC_ANGLES[] = {
        0.7853981633974483, 0.4636476090008061, 0.24497866312686414,
        0.12435499454676144, 0.06241880999595735, 0.031239833430268277,
        0.015623728620476831, 0.007812341060101111, 0.0039062301319669718,
        0.0019531225164788188, 0.0009765621895593195, 0.0004882812111948983,
        0.0002441406201493618, 0.0001220703118936702, 0.00006103515617420877,
        0.0000305175781155261, 0.00001525878906131576, 0.00000762939453110197,
        0.00000381469726560650, 0.00000190734863281019
};

// Точный atan для DoubleDouble
DoubleDouble atan_dd(DoubleDouble x) {
    // Простое приближение для малых x
    if (lt_dd(x, {1e-8, 0.0})) return x;
    // Для больших значений нужно использовать более сложные методы
    // Здесь используем стандартный atan для hi части
    double at = std::atan(x.hi);
    return {at, 0.0}; // Упрощение, для точности 
}

void cordic(DoubleDouble theta, DoubleDouble* cos_val, DoubleDouble* sin_val) {
    DoubleDouble x = CORDIC_K;
    DoubleDouble y = {0.0, 0.0};
    DoubleDouble z = theta;

    if (gt_dd(z, PI_2)) {
        z = sub_dd(z, PI);
        x.hi = -x.hi;
        x.lo = -x.lo;
    } else if (lt_dd(z, neg_dd(PI_2))) {
        z = add_dd(z, PI);
        x.hi = -x.hi;
        x.lo = -x.lo;
    }

    for (int i = 0; i < CORDIC_ITERATIONS; i++) {
        DoubleDouble d = geq_dd(z, {0.0, 0.0}) ? DoubleDouble{1.0, 0.0} : DoubleDouble{-1.0, 0.0};
        DoubleDouble pow2 = div_dd({1.0, 0.0}, {static_cast<double>(1LL << i), 0.0});

        DoubleDouble x_new = sub_dd(x, mul_dd(y, mul_dd(d, pow2)));
        DoubleDouble y_new = add_dd(y, mul_dd(x, mul_dd(d, pow2)));
        DoubleDouble angle = (i < sizeof(CORDIC_ANGLES) / sizeof(CORDIC_ANGLES[0]) ?
                              DoubleDouble{CORDIC_ANGLES[i], 0.0} : atan_dd(pow2));
        z = sub_dd(z, mul_dd(d, angle));

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
        double theta = (2 * i + 1) * PI.hi / (2 * TABLE_SIZE);
        double cos_theta = std::cos(theta);
        DoubleDouble x = mul_dd(PI_2, div_dd(add_dd({1.0, 0.0}, {cos_theta, 0.0}), {2.0, 0.0}));
        trig_table[i].x = x;

        DoubleDouble cos_val, sin_val;
        cordic(x, &cos_val, &sin_val);
        trig_table[i].cos = cos_val.hi;
        trig_table[i].sin = sin_val.hi;

        if (i % (TABLE_SIZE / 50) == 0) {
            std::cout << "=";
            std::cout.flush();
        }
    }

    std::sort(trig_table.begin(), trig_table.end(),
              [](const TableEntry& a, const TableEntry& b) {
                  return lt_dd(a.x, b.x);
              });

    std::cout << "] 100%\nТаблица успешно инициализирована!\n\n";
}

DoubleDouble lagrange_interpolate(DoubleDouble x, bool use_cos = true, int degree = 3) {
    auto it = std::lower_bound(trig_table.begin(), trig_table.end(), x,
                               [](const TableEntry& entry, DoubleDouble val) {
                                   return lt_dd(entry.x, val);
                               });

    int center_idx = std::distance(trig_table.begin(), it);
    int start_idx = std::max(0, center_idx - degree / 2);
    int end_idx = std::min(TABLE_SIZE - 1, start_idx + degree);
    start_idx = end_idx - degree;

    DoubleDouble result = {0.0, 0.0};

    std::cout << "Используемые узлы интерполяции:\n";
    for (int i = start_idx; i <= end_idx; ++i) {
        std::cout << "  x[" << i << "] = " << std::setprecision(16) << trig_table[i].x.hi
                  << ", cos(x) = " << trig_table[i].cos << "\n";
    }
    std::cout << "\n";

    for (int i = start_idx; i <= end_idx; ++i) {
        double y = use_cos ? trig_table[i].cos : trig_table[i].sin;
        DoubleDouble weight = {1.0, 0.0};

        for (int j = start_idx; j <= end_idx; ++j) {
            if (j != i) {
                DoubleDouble num = sub_dd(x, trig_table[j].x);
                DoubleDouble den = sub_dd(trig_table[i].x, trig_table[j].x);
                weight = mul_dd(weight, div_dd(num, den));
            }
        }

        std::cout << "Узел " << i << ": вес = " << std::setprecision(16) << weight.hi << "\n";

        DoubleDouble term = mul_dd({y, 0.0}, weight);
        result = add_dd(result, term);
    }

    return result;
}

void range_reduce(DoubleDouble x_deg, int* quadrant, DoubleDouble* reduced) {
    // Преобразование градусов в радианы: x_rad = x_deg * (PI / 180)
    DoubleDouble one_eighty = {180.0, 0.0};
    DoubleDouble x_rad = mul_dd(x_deg, div_dd(PI, one_eighty));

    // Деление на 2π
    DoubleDouble n = {std::floor(x_rad.hi / TWO_PI.hi), 0.0}; // Упрощение
    DoubleDouble remainder = sub_dd(x_rad, mul_dd(n, TWO_PI));

    // Определение квадранта
    DoubleDouble q = div_dd(remainder, PI_2);
    *quadrant = static_cast<int>(q.hi) % 4;
    if (*quadrant < 0) *quadrant += 4;

    // Приведённый угол
    *reduced = sub_dd(remainder, mul_dd({static_cast<double>(*quadrant), 0.0}, PI_2));
}

DoubleDouble cos_dd(DoubleDouble x_deg, int degree) {
    static bool initialized = false;
    if (!initialized) {
        init_tables();
        initialized = true;
    }

    if (x_deg.hi == 0.0 && x_deg.lo == 0.0) return {1.0, 0.0};

    int quadrant;
    DoubleDouble reduced;
    range_reduce(x_deg, &quadrant, &reduced);

    std::cout << "\nВычисление cos(" << x_deg.hi << "°):\n";
    std::cout << "1. Редукция аргумента:\n";
    std::cout << "   Исходный x (градусы): " << std::setprecision(16) << x_deg.hi << "\n";
    std::cout << "   Исходный x (радианы): " << mul_dd(x_deg, div_dd(PI, {180.0, 0.0})).hi << "\n";
    std::cout << "   Приведенный x (радианы): " << reduced.hi << " (квадрант " << quadrant << ")\n";

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

        DoubleDouble x_dd = {x_deg, 0.0};
        DoubleDouble result = cos_dd(x_dd, degree);
        double x_rad = x_deg * PI.hi / 180.0;
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
    interactive_mode();

    std::cout << "\nПрограмма завершена.\n";
    return 0;
}
