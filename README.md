# Алгоритм вычисления косинуса с повышенной точностью

## Математические основы

### 1. Арифметика двойной точности (DoubleDouble)
Число представляется как пара чисел типа double:
```
value = hi + lo
```
где:
- `hi` содержит старшие значащие биты
- `lo` содержит младшие биты для коррекции точности

**Свойства**:
- Точность: ~30 десятичных знаков (вместо 15 у стандартного double)
- Ошибка операций: O(ε²) вместо O(ε)

### 2. Алгоритм CORDIC
Итерационный метод вычисления тригонометрических функций.

**Ключевые уравнения**:
```
xₙ₊₁ = xₙ - σₙ·yₙ·2⁻ⁿ
yₙ₊₁ = yₙ + σₙ·xₙ·2⁻ⁿ
zₙ₊₁ = zₙ - σₙ·atan(2⁻ⁿ)
```
где σₙ = sign(zₙ)

**Особенности**:
- Требует только сложения, вычитания и битовых сдвигов
- Сходится со скоростью ~1 бит за итерацию
- Требует предварительного вычисления констант atan(2⁻ⁿ)

### 3. Интерполяция Лагранжа
Строит многочлен степени N, точно проходящий через N+1 точку.

**Формула интерполяции**:
```
L(x) = Σ yᵢ · ℓᵢ(x)
```
где ℓᵢ(x) - базисные полиномы:
```
ℓᵢ(x) = Π (x - xⱼ)/(xᵢ - xⱼ) (j ≠ i)
```

**Ошибка интерполяции**:
Для cos(x) и кубической интерполяции (N=3):
```
|E(x)| ≤ |(x-x₀)(x-x₁)(x-x₂)(x-x₃)| / 24
```

## Работа алгоритма

### Этап 1. Предвычисление таблиц
1. Используя CORDIC, вычисляем значения cos/sin на равномерной сетке [0, π/2]
2. Сохраняем результаты в таблицах размером 524288 точек
3. Для каждого значения храним пару (hi, lo)

### Этап 2. Редукция аргумента
1. Приводим входной аргумент x к диапазону [0, 2π]:
   ```
   x' = x - 2π·round(x/2π)
   ```
2. Определяем квадрант:
   ```
   квадрант = floor(2x'/π)
   ```
3. Приводим к базовому интервалу [0, π/2]

### Этап 3. Интерполяция
1. Для приведенного аргумента θ находим 4 ближайших узла таблицы
2. Строим интерполяционный многочлен Лагранжа 3-й степени
3. Вычисляем значение многочлена в точке θ

### Этап 4. Коррекция результата
1. Учитываем квадрант:
- Q1:  cos(θ)
- Q2: -sin(θ)
- Q3: -cos(θ)
- Q4:  sin(θ)
2. Возвращаем результат в формате DoubleDouble

## Оценка точности

Компоненты ошибки:
1. **CORDIC**: ~2⁻⁵⁰ ≈ 10⁻¹⁵
2. **Интерполяция**: ~(π/2·2⁻¹⁹)⁴ ≈ 10⁻¹⁶
3. **Округление**: ~ε² ≈ 10⁻³²

**Итоговая точность**: ~15-16 верных десятичных знаков

## Оптимизации

1. **Таблицы**: Предвычисление исключает повторные расчеты
2. **Редукция**: Быстрое приведение аргумента за O(1)
3. **Интерполяция**: Локальная кубическая аппроксимация
4. **Параллелизм**: Независимость вычислений для разных x

## Ограничения

1. Требует значительной памяти для таблиц (~8MB)
2. Время инициализации ~0.1-1 сек
3. Максимальная точность ограничена точностью double
  - Формат IEEE 754 double precision обеспечивает ~15-17 значащих десятичных цифр
  - Эпсилон машины (ε) ≈ 2.22×10⁻¹⁶




