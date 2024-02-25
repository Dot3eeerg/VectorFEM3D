#import "function_template.typ": *
#import "nstu_template//template.typ": * // Импорт шаблона

// Пример конфигурации
// (многие параметры здесь не обязательны)
//#show: project.with(
//  faculty: "ФПМИ",
//  department: "Кафедра прикладной математики",
//  discipline: "Компьютерная графика",
//  task_type: "Практическое задание №4",
//  task_name: "kek СОЗДАНИЕ КРИВЫХ И ПОВЕРХНОСТЕЙ C ИСПОЛЬЗОВАНИЕМ СПЛАЙНОВ",
//  students: (
//    "Дзюбло Павел",
//    "Данченко Иван"
//  ),
//  teachers: (
//    "Задорожный A.Г.",
//    ""
//  ),
//  //variant: "",
//  group: "ПМ-02",
//  year: "2023",
//)

#show outline.entry.where(
  level: 1
) : it => {
  v(18pt, weak: true)
  set text(size: 16pt)
  strong(it)
}

#show outline.entry.where(
  level: 2
) : it => {
  v(10pt, weak: true)
  set text(size: 12pt)
  strong(it)
}

#show outline.entry.where(
  level: 3
) : it => {
  v(8pt, weak: true)
  set text(size: 10pt)
  strong(it)
}

#outline(title : [Содержание])

#let t01 = $t_(01)$
#let t02 = $t_(02)$
#let t03 = $t_(03)$

#let t12 = $t_(12)$
#let t13 = $t_(13)$
#let t23 = $t_(23)$

#let aVec = $limits(A)^(->)$

#pagebreak()

= Формулировка задания

Решить гиперболическую задачу c помощью векторного метода конечных элементов (далее МКЭ).

Векторный МКЭ для трёхмерной краевой задачи, уравнение гиперболического типа в декартовой системе координат. Четырёхслойная неявная схема по времени. Линейные базисные-вектор функции на параллелепипедах. Краевые условия первого рода.
 
= Постановка нестационарной задачи

Начально-краевая задача определяется уравнением:

$ "rot"(1/mu "rot" aVec) + sigma (diff aVec)/(diff t) + epsilon (diff^2 aVec) / (diff t^2) = limits(J^("ст"))^(->) $

$ lr((limits(A)^(->) times limits(n)^(->)) mid(|) \ _(S_1) = limits(A^g)^(->) times limits(n)^(->)), $

$ aVec_(t=t_0) = aVec_0. $

= Аппроксимация по времени

Обозначим через $aVec = aVec(x, y, z, t_j)$ значение вектор-потенциала на текущем временном слое, а через $limits(A^(arrow.l.double 1))^(->) = aVec(x, y, z, t_(j-1))$, $limits(A^(arrow.l.double 2))^(->) = aVec(x, y, z, t_(j-2))$ и $limits(A^(arrow.l.double 3))^(->) = aVec(x, y, z, t_(j-3))$ - значение вектор-потенциала на трёх предыдущих временных слоях. Тогда в результате аппроксимации по времени получим векторное уравнение:

$ "rot"(frac(1, mu) "rot"limits(A)^(->)) + gamma limits(A)^(->) = limits(F)^(->) $

где коэффициент $gamma$ и вектор-функция $limits(F)^(->)$ определяются схемой аппроксимации по времени. В этом случае используется четырёхслойная неявная схема. Сделаем замену обозначения временных промежутков для более удобной записи:

#v(0.2cm)

$ #grid(
  columns: (3cm, 3cm),
  gutter: 8pt,
  [$t01 = t_j - t_(j-1)$,], [$t12 = t_(j-1) - t_(j-2)$,],
  [$t02 = t_j - t_(j-2)$,], [$t13 = t_(j-1) - t_(j-3)$,],
  [$t03 = t_j - t_(j-3)$,], [$t23 = t_(j-2) - t_(j-3)$.]
) $

#v(0.5cm)

Коэффициент $gamma$ и вектор-функция $limits(F)^(->)$ имеют вид:

#v(0.2cm)

$ gamma = (sigma( t01 t02 + t01 t03 + t02 t03) + 2epsilon(t01 + t02 + t03))/(t01 t02 t03), $

#v(0.5cm)

$ limits(F)^(->) = limits(J^("ст"))^(->) + (sigma t02 t03 + 2epsilon(t02 + t03))/(t01 t12 t13) limits(A^(arrow.l.double 1))^(->) - (sigma t01 t03 + 2epsilon(t01 + t03))/(t02 t12 t23) limits(A^(arrow.l.double 2))^(->) + (sigma t01 t02 + 2epsilon(t01 + t02))/(t03 t13 t23) limits(A^(arrow.l.double 3))^(->). $

#pagebreak()

= Локальные матрицы и вектор

В одномерном случае линейные базисные функции имеют вид:

$ psi_1 = Nu_1(nu) = (nu_(r+1)-nu)/(h_nu), psi_2 = Nu_2(nu) = (nu-nu_r)/(h_nu). $

#let Psi(x) = $limits(psi _#x)^(->)$

Перейдём в трехмерный случай и базисные функции получат вид:

#grid(
  columns: (4cm, 4cm, 4cm, 4cm),
  gutter: 22pt,
  [$#Psi(1) = vec(Y_1(y) dot Z_1(z), 0, 0),$], [$#Psi(2) = vec(Y_2(y) dot Z_1(z), 0, 0),$], [$#Psi(3) = vec(0, X_1(x) dot Z_1(z), 0),$], [$#Psi(4) = vec(0, X_2(x) dot Z_1(z), 0),$],
  
  [$#Psi(5) = vec(0, 0, X_1(x) dot Y_1(y)),$], [$#Psi(6) = vec(0, 0, X_2(x) dot Y_1(y)),$], [$#Psi(7) = vec(0, 0, X_1(x) dot Y_2(y)),$], [$#Psi(8) = vec(0, 0, X_2(x) dot Y_2(y)),$],

  [$#Psi(9) = vec(Y_1(y) dot Z_2(z), 0, 0),$], [$#Psi(10) = vec(Y_2(y) dot Z_2(z), 0, 0),$], [$#Psi(11) = vec(0, X_1(x) dot Z_2(z), 0),$], [$#Psi(12) = vec(0, X_2(x) dot Z_2(z), 0).$]
)

Локальные матрицы жёсткости и масс собираются по формулам:
$ G_(i j) = limits(integral)_(Omega) "rot" #Psi("i") dot "rot" #Psi("j") d Omega, #h(0.5cm) M_(i j) =  limits(integral)_(Omega) gamma #Psi("i") dot #Psi("j") d Omega. $

Компоненты вектора b правой части определяются соотношением:

$ b_i = limits(integral)_(Omega) limits(F)^(->) dot #Psi("i") d Omega $

В программе используется численное интегрирование для вычисления локальных матриц. А вектор правой части вычисляется посредством умножения матрицы масс на вектор $limits(F)^(->)$.

= Входные файлы

*GridParameters* - файл для задания пространственной сетки. Файл состоит из 10 + n-строк (взависимости от количества заданных зон для коэффициента $sigma$), каждая строка состоит из значений, где через пробел разделяются значения (в примере будет показано с использованием символа *|* для более понятного представления читателю):

_1. начало по х | конец по х | кол-во шагов по х | коэф. разрядки по х
2. разбиения слоёв по x
3. начало по у | конец по у | кол-во шагов по у | коэф. разрядки по у
4. разбиения слоёв по y
5. начало по z | конец по z | кол-во шагов по z | коэф. разрядки по z
6. разбиения слоёв по z
7. левая граница | правая граница | нижняя граница | верхняя граница | задняя граница | передняя граница
8. $mu$ | $epsilon$
9. количество разрывов области n
10 - 10 + n. значение $sigma$ | индекс начала зоны по x | индекс конца зоны по x | индекс начала зоны по y | индекс конца зоны по y | индекс начала зоны по z | индекс конца зоны по z_

Пример задания области из будущего теста:

1. 0 5 4 1
2. 0 1 4 5
3. 0 5 4 1
4. 0 1 4 5
5. 0 5 4 1
6. 0 1 4 5
7. 1 1 1 1 1 1
8. 1 1
9. 4
10. 3 0 1 0 2 0 3
11. 1 0 1 2 3 0 3
12. 5 1 3 0 3 0 1
13. 3 1 3 0 3 1 3

*TimeGridParameters* - файл для задания временной сетки. Во временной сетке нужно всего 4 параметра:

_1. начало по t | конец по t | кол-во шагов по t | коэф. разрядки по t_

Пример:

1. 0 10 20 1.1

= Описание и возможности программы

1. В программе можно выбрать как будут собираться временные слои: физически (собирая из начальных условий двухслойную неявную схему, далее трёхслойную и потом четырёхслойную) или сгенерировать точные значения на всех трёх слоях.
2. СЛАУ можно решать с помощью BCGStab-LU или решение методом LU-разложения.

== Немного про классы в программе:

1. Основным классом является *FEM*, в котором происходит генерация портрета глобальной матрицы, сборка локальных матриц и векторов, учёт краевых условий.
2. Генерация пространственной сетки происходит в классе *Grid*, а временной в *TimeGrid*
3. Для реализации численного интегрирования были созданы 2 record struct-а: *SegmentGaussOrder9*, который является методом Гаусса-9, и *TriLinearVectorBasis* для удобного вычисления базисных вектор-функций.
4. Само численное интегрирование реализовано в классе *Integration*, в который передаются все квадратуры.
5. Оставшиеся классы реализуют структуры для удобной работы с основными классами. Например есть 2 класса *SparseMatrix* и *Vector*, которые используются для хранения глобальной матрицы и вектора соответственно.

#pagebreak()

= Пример работы программы

#v(0.2cm)

#grid(
  columns: (2cm, 3cm),
  gutter: 10pt,
  [$dot limits(A)^(->) = vec(y^2, x, z)$,], [$limits(F)^(->) = vec(-1, 0, 0)$,],

  [$dot "Пространственная сетка: " x, y, z in [0, 5], h_(x, y, z) = 1.25$], [], 
  
  [$dot sigma_1 = 3, x in [0, 1], y in [0, 4], z in [0, 5],$], [],

  [$dot sigma_2 = 1, x in [0, 1], y in [4, 5], z in [0, 5],$], [],

  [$dot sigma_3 = 5, x in [1, 5], y in [0, 5], z in [0, 1],$], [],

  [$dot sigma_4 = 3, x in [1, 5], y in [0, 5], z in [1, 5],$], [],

  [$dot mu = 1$,], [$epsilon = 1$,],
  

  [$dot "Временная сетка: " t in [0, 10], 20 "разбиений", k = 1.1$], [],

  [$dot "Краевые условия первого рода заданы на всех границах"$]
)

#v(0.2cm)

#let results = csv("csv/1.csv")

#show figure.where(
  kind: table 
): set figure.caption(position: top)


// ВСТАВКА CSV ФАЙЛА, НАДО В ФУНКЦИЮ ПЕРЕПИСАТЬ
#align(center)[
  #figure(
    supplement: none,
  table(
    columns: 2,
    fill: (_, row) => if calc.even(row) { luma(220) } else { white }, 
    [$t_i$], [Error],
    ..results.flatten(),
  ),
  caption: [Полученный результат])
]

#pagebreak()

= Исследования

== Порядок аппроксимации по пространству

=== Кубическая функция

#v(0.2cm)

#grid(
  columns: (2cm, 3cm),
  gutter: 10pt,
  [$dot limits(A)^(->) = vec(y^3, x^3, z^3)$,], [$limits(F)^(->) = vec(-6y, -6x, 0)$,],

  [$dot "Пространственная сетка: " x, y, z in [0, 5], h_(x, y, z) = 1.25$], [], 
  
  [$dot mu = 1$,], [],
  
  [$dot "Краевые условия первого рода заданы на всех границах"$]
)

#let results = csv("csv/2.csv")

#show figure.where(
  kind: table 
): set figure.caption(position: top)

#align(center)[
  #figure(
    supplement: none,
  table(
    columns: 4,
    fill: (_, row) => if calc.even(row) { luma(220) } else { white }, 
    [Edge number], [Approx value], [Real value], [Error],
    ..results.flatten(),
  ),
  caption: [Полученный результат])
]

#pagebreak()

=== Полином 4-ой степени

#v(0.2cm)

#grid(
  columns: (2cm, 3cm),
  gutter: 10pt,
  [$dot limits(A)^(->) = vec(y^4, x^3, z^3)$,], [$limits(F)^(->) = vec(-12y^2, -6x, 0)$,],

  [$dot "Пространственная сетка: " x, y, z in [0, 5], h_(x, y, z) = 1.25,$], [], 
  
  [$dot mu = 1$,], [],
  
  [$dot "Краевые условия первого рода заданы на всех границах"$]
)

#let results = csv("csv/3.csv")

#show figure.where(
  kind: table 
): set figure.caption(position: top)

#align(center)[
  #figure(
    supplement: none,
  table(
    columns: 4,
    fill: (_, row) => if calc.even(row) { luma(220) } else { white }, 
    [Edge number], [Approx value], [Real value], [Error],
    ..results.flatten(),
  ),
  caption: [Полученный результат])
]

#v(0.5cm)

Порядок аппроксимации = 3.

== Порядок сходимости по пространству

#v(0.2cm)

#grid(
  columns: (2cm, 3cm),
  gutter: 10pt,
  [$dot limits(A)^(->) = vec(y^4, x^3, z^3)$,], [$limits(F)^(->) = vec(-12y^2, -6x, 0)$,],

  [$dot "Пространственная сетка: " x, y, z in [0, 5], h_(x, y, z) = 1.25,$], [], 
  
  [$dot mu = 1$,], [],
  
  [$dot "Краевые условия первого рода заданы на всех границах"$]
)

#align(center)[
  #figure(
    supplement: none,
  table(
    columns: 3,
    fill: (_, row) => if calc.even(row) { luma(220) } else { white }, 
    [Шаг], [Погрешность], [Порядок сходимости],
    [$h$], [$0.07208790544253746$], [-],
    [$h/2$], [$0.01921482055319743$], [1.9],
    [$h/4$], [$0.004815259951729759$], [1.99]
  ),
  caption: [Полученный результат])
]

#v(0.05cm)

Порядок сходимости по пространству $approx 2$

#v(0.3cm)

== Порядок аппроксимации по времени

=== Кубическая функция

#v(0.2cm)

#grid(
  columns: (2cm, 3cm),
  gutter: 10pt,
  [$dot limits(A)^(->) = vec(t^2, t, t^3)$,], [$limits(F)^(->) = vec(2 sigma t + 2 epsilon, sigma, 3 sigma t^2 + 6 epsilon t)$,],

  [$dot "Пространственная сетка: " x, y, z in [0, 5], h_(x, y, z) = 1.25$], [], 
  
  [$dot sigma_1 = 3, x in [0, 1], y in [0, 4], z in [0, 5],$], [],

  [$dot sigma_2 = 1, x in [0, 1], y in [4, 5], z in [0, 5],$], [],

  [$dot sigma_3 = 5, x in [1, 5], y in [0, 5], z in [0, 1],$], [],

  [$dot sigma_4 = 3, x in [1, 5], y in [0, 5], z in [1, 5],$], [],

  [$dot mu = 1$,], [$epsilon = 1$,],
  

  [$dot "Временная сетка: " t in [0, 10], 20 "разбиений", h = 0.5$], [],

  [$dot "Краевые условия первого рода заданы на всех границах"$]
)

#v(0.2cm)

#let results = csv("csv/4.csv")

#show figure.where(
  kind: table 
): set figure.caption(position: top)


#align(center)[
  #figure(
    supplement: none,
  table(
    columns: 2,
    fill: (_, row) => if calc.even(row) { luma(220) } else { white }, 
    [$t_i$], [Error],
    ..results.flatten(),
  ),
  caption: [Полученный результат])
]

#pagebreak()

=== Полином 4-ой степени

#v(0.2cm)

#grid(
  columns: (2cm, 3cm),
  gutter: 10pt,
  [$dot limits(A)^(->) = vec(t^2, t, t^4)$,], [$limits(F)^(->) = vec(2 sigma t + 2 epsilon, sigma, 3 sigma t^2 + 6 epsilon t)$,],

  [$dot "Пространственная сетка: " x, y, z in [0, 5], h_(x, y, z) = 1.25$], [], 
  
  [$dot sigma_1 = 3, x in [0, 1], y in [0, 4], z in [0, 5],$], [],

  [$dot sigma_2 = 1, x in [0, 1], y in [4, 5], z in [0, 5],$], [],

  [$dot sigma_3 = 5, x in [1, 5], y in [0, 5], z in [0, 1],$], [],

  [$dot sigma_4 = 3, x in [1, 5], y in [0, 5], z in [1, 5],$], [],

  [$dot mu = 1$,], [$epsilon = 1$,],
  

  [$dot "Временная сетка: " t in [0, 10], 20 "разбиений", h = 0.5$], [],

  [$dot "Краевые условия первого рода заданы на всех границах"$]
)

#v(0.2cm)

#let results = csv("csv/5.csv")

#show figure.where(
  kind: table 
): set figure.caption(position: top)


#align(center)[
  #figure(
    supplement: none,
  table(
    columns: 2,
    fill: (_, row) => if calc.even(row) { luma(220) } else { white }, 
    [$t_i$], [Error],
    ..results.flatten(),
  ),
  caption: [Полученный результат])
]

#v(0.2cm)

Порядок аппроксимации = 4.

#pagebreak()

== Порядок сходимости по времени

#v(0.2cm)

#grid(
  columns: (2cm, 3cm),
  gutter: 10pt,
  [$dot limits(A)^(->) = vec(t^2, t, t^4)$,], [$limits(F)^(->) = vec(2 sigma t + 2 epsilon, sigma, 3 sigma t^2 + 6 epsilon t)$,],

  [$dot "Пространственная сетка: " x, y, z in [0, 5], h_(x, y, z) = 1.25$], [], 
  
  [$dot sigma_1 = 3, x in [0, 1], y in [0, 4], z in [0, 5],$], [],

  [$dot sigma_2 = 1, x in [0, 1], y in [4, 5], z in [0, 5],$], [],

  [$dot sigma_3 = 5, x in [1, 5], y in [0, 5], z in [0, 1],$], [],

  [$dot sigma_4 = 3, x in [1, 5], y in [0, 5], z in [1, 5],$], [],

  [$dot mu = 1$,], [$epsilon = 1$,],
  
  [$dot "Временная сетка: " t in [0, 10], 20 "разбиений", h = 0.5$], [],

  [$dot "Краевые условия первого рода заданы на всех границах"$]
)

#v(0.2cm)

#align(center)[
  #figure(
    supplement: none,
  table(
    columns: 3,
    fill: (_, row) => if calc.even(row) { luma(220) } else { white }, 
    [Шаг], [Погрешность], [Порядок сходимости],
    [$h$], [$0.17562942397686412$], [-],
    [$h/2$], [$0.013906410305277344$], [3.65],
    [$h/4$], [$0.001020206923985269$], [3.76]
  ),
  caption: [Полученный результат])
]

#v(0.2cm)

Порядок сходимости по времени $approx 4$.

#pagebreak()

= Листинг

== FEM.cs
#include-code("code/FEM.cs")

== SLAE.cs
#include-code("code/SLAE.cs")

== Grid.cs
#include-code("code/Grid.cs")

== Basis.cs
#include-code("code/Basis.cs")

== Integration.cs
#include-code("code/Integration.cs")

== QuadratureNode.cs
#include-code("code/QuadratureNode.cs")