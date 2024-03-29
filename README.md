# Модель двумерного газа

Модель поведения частиц в двумерном газе с использованием методов молекулярной динамики. 

## Описание

Модель молекулярной динамики для симуляции двумерной системы частиц, взаимодействующих через потенциал [Леннарда-Джонса](https://en.wikipedia.org/wiki/Lennard-Jones_potential) в периодических граничных условиях. Для решения уравнений движения использовался [алгоритм
Стремера-Верле](https://en.wikipedia.org/wiki/Verlet_integration).

Реализованы:
- расчёт кинетической, потенциальной и полной энергии
- определение и контроль температуры системы
- расчёт среднеквадратичного смещения
- расчёт парной корреляционной функции
- сохранение позиций частиц и наблюдаемых в файл
- визуализация наблюдаемых, всей системы и траектории отдельной частицы

Подробное описание методов можно найти в работе [Vollmayr-Lee, Katharina, Introduction to Molecular Dynamics Simulations.](https://www.researchgate.net/publication/318567658_Introduction_to_Molecular_Dynamics_Simulation).


![Main window](https://github.com/Onishenko-sci/Molecular-Dinamics-Gas/blob/main/Window.png)


## Использование


### 0. Зависимости
Убедитесь, что у вас установлены необходимые зависимости, такие как g++ для компиляции C++ кода и Python для выполнения скрипта визуализации. Установите необходимые для визуализации библиотеки:

```bash
pip install numpy matplotlib tk scipy
```
    
### 1. Клонирование репозитория

1. Откройте терминал и выполните следующую команду для клонирования репозитория:

    ```bash
    git clone https://github.com/Onishenko-sci/Molecular-Dinamics-Gas.git
    ```


2. Перейдите в каталог проекта:

    ```bash
    cd Molecular-Dinamics-Gas
    ```

### 2. Компиляция и Запуск

Для быстрого старта, выполните следующие шаги:

1. Перейдите в каталог с Makefile:

    ```bash
    cd source
    ```
    
2. Для изменния параметров модели внесите изменения в ```main.cpp``` файл проекта.

3. Используйте команду make с ключом all:

    ```bash
    make all
    ```

Эта команда скомпилирует исходный код программы и выполнит программу для симуляции молекулярной динамики, а затем отобразит результаты с использованием скрипта Python. 

### 3. Просмотр Результатов

Если у вас уже есть скомпилированный файл и вы хотите только отобразить результаты, выполните:

```bash
make show
```

Это выполнит только  Python скрипт для визуализации.

