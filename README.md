# AliceFlow-2D

AliceFlow 2D - двумерный гидродинамический решатель в переменных скорость, давление, температура.
..* Произвольная расчётная область с вырезами на плоскости.
..* Прямоугольная структурированная расчётная сетка. Криволинейная граница расчётной области аппроксимируется ступеньками.
..* Только стационарный решатель.
..* Противопоточная схема для конвекции. 
..* Монотонизатор Рхи и Чоу для давления 1983.
..* SIMPLE алгоритм для связи уравнений неразрывность и скорость-давление 1972.
..* Приближение Обербека-Буссинеска для конвекции.
..* Программа полностью распараллелена на два потока центрального процессора.
..* ICCG решатель (метод сопряжённыъ градиентов предобусловленный неполным разложением Холецкого с нулевым заполнением) для поправки давления ван дер Ворста и Мейджеринка.


Примеры задач рассчитанных в AliceFlow 2D
![alt_text](https://github.com/kirill7785/AliceFlow-2D/blob/main/pic/Скорость%202021.jpg)
Cкорость, м/с

![alt_text](https://github.com/kirill7785/AliceFlow-2D/blob/main/pic/Давление%202021.jpg)
Давление, Па

![alt_text](https://raw.githubusercontent.com/kirill7785/AliceFlow-2D/main/pic/Скорость%20в%20сужающемся%20канале.bmp)
Cкорость, м/с

![alt_text](https://raw.githubusercontent.com/kirill7785/AliceFlow-2D/main/pic/Давление%20в%20сужающемся%20канале.bmp)
Давление, Па

![alt_text](https://raw.githubusercontent.com/kirill7785/AliceFlow-2D/main/pic/змеевик%20скорость.bmp)
Cкорость, м/с

![alt_text](https://raw.githubusercontent.com/kirill7785/AliceFlow-2D/main/pic/змеевик%20Давление.bmp)
Давление, Па

Выращивание кристаллов методом Бриджмена
![alt_text](https://raw.githubusercontent.com/kirill7785/AliceFlow-2D/main/pic/Скорость%20Бриджмен%202.bmp)
Cкорость, м/с
![alt_text](https://raw.githubusercontent.com/kirill7785/AliceFlow-2D/main/pic/Давление%20Бриджмен%202.bmp)
Давление, Па
![alt_text](https://raw.githubusercontent.com/kirill7785/AliceFlow-2D/main/pic/Температура%20Бриджмен%202.bmp)
Температура, С
