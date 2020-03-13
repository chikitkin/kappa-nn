KAPPA (Kinetic Approach to Physical Processes in Atmospheres)
=============================================================

Используемые сторонние библиотеки
-----------------------------------

1. [yaml-cpp](https://github.com/jbeder/yaml-cpp)
2. [Armadillo](http://arma.sourceforge.net/)

Инструкции по установке (Windows)
---------------------------------

Создаете пустой C++ проект в Visual Studio.

Armadillo и yaml-cpp: лежат на Гугл-диске (Проект/Библиотеки/libs.zip) (когда распакуете libs.zip в какую-нибудь папку, там будет две папки - Armadillo и yamllib).
Также скачиваете с Гугл-диска particles-v2.yaml (Проект/БД/particles-v2.yaml; только надо его переименовать в particles.yaml)

Настройки Visual Studio: все настройки, если отдельно не оговорено, для всех конфигураций сразу (и Debug, и Release)

	1. Открываете свойства вашего Solution, в разделе Configuration Properties (если Visual Studio русскоязычная, раздел "Свойства конфигурации" меню "Свойства"), в пункте "Платформа" выбираете x64
	2. В раздел "Header files" ("Заголовочные файлы") в дереве проекта добавить kappa.hpp (из папки Kappa/src)
	3. В раздел "Source files" ("Файлы исходного кода") в дереве проекта добавить все файлы из папок Kappa/src/approximations, Kappa/src/Particles, а также файлы interaction.cpp и numeric.cpp
	4. В разделе C/C++ в пункте Include Directories ("Дополнительные каталоги включаемых файлов") добавляете папки Kappa/src, Armadillo/include и yamllib/include
	5. В разделе C/C++ - Preprocessor Definitions (C/C++ - Препроцессор - Определения препроцессора) добавляете ARMA_USE_LAPACK и ARMA_USE_BLAS (каждое - на новой строке)
	6. В разделе Linker - General (Компоновщик - Общие) в меню Additional Library Directories (дополнительные каталоги библотек) добавляете папки Armadillo/examples/lib_win64 и yamllib
	7. В разделе Linker - Input (Компоновщик - Ввод) в меню Additional Dependencies (Дополнительные зависимости) добавляете lapack_win64_MT.lib и blas_win64_MT.lib (каждое - на новой строке)
	8. В разделе Linker - Input в меню Additional Dependencies добавляете libyaml-cppmdd.lib (для конфигурации Debug) и libyaml-cppmd.lib (для конфигурации Release)

После этого в раздел "Source files" можно добавить любой файл из папки Kappa/tests и компилировать проект.

Для работы программы также надо взять файлы lapack_win64_MT.dll и blas_win64_MT.dll из Armadillo/examples/lib_win64, скопировать их в директорию с запускаемым exe-файлом, а также скопировать в эту директорию файл particles.yaml и все должно заработать.


Тесты
-----

Kappa/tests/approximations - ZintTest.cpp (выводит стат. суммы, осредненную энергию и теплоемкости в однотемпературном приближении) и ZrotTest.cpp (выводит вращательную теплоемкость для моделей жесткого и нежесткого ротатора, а также по упрощенной формуле, для разных колебательных уровней).
