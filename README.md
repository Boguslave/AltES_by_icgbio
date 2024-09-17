## AltES_by_icgbio
# Installation
Install Python 3.x and corresponding versions of NumPy, SciPy, Pandas, Matplotlib, Seaborn and itertools.
Download and install STAR (version 2.5 or later) for FASTQ input.
Download rMATS v4.0.2 (turbo)
Obtain STAR genome index for genome by either of the following two ways
Download pre-built STAR indexes if using Human (hg38, hg19) or Mouse (mm10).
Build your own STAR index following STAR manual from genome fasta sequence
Untar rMATS and STAR indexes.
# Using AltES:
```
python altes.py rmats.txt output.txt delta
```
rmats.txt 	            Directory and name of rMATs output file (expected SE.MATS.JC.txt)
output.txt              Directory and name of output text file
delta                   Float number for accuracy (pvalue). Default expected 0.000001

# How this tool working
Первый шаг обработки заключается в преобразовании данных RNA-seq в данные величины экспрессии каждого экзона с помощью пакета rMATs. Подробности по этому шагу можно найти в приложении 1. 
Далее, по данным rMATS для каждого гена компилируем матрицу экспрессии АС экзонов A (рис. 1) размерности    , где m – количество образцов (60), а n – количество АС экзонов, содержащую данные об экспрессии событий пропуска и вставки экзонов в каждом конкретном образце. 

<img src="https://latex.codecogs.com/svg.latex?\Large&space;A = \begin{pmatrix}a_{11} & a_{12} & \cdots & a_{1n} \\ a_{21} & a_{22} & \cdots & a_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ a_{m1} & a_{m2} & \cdots & a_{mn} \end{pmatrix}
" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />

Также введем матрицу S, содержащую данные об экспрессии пропуска экзонов

<img src="https://latex.codecogs.com/svg.latex?\Large&space;S = \begin{pmatrix}s_{11} & s_{12} & \cdots & s_{1n} \\ s_{21} & s_{22} & \cdots & s_{2n} \\ \vdots & \vdots & \ddots & \vdots \\ s_{m1} & s_{m2} & \cdots & s_{mn} \end{pmatrix}
" title="\Large x=\frac{-b\pm\sqrt{b^2-4ac}}{2a}" />

Каждая строка соответствует АС экзону, каждый столбец - эксперименту. В столбцах расположены частоты (число) пропусков экзонов в каждом эксперименте.
Для каждой пары подсчитывается усредненная по весу корреляция Спирмена между строками включения и пропуска экзонов, предназначенная для поиска нелинейных зависимостей между переменными, представленных целыми числами.

Формат выдачи AltES
 
Столбцами таблицы являются:
•	GeneID по аннотации
•	GeneSymbol – краткое имя гена
•	id1, id2 – id экзона в выдаче rMATs
•	p-критерий для корреляции экзонов
•	Координаты экзонов
•	Их ψ
