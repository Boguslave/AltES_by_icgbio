## AltES_by_icgbio
# Installation:
Install Python 3.x and corresponding versions of NumPy, SciPy, Pandas, Matplotlib, Seaborn and itertools.

Download and install STAR (version 2.5 or later) for FASTQ input.

Download rMATS v4.0.2 (turbo)

Obtain STAR genome index for genome by either of the following two ways

Download and unpacking AltES

# Preprocessing data:

Download pre-built STAR indexes if using Human (hg38, hg19) or Mouse (mm10).

Build your own STAR index following STAR manual from genome fasta sequence

Untar rMATs and STAR indexes.

Run rMATs for STAR indixes (https://rmats.sourceforge.io/user_guide.htm)

# Using AltES:
```
python altes.py rmats.txt output.txt delta
```
rmats.txt&nbsp;&nbsp;&nbsp;&nbsp; 	            Directory and name of rMATs output file (expected SE.MATS.JC.txt)

output.txt&nbsp;&nbsp;&nbsp;&nbsp;              Directory and name of output text file

delta&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;                   Float number for accuracy (pvalue). Default expected 0.000001

# How this tool working:
Первый шаг обработки заключается в преобразовании данных RNA-seq в данные величины экспрессии каждого экзона с помощью пакета rMATs. Подробности по этому шагу можно найти в приложении 1. 

Далее, по данным rMATS для каждого гена компилируем матрицу экспрессии АС экзонов A (рис. 1) размерности    , где m – количество образцов (60), а n – количество АС экзонов, содержащую данные об экспрессии событий пропуска и вставки экзонов в каждом конкретном образце. 

<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\Large&space;A%20=%20\begin{pmatrix}a_{11}%20&%20a_{12}%20&%20\cdots%20&%20a_{1n}%20\\%20a_{21}%20&%20a_{22}%20&%20\cdots%20&%20a_{2n}%20\\%20\vdots%20&%20\vdots%20&%20\ddots%20&%20\vdots%20\\%20a_{m1}%20&%20a_{m2}%20&%20\cdots%20&%20a_{mn}%20\end{pmatrix}" width="40%">
</p>

Также введем матрицу S, содержащую данные об экспрессии пропуска экзонов

<p align="center">
<img src="https://latex.codecogs.com/svg.latex?\Large&space;S%20=%20\begin{pmatrix}s_{11}%20&%20s_{12}%20&%20\cdots%20&%20s_{1n}%20\\%20s_{21}%20&%20s_{22}%20&%20\cdots%20&%20s_{2n}%20\\%20\vdots%20&%20\vdots%20&%20\ddots%20&%20\vdots%20\\%20s_{m1}%20&%20s_{m2}%20&%20\cdots%20&%20s_{mn}%20\end{pmatrix}" width="40%">
</p>

Каждая строка соответствует АС экзону, каждый столбец - эксперименту. В столбцах расположены частоты (число) пропусков экзонов в каждом эксперименте.

После производится подсчет суммарной экспрессии включения и исключения по каждому образцу и производится фильтрация случаев, в которых некорректные результаты выдает метрика корреляции

Для каждой пары подсчитывается усредненная по весу корреляция Спирмена между строками включения и пропуска экзонов, предназначенная для поиска нелинейных зависимостей между переменными, представленных целыми числами.

Формат выдачи AltES
 
Столбцами таблицы являются:

•	GeneID по аннотации

•	GeneSymbol – краткое имя гена

•	id1, id2 – id экзона в выдаче rMATs

•	p-критерий для корреляции экзонов

•	Координаты экзонов

•	Их ψ
