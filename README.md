# Домашнее задание №4: анализ Hi-C данных с помощью cooler и cooltools

В данной работе был выполнен анализ Hi-C данных в формате `.mcool`. Для анализа были выбраны два файла с портала 4DN. Работа выполнялась с помощью пакетов `cooler` и `cooltools`, которые позволяют получать информацию о Hi-C матрицах, извлекать таблицы бинов и контактов, работать со сбалансированными матрицами, строить зависимость числа контактов от геномного расстояния, а также рассчитывать insulation score и находить границы TAD.

В ходе работы были выполнены следующие этапы:

1. получение информации и атрибутов Hi-C матриц;
2. извлечение таблиц бинов и контактов;
3. проверка различий между сырыми и сбалансированными контактами;
4. построение зависимости контактной частоты от расстояния;
5. расчет insulation score для выбранного участка;
6. определение границ TAD;
7. сравнение результатов для двух Hi-C файлов;
8. анализ зависимости числа найденных TAD boundaries от размера окна.


Скачиваю два .mcool файла: 4DNFIQWCSCVX и 4DNFIGATB8OL

Для начала были просмотрены доступные разрешения внутри .mcool файлов с помощью команды cooler ls.

как видим у обоих доступны вплоть до:
data/sample1.mcool::/resolutions/10000000
data/sample2.mcool::/resolutions/10000000

Для анализа был выбран участок chr2:10000000-11000000 длиной 1 Mb. Анализ выполнялся на разрешении 10000 bp. Также были выбраны несколько размеров окна для расчета insulation score: 50 kb, 100 kb, 150 kb, 200 kb и 300 kb.

Получаем вот такой результат:
```
{
    "bin-size": 10000,
    "bin-type": "fixed",
    "creation-date": "2024-12-06T12:12:50.403600",
    "format": "HDF5::Cooler",
    "format-url": "https://github.com/open2c/cooler",
    "format-version": 3,
    "generated-by": "cooler-0.9.3",
    "genome-assembly": "unknown",
    "nbins": 308837,
    "nchroms": 24,
    "nnz": 285987146,
    "storage-mode": "symmetric-upper",
    "sum": 604633594
}
{
    "bin-size": 10000,
    "bin-type": "fixed",
    "creation-date": "2024-03-02T09:47:35.430518",
    "format": "HDF5::Cooler",
    "format-url": "https://github.com/mirnylab/cooler",
    "format-version": 3,
    "generated-by": "cooler-0.8.3",
    "genome-assembly": "unknown",
    "nbins": 272564,
    "nchroms": 21,
    "nnz": 128699895,
    "storage-mode": "symmetric-upper",
    "sum": 212791328
}
```

С помощью команды cooler info была получена основная информация о Hi-C матрицах: размер бина, число бинов, число ненулевых пикселей, сборка генома и другие параметры. Команда cooler attrs использовалась для просмотра дополнительных атрибутов cooler-объектов.

```bash
cooler dump -t bins -H "$URI1" | head -n 20 > results/tables/sample1_bins_head.tsv
cooler dump -t bins -H "$URI2" | head -n 20 > results/tables/sample2_bins_head.tsv
```

С помощью cooler dump была получена таблица бинов. В таблице присутствуют столбцы chrom, start, end и weight. Столбцы chrom, start и end задают координаты бина в геноме, а weight содержит веса, которые используются для балансировки Hi-C матрицы.

```bash
cooler dump --join -H -r "$REGION" "$URI1" | head -n 50 > results/tables/sample1_raw_contacts_head.tsv
cooler dump --join -H -r "$REGION" "$URI2" | head -n 50 > results/tables/sample2_raw_contacts_head.tsv
```

Была получена таблица контактов для выбранного участка. В сырой таблице контактов основным численным столбцом является count. Этот столбец содержит исходное число контактов между двумя бинами, поэтому такие значения не являются сбалансированными.

```bash
cooler dump --join --balanced -H -r "$REGION" "$URI1" | head -n 50 > results/tables/sample1_balanced_contacts_head.tsv
cooler dump --join --balanced -H -r "$REGION" "$URI2" | head -n 50 > results/tables/sample2_balanced_contacts_head.tsv
```

Была получена таблица сбалансированных контактов. В этой таблице появляется столбец balanced. Он содержит значения контактов, нормированные с учетом весов бинов. Поэтому таблица, полученная с параметром --balanced, является сбалансированной.
```bash
source config.sh
export S1 S2 MC1 MC2 RES URI1 URI2 CHR START END REGION THREADS WINDOWS MAIN_WINDOW
python scripts/open_balanced_matrix.py
cat results/info/sample1_balanced_matrix_info.txt
cat results/info/sample2_balanced_matrix_info.txt
```

Hi-C матрицы были открыты с помощью Python-интерфейса cooler. При использовании balance=True извлекается сбалансированная матрица контактов, так как к значениям контактов применяются веса бинов. Для выбранной хромосомы были получены матрицы внутрихромосомных контактов, а для выбранного участка длиной 1 Mb были сохранены отдельные сбалансированные матрицы.
Был создан view-файл для выбранной хромосомы. Этот файл задает область генома, в пределах которой будут выполняться расчеты cooltools. В данном случае анализ ограничен выбранной хромосомой.

```bash
cooltools expected-cis \
  --nproc "$THREADS" \
  --view data/view_chr.tsv \
  -o results/expected/sample1_expected_cis.tsv \
  "$URI1"

cooltools expected-cis \
  --nproc "$THREADS" \
  --view data/view_chr.tsv \
  -o results/expected/sample2_expected_cis.tsv \
  "$URI2"

head results/expected/sample1_expected_cis.tsv
head results/expected/sample2_expected_cis.tsv
```

Был построен график зависимости средней контактной частоты от геномного расстояния в логарифмических координатах. На графике видно, что контактная частота уменьшается с увеличением расстояния. Это ожидаемо для Hi-C данных, потому что близко расположенные участки хромосомы обычно контактируют чаще, чем удаленные.

#image[figures/contacts_vs_distance.png]

```bash
cooltools insulation \
  --threshold Li \
  --nproc "$THREADS" \
  --view data/view_chr.tsv \
  -o results/insulation/sample1_insulation.tsv \
  "$URI1" $WINDOWS

cooltools insulation \
  --threshold Li \
  --nproc "$THREADS" \
  --view data/view_chr.tsv \
  -o results/insulation/sample2_insulation.tsv \
  "$URI2" $WINDOWS

head results/insulation/sample1_insulation.tsv
head results/insulation/sample2_insulation.tsv
```

С помощью cooltools insulation был рассчитан insulation score для нескольких размеров окна: 50 kb, 100 kb, 150 kb, 200 kb и 300 kb. Также были определены границы TAD. В таблице результатов присутствуют значения insulation score, сила границы boundary_strength и логический столбец is_boundary, показывающий, является ли данный бин границей TAD.

На основе результатов cooltools insulation были созданы два BED-файла с границами TAD. Для основного анализа использовалось окно 100 kb. В BED-файле четыре столбца: chrom, start, end и score. В поле score записана сила границы TAD, то есть значение boundary_strength_100000.

Был построен график insulation score для выбранного участка длиной 1 Mb. Вертикальными пунктирными линиями отмечены найденные границы TAD. Минимумы insulation score соответствуют позициям, которые разделяют соседние домены. При сравнении двух образцов можно увидеть, совпадают ли положения границ TAD и насколько похожа локальная доменная организация участка.

Было подсчитано число найденных границ TAD для разных размеров окна: 50 kb, 100 kb, 150 kb, 200 kb и 300 kb. Количество найденных границ зависит от размера окна. Малые окна чувствительнее к локальным изменениям Hi-C карты и могут находить больше мелких границ. Большие окна сглаживают локальные особенности и отражают более крупномасштабную доменную структуру.

```bash
bedtools window \
  -w "$RES" \
  -a beds/sample1_TAD_boundaries_w100000.bed \
  -b beds/sample2_TAD_boundaries_w100000.bed \
  > results/tables/common_boundaries_sample1_sample2_w100000.tsv
```

Границы TAD между двумя файлами были сравнены с помощью bedtools window. Границы считались совпадающими, если расстояние между ними не превышало одно разрешение матрицы, то есть 10 kb. Совпадающие границы указывают на сходство доменной организации в выбранном участке, а различающиеся границы могут быть связаны с биологическими различиями между образцами или техническими особенностями данных.

В результате анализа были получены таблицы с информацией о матрицах, таблицы бинов и контактов, графики зависимости контактов от расстояния, графики insulation score, таблицы с числом границ TAD при разных размерах окна и BED-файлы с границами TAD.