# WFA

An `offset` is the number of steps made along a diagonal
`(v,h)` is a position on the matrix

```
v = offset - k
h = offset
A_k = (text_len-query_len) The central diagonal
```

DP Matrix

```
              <--------- h --------->
     |---|---|---|---|---|---|---|---|
     |   |   | 0 | 1 | 2 | 3 | 4 | 5 |
     |---|---|---|---|---|---|---|---|
     |   |   | G | A | G | A | T | A |
     |---|---|---|---|---|---|---|---|
 ^   | 0 | G |   |   |   |   |   |   |
 |   |---|---|---|---|---|---|---|---|
 |   | 1 | A |   |   |   |   |   |   |
 |   |---|---|---|---|---|---|---|---|
 V   | 2 | C |   |   |   |   |   |   |
 |   |---|---|---|---|---|---|---|---|
 |   | 3 | A |   |   |   |   |   |   |
 |   |---|---|---|---|---|---|---|---|
 |   | 4 | C |   |   |   |   |   |   |
 |   |---|---|---|---|---|---|---|---|
 v   | 5 | A |   |   |   |   |   |   |
     |---|---|---|---|---|---|---|---|
```


### Citation

**Santiago Marco-Sola, Juan Carlos Moure, Miquel Moreto, Antonio Espinosa**. ["Fast gap-affine pairwise alignment using the wavefront algorithm."](https://doi.org/10.1093/bioinformatics/btaa777) Bioinformatics, 2020.
