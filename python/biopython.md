# Sequence objects
Sequences are essentially strings of letters like `AGTACACTGGT`, which seems very natural since this is the most common way that sequences are seen in biological file formats.

The most important difference between `Seq` objects and standard Python strings is they have different methods. Although the `Seq` object supports many of the same methods as a plain string, its `translate()` method differs by doing biological translation, and there are also additional biologically relevant methods like `reverse_complement()`.

## Sequences act like strings
In most ways, we can deal with `Seq` objects as if they were normal Python strings, for example getting the length, or iterating over the elements:
```python
from Bio.Seq import Seq
my_seq = Seq("GATCG")
for index, letter in enumerate(my_seq):
    print("%i %s" % (index, letter))
```
output:
```
0 G
1 A
2 T
3 C
4 G
```

```python
print(len(my_seq))
```
output:
```
5
```

To access elements of the sequence in the same way as for strings:
(remember that Python counts from zero)
```python
>>> print(my_seq[0]) # first letter
G
>>> print(my_seq[2]) # third letter
T
>>> print(my_seq[-1]) # last letter
G
```

The `Seq` object has a `.count()` method, just like a string. Note that this means that like a Python string, this gives a non-overlapping count:
```python
>>> from Bio.Seq import Seq
>>> "AAAA".count("AA")
2
>>> Seq("AAAA").count("AA")
2
>>> "TTTTTT".count("T")
6
>>> Seq("TTTTTT").count("TT")
3
```

For some biological uses, you may actually want an overlapping count. When searching for single letters, this makes no difference:
```python
>>> from Bio.Seq import Seq
>>> my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
>>> len(my_seq)
32
>>> my_seq.count("G")
9
>>> 100 * (my_seq.count("G") + my_seq.count("C")) / len(my_seq)
46.875
```

While you could use the above snippet of code to calculate a GC%, note that the `Bio.SeqUtils` module has several GC functions already built. For example:
```python
>>> from Bio.Seq import Seq
>>> from Bio.SeqUtils import gc_fraction
>>> my_seq = Seq("GATCGATGGGCCTATATAGGATCGAAAATCGC")
>>> gc_fraction(my_seq)
0.46875
```

Note that using the `Bio.SeqUtils.gc_fraction()` function should automatically cope with mixed case sequences and the ambiguous nucleotide S which means G or C.

Also note that just like a normal Python string, the Seq object is in some ways “read-only”. If you need to edit your sequence, for example simulating a point mutation, look at `MutableSeq` object.


## Slicing a sequence



source: https://biopython.org/docs/latest/Tutorial/index.html

