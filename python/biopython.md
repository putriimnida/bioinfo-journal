# Sequence objects
Sequences are essentially strings of letters like `AGTACACTGGT`, which seems very natural since this is the most common way that sequences are seen in biological file formats.

The most important difference between `Seq` objects and standard Python strings is they have different methods. Although the `Seq` object supports many of the same methods as a plain string, its `translate()` method differs by doing biological translation, and there are also additional biologically relevant methods like `reverse_complement()`.

## Sequences acts like strings
In most ways, we can deal with Seq objects as if they were normal Python strings, for example getting the length, or iterating over the elements:
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



