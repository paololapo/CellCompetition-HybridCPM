# Notes about the implementation
## Cell types
The simulation is implemented such that the cells are organized into 11 types, using type `0` for the medium. <br>
The original paper studies the competition between *wild-type* and *scribble* MDCK cells: 
* *Wild cells:* cells of type `[1, 3, 5, 7]` refer to the wild-type phenotype (winners). Type `9` refers to apoptotic wild-type cells.
* *Scribble cells:* cell of type `[2, 4, 6, 8]` refer to the cells depleted for the polarity protein scribble (losers). Type `10` refers to apoptotic scribble cells.

The different types within the same class represent the reflect progression through the cell cycle. Indeed, each type is updated during mitosis, e.g. 
```python
if self.parentCell.type == 1:
    self.childCell.type = 3
    self.parentCell.type = 3
```
Note that this is organized into a ring structure, i.e.
```python
elif self.parentCell.type == 7:
    self.childCell.type = 1  
    self.parentCell.type = 1
```