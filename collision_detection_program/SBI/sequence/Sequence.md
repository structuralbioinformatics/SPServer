#SEQUENCE LIBRARY

##SEQUENCE
The **SEQUENCE** object is composed basically by an _ID_ and the _SEQUENCE_.

The _str()_ function returns a tab format. Other formats should be called through the _format()_ function

###ATTRIBUTES
####id:
Corresponds to the first section of a space/tab separated given data, the rest of that data will be
stored, if required, in a _info_ attribute.
####sequence:
Sequence, as **string**
####info:
Corresponds to all the information given after the first space/tab
####is_gapped:
**Boolean** to check if the sequence contains gaps (gaps are understood as '-' or 'x')

###ESPECIAL PROPERTIES: CONTAINER PROPERTIES
The **SEQUENCE** object can be accesses as a container. That means that it can behave as an _iterable_. 
That behaviour is applied over the _sequence_ atribute. Thus:
```python

	a = Sequence(seqID='prot1',sequence='ARTGMP')
	print a[0]  # A
	for aa in a:
		print aa
		if aa == 'T': break
	# ART

	print len(a) # 6

```

###ESPECIAL PROPERTIES: COMPARISSON
**SEQUENCE(s)** can compare against each other in the following fashion:
```python

	a = Sequence(seqID='prot1',sequence='ARTGMP')
	b = Sequence(seqID='prot2',sequence='ARTGMP')
	c = Sequence(seqID='prot3',sequence='MKLARTGMPTYP')
	a == b: True
	a > b:  len(a) > len(b): False
	a < c:  len(a) < len(c): True
	a.contained(c):    True
	a.contains(b):     True
	a.contains('RTG'): True

```
All the _rich comparison methods_ are implemented

###FUNCTIONS


