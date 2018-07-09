network.file <- "../data/network.txt"
expression <- "../data/expression.txt"
description.file <- "../data/description.csv"
metadata <- "../data/metadata.txt"

Necorr(network.file=network.file, expression=expression,
       description.file=description.file,
       condition="cond1",metadata=metadata, name="test", NSockets = 4)
