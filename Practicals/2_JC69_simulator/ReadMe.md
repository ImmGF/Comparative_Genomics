Data:

-  S - input sequence
-  a - the alpha parameter of the JC69 model
-  d - time quantum
-  K - number of simulation steps

Write a simulator for the JC69 model that will generate the resulting sequence of S by applying K substitution steps, where one step is the transformation of the entire sequence over time d.

For each sequence created in the i-th step, calculate how the frequency distribution of nucleotides differs from the stationary r (Euclidean distance or city metric is enough). Draw a graph.

How will this graph change when we scale d to d / 2 and double the number of steps?
