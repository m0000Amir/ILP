# A Problem of Optimal Location of Given Set of Base Stations in Wireless Networks with Linear Topology
Here we will formulate our placement problem in the form of a mixed – integer programming model. 

**The problem description**.

We have a line segment of length `L` with the ends points `a_1` and `a_n`. Each point `a_i` ,`i = 1 ... n` is defined with its one-dimensional coordinate `l_i`. Also we have a set of stations as `S = {s_j}`, `j=1 ... m`. 

Each station has two parameters: a coverage radius `r_j` and a communication radius `R_j`. 
Then we can define a station `s_j` as a set of parameters: `s_j = {r_j, R_j\}`.

There are special stations `s_0` and `s_{n+1}` which are gateways. These stations are already placed at the ends `a_0` and `a_{n+1}`. 
For those stations  `r_0 = r_{n+1} = 0`. 

We need **to maximaze total coverage** at the line segment.
