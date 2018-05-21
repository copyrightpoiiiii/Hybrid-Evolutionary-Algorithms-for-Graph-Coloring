# Hybrid Evolutionary Algorithms for Graph Coloring
Based on this paper, we can implement an algorithm to solve the Graph Coloring Problem.<br>
First, we can enumerate an answer by dichotomy.<br>
Then, we use the hybird algorithm to check the answer.<br>
In this algorithm, we initialize a population and generate a offspring.<br>
Then we can use tabu search to optimize it.<br>
<br>
However, debugging this code is too hard for me. <br>
I hope you can help me.<br>
<br>
upd1:<br>
Now, the code can run normally.<br>
However, I found the LS operator had no obvious effect on the optimization of the answer.<br>
<br>
upd2:<br>
Through some tiny changes, the LS operator can make some remarkable optimizations.<br>
But I am not sure whether this operation can accept a poor solution.
