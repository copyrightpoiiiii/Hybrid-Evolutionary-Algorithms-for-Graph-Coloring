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
But I am not sure whether this operation can accept a poor solution.<br>
<br>
upd3:<br>
Now, there are two main problems needs to be solved.<br>
First, the time complexity of LS operator is too high.<br>
Second, the algorithm to solve the distance between two answers may give an approximate solution.<br> 
<br>
upd4:<br>
Great optimization, faster running speed and faster convergence speed.<br>
<br>
upd5:<br>
In UPD_HGC.cpp, we use some different and better way to Local_Search.<br>
However, we should accelerate the convergence speed.<br>
