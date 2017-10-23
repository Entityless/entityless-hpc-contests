# UPDATE 20190514

This repo may contain lots of codes in bad style and toxic commit histors, written when I am young, simple and naive.

This repo is no longer maintained.

# dnn
代码主要基本框架和CPU矩阵乘法是[EnigmaHuang](https://github.com/EnigmaHuang)写的

ASC16 决赛赛题，BP神经网路在天河二号上的实现（CPU+MIC）。

天河二号（曾经）每个节点有3张MIC。

这份代码中，每个节点跑一个进程，每个进程offload到三张MIC上。

原版代码是单精度浮点的，这份代码是双精度的。因为我们遇到了精度问题，在跑决赛的算例的时候，单精度的代码那没办法满足三个算例的精度的需求。双精度的代码好像是能够达到精度要求的。反正我们没跑完，拿了0分。

主要的代码都在dnn_kernel.cpp里面。当时也是初生牛犊不怕虎，刚写一年多代码，就这样胡搞一通。从代码工程学的角度来看，这份代码简直是令人发指，真不知道当时自己怎么能够写下去的。

代码应用了多重流水线来掩盖通信的开销。同时，CPU也开了很多条线程，每条线程都有不同的作用。对于MPI_Allreduce，我采用了手动实现sendrecv和计算来代替，计算部分使用了多线程实现，并且实现了分段sendrecv和叠加（在sendrecv一段数据的同时，能够多线程sum）。worker_count指的是allreduce的worker数量。

| 线程号 | 作用 |
| --- | --- |
| 0 | offload到MIC:0，然后负责通信 |
| 1 | offload到MIC:1，然后负责MIC:0的PIC-e传输 |
| 2 | offload到MIC:2，然后负责MIC:1的PIC-e传输 |
| 3 | 负责MIC:2的PIC-e传输 |
| 4 | Master，总管所有flag |
| 5 ~ 5 + worker_count - 1 | 手动Allreduce的worker |
| 5 + worker_count | IO流水线 |
| 其余 | CPU协同计算（不建议） |


代码使用了一定的 dirty code 来实现CPU和MIC的同步。在整个程序的运行中，只会offload一次。CPU和MIC通过flag进行通信。CPU+MIC的问题在于，通信只能够从CPU这边发起。所以，MIC上面有一个字段是flag。CPU有一条线程，专门用来检测。当时我还菜到都不知道有volatile这个关键字，所以flag的用法极为奇葩。和flag相关的地方使用了GetSecondELement函数，能不能看懂就随缘吧。（反正我是不想再去review我大二下学期写的flag机制

MIC卡有一点非常坑，就是，如果有两条线程同时尝试和MIC进行数据交换，那么程序可能会挂掉。因此，本程序对于三张MIC，每张MIC由一条线程负责数据传输（包括in和out）。
