13期华为算法挑战赛 多核信息处理系统 multi-core message processing system

提供一个 generator 的链接(实际上是来自 algotester 某场比赛的原题)

主要想法是 先贪心地提高 affinity score, 然后模拟退火/爬山算法优化

模拟退火的实现: 任意取出一个 user, 调换 task 的顺序, 优先保证亲和性, 并计算得分

这里采用了链表的实现, 同时采用了复杂度优化, 把 msgtype 相同的合并成一个节点, 加快链表遍历速度

在预处理上, 考虑资源的紧缺程度 urgent (在 generator 代码中有体现), 范围从 0.0 - 3.0

如果 urgent > 1.2, 说明资源充足, 可以舍弃一些 core (牺牲部分 capability score) 换取更高的 affinity score

如果 urgent < 1.0, 说明资源紧缺, 可以舍弃一些 task, 优先保证其余的 task 能满足, 以换取更高的 capability score

具体实现, 可以拟合一个得分函数, 对每个不同数量的 core 计算得分, 选取得分最大的一个 amount.
