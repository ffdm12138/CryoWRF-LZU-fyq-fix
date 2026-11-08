# CryoWRF Bug 修复案卷

这里汇总本地 CryoWRF 耦合稳定性工作的公开中文案卷。编号沿用既有案卷。详细条目按现象、定位、修复和验收范围记录。源码与公开发布校验见[发布回执](发布回执.md)。

| 编号 | 问题 | 证据状态 |
| --- | --- | --- |
| 1 | 零质量雪层触发应力断言 | 已修复。 |
| 2 | 冰面早退跳过 2 米诊断 | 已修复，适用范围限雪面接口。 |
| 3 | MYNN 稳定度查表边界及非有限输入 | 防御性修复和原风险窗穿越已有证据。 |
| 4 | 海侧温湿度诊断未传回地表驱动 | 已修复并核对接口。 |
| 5 | 海冰混合消费未产生的 `QZ0_SEA` | 已修复并完成本地验收。 |
| 6 | [续算状态不完整](6.cryowrf_restart_sn_bug.md) | 三组同构建通过。历史轨迹复现失败另列。 |
| 7 | 嵌套启动中的 SNOWPACK 有序格点布局与生命周期 | 已在本地复现并修复。 |
| 8 | 自适应时间步下 SNOWPACK 时钟不同步 | 原修复在有限窗口验证。 |
| 9 | [耦合积分中的非有限值](9.cryowrf_mynn_nan_instrumentation_report.md) | 局部重放有证据，长期闭环开放。 |
| 10 | [续算雪层容量截断](10.cryowrf_wrfrst_layer_truncation_bug_report.md) | 写端拒绝超容量，随本分支同步。 |
| 11 | [高海拔雪层增殖](11.cryowrf_snowpack_layer_proliferation_bug_report.md) | 层预算调查开放。 |
| 12 | [嵌套插值发送缓冲不足](12.cryowrf_rsl_nest_buffer_budget_bug_report.md) | 按需扩容，随本分支同步。 |
| 13 | [吹雪微物理写回非有限值](13.cryowrf_nest_bdy_v_nan_bug_report.md) | 写回修复随本分支同步，全年穿越独立验证。 |
| 14 | [浅雪吹雪恒零](14.cryowrf_zero_blowing_snow.md) | 物理模式十日窗口通过。公开分支的前置差异见回执。 |

[官方简化冷初始化](官方简化冷初始化.md)记录初始化修订及已验收边界。
