%% 碳排放流 Carbon Emission Flow (CEF)

% 复现论文《电力系统碳排放流的计算方法初探》[1]。
% [1] 周天睿,康重庆,徐乾耀,等.电力系统碳排放流的计算方法初探[J].电力系统自动化,2012,36(11):44-49.

% 作者：罗清局
% 邮箱：luoqingju@qq.com
% 华南理工大学电力学院
% 综合智慧能源系统优化运行与控制团队 ISESOOC 爱思科

% 依赖库 MATPOWER https://matpower.org/

% 作者水平有限，难免会有错误及不足之处，敬请读者批评指正！

clc
clear

define_constants; % MATPOWER：Defines useful constants for indexing data

mpc = case14;

Pg = [120 40 60 19 20]'; % 设置发电机有功出力
mpc.gen(:, PG) = Pg; % 修改MATPOWER算例的发电机有功出力

mpopt = mpoption('verbose', 0, 'out.all', 0); % MATPOWER：不打印计算结果
res = rundcpf(mpc, mpopt); % MATPOWER：计算直流潮流
if res.success ~= 1
    error('----------直流潮流计算失败！----------')
end

N = size(res.bus, 1); % 节点数（母线数）
K = size(res.gen, 1); % 发电机数

Pd = res.bus(:, PD); % 节点负荷
gen_bus = res.gen(:, GEN_BUS); % 发电机节点

fbus = res.branch(:, F_BUS); % 线路 "from" 端节点
tbus = res.branch(:, T_BUS); % 线路 "to" 端节点

Pl_from = res.branch(:, PF); % 线路 "from" 端功率
Pl_to = res.branch(:, PT); % 线路 "to" 端功率

Pl_from(Pl_from < 0) = 0; % 反向的功率置零
Pl_to(Pl_to < 0) = 0; % 反向的功率置零

idx_PF = Pl_from > 0; % 线路 "from" 端功率索引
PB_F_Mat = sparse(fbus(idx_PF), tbus(idx_PF), Pl_from(idx_PF), N, N); % 线路 "from" 端潮流分布矩阵

idx_PT = Pl_to > 0; % 线路 "to" 端功率索引
PB_T_Mat = sparse(tbus(idx_PT), fbus(idx_PT), Pl_to(idx_PT), N, N); % 线路 "to" 端潮流分布矩阵

% 支路潮流分布矩阵(branch power flow distribution matrix) N 阶方阵
PB_Mat = PB_F_Mat + PB_T_Mat;

% 机组注入分布矩阵(power injection distribution matrix) K×N 阶矩阵
PG_Mat = sparse(1:K, gen_bus, Pg, K, N);

% 负荷分布矩阵(load distribution matrix) M×N 阶矩阵 M 为负荷数
% 为了简化，假设每个节点都存在负荷，不存在负荷的节点按照负荷为零处理
% 所以，负荷分布矩阵变为 N×N 阶矩阵
PL_Mat = sparse(1:N, 1:N, Pd, N, N);

PZ_Mat = [PB_Mat; PG_Mat];

% 节点有功通量矩阵(nodal active power flux matrix) N 阶对角阵
PN_Mat = sparse(1:N, 1:N, ones(1, N+K)*PZ_Mat, N, N); % 论文中的公式2

% 发电机组碳排放强度向量(unit carbon emission intensity vector)
EG_Vec = [875 525 0 520 0]'; % % 论文中的公式14

% 节点碳势向量(nodal carbon intensity vector)
EN_Vec = (PN_Mat - PB_Mat') \ (PG_Mat' * EG_Vec);  % 论文中的公式13

% 支路碳流率分布矩阵 (branch carbon emission flow rate distribution matrix) N 阶方阵
% 论文中的公式5可能存在笔误，应该是 RB = diag(EN)*PB
RB_Mat = sparse(1:N, 1:N, EN_Vec, N, N) * PB_Mat;
RB_Mat = RB_Mat./1000; % kgCO2/h ==> tCO2/h

% 负荷碳流率向量(load carbon emission rate vector)
RL_Vec = PL_Mat * EN_Vec; % 论文中的公式7
RL_Vec = RL_Vec./1000; % kgCO2/h ==> tCO2/h

% 支路碳流密度(branch carbon emission flow intensity)
EB_Mat = sparse(1:N, 1:N, EN_Vec, N, N) * spones(PB_Mat);

% 机组注入碳流率
IN_Vec = PG_Mat' * EG_Vec;
IN_Vec = IN_Vec./1000; % kgCO2/h ==> tCO2/h

%% 整理结果

L = size(res.branch, 1); % 线路数
Pl = res.branch(:, PF); % 线路功率
EB_Vec = zeros(L, 1); % 支路碳流密度
RB_Vec = zeros(L, 1); % 支路碳流率
for i = 1:L
    if Pl(i) > 0
        EB_Vec(i) = EB_Mat(fbus(i), tbus(i));
        RB_Vec(i) = RB_Mat(fbus(i), tbus(i));
    else
        EB_Vec(i) = EB_Mat(tbus(i), fbus(i));
        RB_Vec(i) = RB_Mat(tbus(i), fbus(i));
    end
end

Table2 = [(1:N)', full(diag(PN_Mat)), EN_Vec]; % 论文中的表2
Table3 = [fbus, tbus, Pl, EB_Vec, RB_Vec]; % 论文中的表3
Table4 = [(1:N)', RL_Vec, IN_Vec]; % 论文中的表4

fprintf('\n')
disp('表2 节点有功通量与节点碳势')
disp('节点    节点有功通量(MW)    节点碳势(gCO2/kWh)')
disp(Table2)
fprintf('\n')

fprintf('\n')
disp('表3 支路有功潮流与碳流率')
disp('起始节点    终止节点    支路有功潮流(MW)    支路碳流密度(gCO2/kWh)    碳流率(tCO2/h)')
disp(Table3)
fprintf('\n')

fprintf('\n')
disp('表4 负荷碳流率和机组注入碳流率')
disp('节点    负荷碳流率(tCO2/h)    机组注入碳流率(tCO2/h)')
disp(Table4)
fprintf('\n')