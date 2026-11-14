# VCF 变异工具集（FastAPI）

两个相互独立、共用同一套校验规则的流程：

1. **双 VCF 变异核对**：同一参考、同一样本下，核对两套变异识别程序输出的单样本 VCF，
   按**归一化后的等位变异**列出两侧一致、仅单侧出现与剂量不符的结果，并保留每条结论的
   VCF 原始行号，供人工复查。
2. **参考坐标转换（liftover）**：把按**旧参考**调用的单样本 VCF，依据制表符坐标映射表
   整批转换到**新参考**，输出标注源行号与 GT 的目标 VCF；反向块自动反向互补并重新锚定。

## 双 VCF 核对

### 功能与判定规则

- **输入限制**
  - 单条染色体参考 FASTA，纯 A/C/G/T 序列，不超过 **1 MiB**；
  - 两份单样本 VCF（恰好 10 列、1 个样本），各不超过 **256 KiB**。
- **接受的变异**
  - A/C/G/T 构成的**单碱基替换**（SNV）；
  - 事件长度不超过 **20 bp** 的**纯插入 / 纯缺失**；
  - GT 必须是**非缺失、未分相的二倍体**（`0/1`、`1/1`，多 ALT 可为 `1/2`）。
  - 其他形态（MNV、符号化 ALT、`*`、分相 `|`、缺失 `.`、非二倍体等）一律判失败。
- **归一化**：对照参考序列对 indel 做左对齐（与 `bcftools norm` / `vt normalize` 等价），
  重复区域中位置写法不同的同一插入 / 缺失归为同一等位变异；SNV 位置不变。
- **剂量计算**：多 ALT 记录按 GT 分别计算每个 ALT 的 0/1/2 剂量，剂量为 0 的变异忽略；
  同一文件内同一归一变异的多条记录剂量累加，超过二倍体上限 2 即判定为**互相矛盾**。
- **整次失败**（返回 422，只显示错误、不展示任何局部结论）：
  - 任一记录 REF 与参考序列不符（含坐标越界、CHROM 不一致）；
  - 同一文件中同一归一变异存在互相矛盾的剂量；
  - 以及任何不满足上述输入规则的情形。
- **结果分组**（均按归一变异排序，并标注两侧原始行号与原始 `POS:REF>ALT` 写法）：
  两侧一致 / 仅 VCF A 或 B 出现 / 剂量不符。

## 目录结构

```
app/            FastAPI 应用与核心逻辑（bio.py / convert.py 仅用标准库）
  main.py       路由与上传校验（/compare 核对、/convert 转换）
  bio.py        FASTA/VCF 解析、左对齐归一化、剂量与分组
  convert.py    映射块解析、正反向坐标转换、反向互补重锚定、目标 VCF 生成
  templates/    上传页与结果页（Jinja2）
tests/          单元 / 端到端测试（bio 41、convert 43、web 16）
samples/        可直接试用的演示数据
scripts/        演示数据生成器（坐标由参考自动计算）
Dockerfile, docker-compose.yml, requirements.txt
```

## 参考坐标转换（旧参考 → 新参考）

### 输入

- 旧、新两份**单条染色体** FASTA（纯 A/C/G/T，各不超过 1 MiB）；
- 制表符映射表（不超过 256 KiB）：每行一个映射块，5 列为
  `旧起点 ⇥ 旧终点 ⇥ 新起点 ⇥ 新终点 ⇥ 方向`，坐标均为 **1 基闭区间**，
  方向取 `+`（正向）或 `-`（反向）；首行允许写列名表头
  `old_start old_end new_start new_end strand`；
- 按旧参考调用的**单样本** VCF（不超过 256 KiB），每行**恰好一个 ALT**。

映射表约束：每个块两侧区间**等长**；同一侧各块之间**互不重叠**（端点相接允许）。

### 转换规则

- 沿用核对工具的文件大小限制、旧参考 **REF 校验**与**未分相二倍体 GT** 规则；
  仅接受单碱基替换，或事件长度不超过 **20 bp** 的纯插入 / 纯缺失。
- 变异的旧参考 REF 区间必须被**同一个映射块完整覆盖**；跨块或落在映射缺口内即失败。
- 正向块坐标一一对应；**反向块**先把坐标镜像到新链，再把 REF/ALT **反向互补**，
  插缺按 VCF 锚定约定在新链上重新锚定，最后统一对照新参考做 REF 校验与左对齐归一化。
- 同一目标归一变异**合并剂量**：目标 VCF 用 `INFO/SRC` 标注全部源 VCF 物理行号，
  GT 写合并后的 `0/1` 或 `1/1`（剂量 0 的 `0/0` 源记录不产生变异）。
- **整批失败**（HTTP 422，只显示错误、不生成任何部分文件）：变异跨块、无映射、
  目标 REF 与新参考不符、合并剂量超过二倍体上限 2，以及任何不满足输入规则的情形。

### 试用

浏览器打开 <http://localhost:8000/convert>，依次选择
`samples/old_ref.fa`、`samples/new_ref.fa`、`samples/liftover.map.tsv`、
`samples/source.vcf`，提交后浏览器直接下载 `target.vcf`。

也可以走 JSON/文件流接口（成功返回 VCF 文件流，失败返回 422 JSON）：

```bash
curl -OJ -F old_fasta=@samples/old_ref.fa \
        -F new_fasta=@samples/new_ref.fa \
        -F mapping=@samples/liftover.map.tsv \
        -F vcf=@samples/source.vcf \
        http://localhost:8000/api/convert
```

演示数据由 `python3 scripts/make_liftover_demo.py` 生成，覆盖正向 / 反向块上的
SNV、插入、缺失，以及重排后按新坐标排序与纯合剂量输出。

## 用 Docker Compose 启动（推荐）

```bash
docker compose up --build
```

浏览器入口：

- 上传页面：<http://localhost:8000/>
- 交互式 API 文档（Swagger UI）：<http://localhost:8000/docs>

停止：`docker compose down`。容器以只读根文件系统、`no-new-privileges` 运行，
内置 `/` 健康检查。

## 本地运行（无 Docker 时）

```bash
python3 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

## 试用

打开 <http://localhost:8000/>，分别选择：

- 参考：`samples/demo.fa`
- VCF A：`samples/caller_a.vcf`
- VCF B：`samples/caller_b.vcf`

演示数据覆盖：SNV 一致、重复区内插入 / 缺失的两种等价写法归一、多 ALT 按 GT 拆分、
`0/0` 剂量 0 被忽略、剂量不符（1/1 vs 0/1）、仅单侧出现。

也可以走 JSON 接口：

```bash
curl -F 'fasta=@samples/demo.fa' \
     -F 'vcf_a=@samples/caller_a.vcf' \
     -F 'vcf_b=@samples/caller_b.vcf' \
     http://localhost:8000/api/compare
```

## 测试

```bash
# 核心逻辑（仅标准库）
python3 -I -m unittest discover -s tests -t . -p 'test_bio.py' -v

# 全部测试（Web 层需要先安装 requirements.txt）
python3 -I -m unittest discover -s tests -t . -v
```
