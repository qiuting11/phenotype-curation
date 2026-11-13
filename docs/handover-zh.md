# 交接文档：表型数据整理项目

**日期**: 2026-11-06
**仓库**: https://github.com/qiuting11/phenotype-curation
**代码位置**: `E:\ncbi\phenotype-curation\`
**数据**: 19个宏基因组学项目，9142个ENA runs，6592个唯一biosamples

---

## 一、项目背景

这个 pipeline 做的事情：给定一批宏基因组学研究的 BioProject 编号，自动生成一份标准化的**表型CSV**（phenotype_clean.csv），其中每个 biosample 都标注了它在研究中的时间阶段和是否处于活跃干预期（`is_intervention`）。

`is_intervention` 是下游差异丰度分析的关键变量——你需要知道一个样本是来自基线/洗脱期（`no`）还是活跃干预期（`yes`）才能做正确的纵向比较。ENA 元数据从不提供这一信息，必须从论文中手动提取。

---

## 二、项目文件结构

```
phenotype-curation/
├── cli.py                      # 入口: phenotype-curation {download,register,classify,build,validate,supplement}
├── core/
│   ├── types.py                # 数据类型定义
│   ├── http_client.py          # API 请求（限速+重试）
│   ├── normalizer.py           # 字段名标准化（14字段同义词表）
│   ├── schemas.py              # JSON Schema 验证
│   └── supplement_parser.py    # 补充材料解析（xlsx/csv/pdf）
├── phases/
│   ├── phase1_fetch.py         # 从 ENA 下载元数据
│   ├── phase2_registry.py      # 自动生成项目注册表
│   ├── phase3_classify.py      # 论文分类（LLM辅助）
│   ├── phase4_build.py         # 构建表型表
│   └── phase5_validate.py      # 7项校验
├── schemas/                    # 4个JSON Schema数据契约
└── templates/
    └── new_project_checklist.md # 新项目接入指南
```

---

## 三、核心数据表说明

### 3.1 表型表（phenotype_clean.csv）

每行 = 一个 biosample。36列。这是最终输出，也是需要校验的核心文件。

**必备字段（always populated）：**

| 列名 | 含义 | 示例 | 来源 |
|------|------|------|------|
| `bioproject` | BioProject 编号 | `PRJNA594393` | ENA |
| `biosample_id` | BioSample 编号 | `SAMN13514848` | ENA |
| `run_accession` | Run 编号（多个用`;`分隔） | `SRR11181927;SRR11181928` | ENA |
| `subject_id` | 受试者ID（解析后的） | `K001`, `S001` | 样本名解析 |
| `sample_id_local` | 原始样本名 | `K001_BL` | ENA sample_alias |

**干预分类字段：**

| 列名 | 含义 | 示例 | 来源 |
|------|------|------|------|
| `intervention_type` | 干预大类 | `fmt`, `diet_prebiotic`, `immunotherapy`, `observational` | Phase 3 registry |
| `intervention_agent` | 具体干预物 | `FOS (fructo-oligosaccharides)` | Phase 3 registry |
| `intervention_group` | 样本所属分组 | `FOS_Cohort_1`, `Placebo`, `Control` | ENA attributes |
| `is_intervention` | **是否处于活跃干预期** | `yes`, `no`, `unknown` | Phase 4 决策树（见第四节） |
| `disease` | 疾病/健康状态 | `healthy`, `NSCLC`, `melanoma`, `T2D` | Phase 3 registry |

**时间轴字段：**

| 列名 | 含义 | 示例 | 来源 |
|------|------|------|------|
| `timepoint_label` | 时间阶段标签 | `baseline`, `feeding_1`, `washout`, `on_treatment` | Phase 4 解析 |
| `timepoint_index` | 时间点数值 | `0`, `14`, `28` | 样本名解析 |
| `timepoint_unit` | 时间单位 | `day`, `week`, `timepoint_index` | 样本名推断 |

**宿主元数据字段：**

| 列名 | 含义 | 示例 | 来源 |
|------|------|------|------|
| `host_age` | 年龄 | `45` | ENA attributes |
| `host_sex` | 性别 | `male`, `female` | ENA attributes |
| `host_bmi` | BMI | `24.5` | ENA attributes |
| `ethnicity` | 民族/种族 | （多为空） | 补充材料 |
| `clinical_outcome` | 临床结局 | （多为空） | 补充材料 |

**样本配对字段：**

| 列名 | 含义 | 示例 |
|------|------|------|
| `paired_pre_sample` | 同一受试者的前一时间点 | `SAMN13514848` |
| `paired_post_sample` | 同一受试者的后一时间点 | `SAMN13514851` |
| `sample_seq` | 受试者内样本序列号 | `1`, `2`, `3` |
| `total_samples` | 受试者总共的样本数 | `5` |
| `pairing_confidence` | 配对可信度 | `confirmed`, `candidate`, `unpaired` |

**审计/溯源字段：**

| 列名 | 含义 | 示例 |
|------|------|------|
| `phenotype_source` | 表型数据来源 | `ena`, `ncbi_biosample`, `supplement_xlsx` |
| `intervention_source` | 干预分类来源 | `registry`, `paper_figure_1`, `ena_metadata` |
| `subject_id_confidence` | 受试者ID解析可信度 | `parsed`, `parsed_cohort_subject`, `candidate_only_not_verified` |
| `subject_id_status` | ID分配状态 | `parsed_candidate` |
| `context_role` | 样本的角色 | `donor`, `mum`, `baby` |
| `context_id` | 关联ID | 供体ID或母亲ID |
| `family_id` | 家庭组ID | 配对研究中使用的 |
| `notes` | 备注 | 协议偏差、映射问题等 |

---

### 3.2 干预注册表（intervention_registry.csv）

每行 = 一个项目的论文验证分类。15列。**这是整个 pipeline 的核心——Phase 4 依赖它来决定 is_intervention。**

| 列名 | 含义 | 必须人工校验？ | 示例 |
|------|------|:---:|------|
| `bioproject` | 项目编号 | 否 | `PRJNA594393` |
| `study_design` | 研究设计类型 | **是** | `RCT`, `RCT_crossover`, `observational_cohort`, `dose_escalation` |
| `intervention_type` | 干预类型 | **是** | `fmt`, `diet`, `immunotherapy`, `observational` |
| `arms` | 研究组（`;`分隔） | **是** | `FMT;Placebo` |
| `control_type` | 对照组类型 | **是** | `placebo`, `habitual_diet`, `run_in_baseline` |
| `randomization` | 是否随机 | 是 | `yes`, `no` |
| `blinding` | 盲法类型 | 是 | `double_blind`, `open_label` |
| `is_one_time` | **是否一次性干预** | **关键！** | `yes` → 所有样本 is_intervention=no |
| `sample_timing` | 采样时间段概括 | 是 | `run_in_2wk_feeding_4wk_washout_4wk` |
| `sample_timing_schedule` | 采样时间点列表 | **关键！** | `day -14, 0, 7, 14, 28, 56` |
| `days_since_intervention` | 干预后天数数组 | 有用 | `[0, 14, 28, 56]` |
| `primary_pmid` | PubMed ID | 否 | `32620143` |
| `primary_pmc` | PMC ID | 否 | `PMC7386241` |
| `classification_confidence` | 分类可信度 | 元数据 | `verified`, `good`, `moderate`, `low` |
| `notes` | 备注 | 自由 | 关键设计细节 |

---

### 3.3 项目注册表（project_registry.csv）

每行 = 一个项目。Phase 2 自动生成，Phase 3 前人工复核。

| 列名 | 含义 | 来源 |
|------|------|------|
| `bioproject_accession` | 项目编号 | ENA |
| `keep_flag` | 是否保留 | 人工 |
| `scope_basis` | 纳入原因 | 人工 |
| `primary_pmid` | PubMed ID | Phase 2 自动/人工 |
| `run_n_local` | 本地下载的 run 数 | Phase 1 |
| `biosample_n_local` | 本地 biosample 数 | Phase 1 |
| `has_intervention_metadata` | ENA 元数据中有干预信息？ | Phase 2 自动检测 |
| `has_timepoint_metadata` | ENA 元数据中有时间点信息？ | Phase 2 自动检测 |
| `raw_group_values` | 检测到的分组值 | Phase 2 自动提取 |
| `study_design_candidate` | 自动推测的研究设计 | Phase 2 自动 |
| `manual_review_status` | 人工审核状态 | 人工 |
| `supplement_status` | 补充材料状态 | 人工 |

---

## 四、is_intervention 决策逻辑

这是整个 pipeline 最核心的部分。每个样本的 `is_intervention` 通过以下 4 条规则确定：

```
Rule 1: 观察性研究？ → 全部 no
Rule 2: 一次性干预？  → 全部 no（样本都在干预后，不是"活跃干预期"）
Rule 3: 有 phase_map？→ 按 day 范围映射（用于FOS/PDX剂量递增等）
Rule 4: 其他          → 按 timepoint_label 推断
```

**is_one_time 的判断标准（非常重要，容易出错）：**

- **一次性干预 = yes**：干预只发生一次，所有样本在之后追踪。如：单次 FMT 移植、HSCT 干细胞移植、阴道播种。
- **持续性干预 = no**：干预持续进行（每日/每周给药、持续饮食）。如：膳食纤维补充、抗PD-1免疫治疗。

**timepoint_label → is_intervention 映射：**

| label | is_intervention | 示例场景 |
|-------|:---:|----------|
| `baseline`, `screening`, `run_in`, `pre` | `no` | 干预前 |
| `feeding_1`, `feeding_2`, `on_treatment`, `during_diet` | `yes` | 活跃干预期 |
| `washout`, `followup`, `post_treatment`, `post` | `no` | 干预后 |
| `longitudinal_unspecified` | `unknown` | 无法判定 |

---

## 五、19个项目清单与现状

| # | BioProject | 设计 | 干预 | 样本数 | is_intervention 状态 |
|---|-----------|------|------|--------|---------------------|
| 1 | PRJEB22863 | RCT | anti-PD-1 | 42 | ✅ T0=no, T1/T2=yes |
| 2 | PRJEB33500 | RCT_crossover | Mediterranean diet | 109 | ⚠️ 需补充数据确认 |
| 3 | PRJEB39610 | observational | observational | 148 | ✅ 全部 no |
| 4 | PRJEB48663 | observational | observational | 51 | ✅ 全部 no |
| 5 | PRJEB75958 | RCT | FMT capsules | 73 | ✅ baseline=no, during=yes, post=no |
| 6 | PRJNA397450 | observational | observational | 121 | ✅ 全部 no |
| 7 | PRJNA541981 | observational | anti-PD-1/CTLA-4 | 47 | ⚠️ 需补充数据确认 |
| 8 | PRJNA594393 | dose_escalation | FOS (K001) | 947 | ✅ phase_map 精确映射 |
| 9 | PRJNA594610 | dose_escalation | PDX (K002) | 947 | ✅ phase_map 精确映射 |
| 10 | PRJNA594613 | dose_escalation | FOS→PDX crossover | 947 | ✅ phase_map 精确映射 |
| 11 | PRJNA594620 | dose_escalation | PDX→FOS crossover | 947 | ✅ phase_map 精确映射 |
| 12 | PRJNA637785 | RCT | FMT (single dose) | 158 | ✅ is_one_time=yes, 全部no |
| 13 | PRJNA701480 | RCT | vaginal seeding | 49 | ✅ is_one_time=yes, 全部no |
| 14 | PRJNA737472 | RCT | FMT / ENMA | 168 | ⚠️ 需补充数据确认 |
| 15 | PRJNA762360 | observational | anti-PD-1 | 63 | ⚠️ 横断面，全部unknown |
| 16 | PRJNA895487 | observational | observational | 28 | ✅ 全部 no |
| 17 | PRJNA904790 | observational | FMT multi-dose | 496 | ⚠️ 需补充数据确认 |
| 18 | PRJNA929328 | RCT | HSCT + FMT | 113 | ⚠️ is_one_time 需确认 |
| 19 | PRJNA929986 | observational | observational | 172 | ✅ 全部 no |

> ⚠️ 标记的项目在 Phase 5 校验中有 warning——`is_intervention=yes` 样本数为 0 且预期应有 yes。需要下载 Strategy 2 样本属性数据（`phenotype-curation download --sample-attrs`）或人工检查论文补充材料。

---

## 六、需要人工校验的内容

### 第一优先级（影响 is_intervention 正确性）

1. **`is_one_time` 分类**：对照每个项目的论文 Methods，确认干预是一次性（is_one_time=yes）还是持续性（is_one_time=no）。判错会导致整批样本的 is_intervention 全错。

2. **`sample_timing_schedule`**：从论文 Table 1 / Figure 1 提取精确的采样时间点（如 day 0, 7, 14, 28, 56）。这对 day-level 精度项目（FOS/PDX）至关重要。

3. **`study_design` 和 `intervention_type`**：确认分类正确（特别是 observational vs RCT 的区分）。

### 第二优先级（影响元数据完整性）

4. **`disease` 标签**：检查 `phase4_build.py` 中 `DISEASE_MAP` 的 19 个项目疾病标签是否准确。

5. **`intervention_agent`**：检查 `INTERVENTION_AGENT_MAP` 中的干预物描述是否准确。

6. **`arms` 和 `control_type`**：确认研究中实际的分组和对照类型。

7. **`context_role` 和 `context_id`**：对于涉及捐赠者/母婴配对的项目，确认角色标注正确。

### 第三优先级（影响下游可用性）

8. **subject_id 解析**：抽查部分项目的 subject_id 是否正确解析（特别是使用非标准命名的项目）。

9. **纵向配对**：检查 `pairing_confidence` 为 `candidate` 的样本是否正确配对。

---

## 七、当前已知不足

### 功能不足

| 问题 | 影响 | 可能解决方向 |
|------|------|-------------|
| `core/types.py` 遮蔽 stdlib 的 `types` 模块 | `python core/*.py` 直接运行失败 | 重命名 `core/types.py` 为 `core/data_types.py` |
| PMID→BioProject 链接 ~50% 覆盖率 | 部分论文无法自动找到数据 | 改进 ENA 搜索或利用 SRA Run Selector |
| Phase 3 仍需人工读论文 | 无法全自动 | 设计上保留的人类判断环节 |
| 无 pytest 单元测试 | 改动后无法自动回归验证 | 添加 test_/ 目录 |
| 没有 CI/CD | 不便于协作 | 配置 GitHub Actions |

### 数据不足

| 问题 | 涉及项目 | 影响 |
|------|---------|------|
| 5个持续干预 RCT 的 `is_intervention=yes` 为 0 | #2,7,14,15,17 | Phase 5 校验 warning |
| 部分项目缺少 per-sample attributes | 多个 | host_age/sex/BMI 和数据来自 ENA filereport 而非 sample XML |
| 横断面研究无法区分时间阶段 | #15 | 全部 `is_intervention=unknown` |
| 补充材料表格未完全提取 | 多个 | ethnicity、clinical_outcome 字段多为空 |

### 设计债务

| 问题 | 说明 |
|------|------|
| Phase maps / overrides 硬编码在代码里 | `PHASE_MAPS`, `TIMEPOINT_OVERRIDES` 等字典嵌在 `phase4_build.py` 中，每次新增项目需要改代码 |
| DISEASE_MAP / INTERVENTION_AGENT_MAP 硬编码 | 理想情况下应从 intervention_registry.csv 读取 |
| subject_id 解析规则是字符串匹配 | `parse_subject_id()` 中 5 种匹配模式是硬编码的正则, 非标准命名需要改代码 |
| 代码和论文结论之间的追溯链不完整 | `traceability_row` schema 已定义但实际 traceability CSV 可能未完整填充 |

---

## 八、给接手同事的操作指南

### 安装

```bash
cd E:\ncbi\phenotype-curation
pip install -e .
```

### 校验现有输出（最重要）

```bash
# 校验表型表
phenotype-curation validate phenotype_clean.csv \
  --registry intervention_registry.csv \
  --out validation_report.json

# 查看报告
cat validation_report.json
```

### 接入新项目

```bash
# 1. 下载
phenotype-curation download --projects PRJNA1234567 --sample-attrs

# 2. 自动注册
phenotype-curation register ena_metadata/ena_metadata_master.csv \
  --sample-attrs-dir ena_metadata/

# 3. 生成分类提示
phenotype-curation classify --registry project_registry.csv

# 4. 【人工】读论文, 填写分类信息到 intervention_registry.csv
#    关键字段: study_design, intervention_type, is_one_time, sample_timing_schedule

# 5. 构建
phenotype-curation build ena_metadata/ena_metadata_master.csv \
  --registry intervention_registry.csv --out phenotype_clean.csv

# 6. 校验
phenotype-curation validate phenotype_clean.csv \
  --registry intervention_registry.csv
```

### 快速检查特定项目

```bash
# 查看某项目的所有样本时间点分布
python -c "
import csv
with open('phenotype_clean.csv') as f:
    for row in csv.DictReader(f):
        if row['bioproject'] == 'PRJNA594393':
            print(f\"{row['subject_id']:6s} {row['timepoint_label']:20s} {row['timepoint_index']:5s} is_intv={row['is_intervention']}\")
"
```

---

## 九、快捷参考

- **JSON Schema**（数据格式定义）：`schemas/` 目录下 4 个 `.json` 文件
- **新项目接入详细指南**：`templates/new_project_checklist.md`
- **API 端点参考**：同上文件末尾
- **所有记忆文件**：`C:\Users\qiuting1\.claude\projects\E--ncbi\memory\`
- **主要入口**：`cli.py`（7个子命令）

### 常用 API 端点

| 用途 | URL |
|------|-----|
| ENA 项目 run 列表 | `https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA123&result=read_run&fields=...&format=tsv` |
| ENA 样本 XML | `https://www.ebi.ac.uk/ena/browser/api/xml/SAMN123456` |
| NCBI BioSample | `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&id=SAMN123&retmode=xml` |
| PMID→PMCID | `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=pubmed&db=pmc&id=12345678` |
| PMC 全文 XML | `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=PMC1234567&rettype=xml` |