# Neo4j Query: 对于DDR来说，TP53和ATM为什么在不同的通路里，他们各自的通路都有哪些基因？哪些是shared？这两条通路有何不同？为什么人类要进化出两条通路？

**Genes:** ATM, TP53  
**Model:** `deepseek-v4-pro`  
**Generated:** 2026-06-05 17:45:46

---

## Cypher Queries (for reproducibility)

```cypher
Q1: DDR pathways (T1+T2) for ATM and TP53
```

```cypher
Q2: ATM and TP53 full gene profiles (NCBI summary + UniProt full text + all GO BP/MF)
```

```cypher
Q3: All core pathway members (genes in >=2 DDR pathways, full summaries, no limit)
```

```cypher
Q4: ATM-only (7 genes) and TP53-only (22 genes) full annotation
```

---

## DeepSeek Analysis

# ATM and TP53 in the DNA Damage Response: Divergent Pathways, Overlapping Networks

## 1. Why Are ATM and TP53 in Different Pathways?

The knowledge graph data reveals that ATM and TP53 reside in **distinct functional layers** of the DNA damage response (DDR) despite their tight regulatory connection. According to the Reactome pathway annotations (Q1):

- **ATM** is predominantly mapped to **DNA‑damage sensing, repair‑focused pathways** and **checkpoint signalling**. For instance, T1‑leaf pathways for ATM include “Processing of DNA double‑strand break ends”, “Nonhomologous End‑Joining (NHEJ)”, “Recruitment and ATM‑mediated phosphorylation of repair and signaling proteins at DNA double strand breaks”, and “Sensing of DNA Double Strand Breaks” (all T1). These positions ATM directly at the site of damage as an **upstream kinase**.

- **TP53** is enriched in **transcriptional regulation and effector pathways** downstream of damage signals. Its T1‑leaf pathways include “TP53 Regulates Transcription of DNA Repair Genes”, “TP53 Regulates Transcription of Caspase Activators and Caspases”, “TP53 Regulates Transcription of Genes Involved in G1 Cell Cycle Arrest”, and multiple “Regulation of TP53 Activity…” nodes (T1/T2). This confirms TP53’s role as a **transcription factor** that orchestrates gene expression programs following activation (often by ATM).

The mechanistic separation is further supported by GO annotations from Q2:
- ATM has GO terms “DNA damage checkpoint signaling”, “double‑strand break repair via nonhomologous end joining”, “signal transduction in response to DNA damage”, and “histone H2AXS139 kinase activity” (T1/T2-level evidence).
- TP53 has GO terms “DNA damage response, signal transduction by p53 class mediator”, “positive regulation of transcription by RNA polymerase II”, and “intrinsic apoptotic signaling pathway in response to DNA damage” (T1/T2).

Thus, ATM is classified as a **sensor‑transducer kinase**, whereas TP53 is classified as a **transcription‑factor effector**, which naturally places them in separate pathway categories within the Reactome hierarchy.

---

## 2. Pathway‑Associated Genes: ATM vs. TP53

The data define “core pathway members” as genes participating in ≥2 DDR‑related pathways (Q3). These gene sets strongly overlap, but distinct subsets exist.

### 2.1 ATM‑specific core pathway members (ATM‑only)
Genes present only in the ATM core set (Q4, `"ATM_only"`), excluding shared genes, are:
- **MSH3** – mismatch repair (MutSβ)
- **MSH6** – mismatch repair (MutSα)
- **MUTYH** – base excision repair (oxidative damage)
- **NEIL1** – base excision repair (oxidized pyrimidines)
- **NEIL3** – base excision repair (interstrand cross‑links)
- **NTHL1** – base excision repair (oxidized pyrimidines)
- **OGG1** – base excision repair (8‑oxoguanine)

These seven genes are all **DNA glycosylases or mismatch‑recognition proteins** that function in BER and MMR—pathways not directly regulated by TP53 transcription according to the provided data.

### 2.2 TP53‑specific core pathway members (TP53‑only)
Genes present only in the TP53 core set (Q4, `"TP53_only"`) are:
- **FKBPL** – p21 stability regulator
- **HSP90AA1**, **HSP90AB1** – molecular chaperones
- **Tubulin genes**: *TUBA1A, TUBA1B, TUBA1C, TUBA3C, TUBA3D, TUBA3E, TUBA4A, TUBA4B, TUBA8, TUBAL3, TUBB1, TUBB2A, TUBB2B, TUBB3, TUBB4A, TUBB4B, TUBB6, TUBB8, TUBB8B* (19 genes)

These are mainly **chaperones** and **cytoskeletal/mitotic spindle components**. While some of these (e.g. HSP90, tubulins) are not classically considered DDR genes, their inclusion in the TP53‑only set indicates that they participate in TP53‑centric regulatory or effector modules, such as p53 protein stabilization (HSP90) and mitotic arrest execution (tubulins).

### 2.3 Shared core pathway members
The vast majority of core genes (684 out of ~700/group) are shared (Q4 summary). Representative shared genes from Q3 lists include:

| Category | Genes |
|----------|-------|
| **Sensor kinases** | ATM, ATR |
| **Mediators/adaptors** | MDC1, TP53BP1, TOPBP1, RHNO1 |
| **Downstream kinases** | CHEK1, CHEK2, PLK3 |
| **DSB repair** | MRE11, RAD50, NBN, BRCA1, BARD1, RBBP8, RAD51, RAD51B/C/D, XRCC2, XRCC3, BRCA2, PALB2, BRIP1, BLM, WRN, DNA2, EXO1, RMI1/2, TOP3A |
| **Checkpoint clamps** | RAD1, RAD9A/B, HUS1, RAD17, RFC2/3/4/5 |
| **Replication factors** | RPA1/2/3, PCNA |
| **Ubiquitin/chromatin** | RNF8, RNF168, HERC2, UBC, UBB, UBA52, RPS27A, PIAS4, H2AX, many histone genes |
| **TP53 regulators** | MDM2, MDM4, CDKN2A/ARF |
| **CDK‑cyclin inhibitors** | CDKN1A (p21), CDKN1B (p27) |
| **Apoptosis/cell‑death effectors** | BAX, BID, BBC3 (PUMA), PMAIP1 (NOXA), CASP2/6/10, APAF1, PIDD1, CRADD, BIRC5, TP53BP2, PPP1R13B |
| **Cell‑cycle machinery** | CDK1, CDK2, CCNA1, CCNA2, CCNB1, CCNE1, CCNE2, CDC25A/C, SFN, GTSE1 |
| **Transcriptional co‑regulators** | EP300, CREBBP, PRMT5, HIPK1/2, HDAC1/2, and many TAF/TFIID subunits |

This extensive sharing reflects the deep integration of ATM‑initiated signalling and TP53‑executed transcriptional responses.

---

## 3. Key Differences Between the ATM and TP53 Pathways

| Aspect | ATM‑centric pathways | TP53‑centric pathways |
|--------|----------------------|------------------------|
| **Primary function** | DNA damage sensor & signal amplifier | Transcriptional effector deciding cell fate |
| **Molecular activity** | Ser/Thr kinase; phosphorylates target proteins at S/T‑Q motifs (UniProt) | Sequence‑specific transcription factor (GO: “DNA‑binding transcription factor activity”) |
| **Immediate substrates** | H2AX, NBN, BRCA1, CHEK2, TP53, MDC1, RAD50 (NCBI summary) | Binds DNA at consensus sites; activates genes such as CDKN1A, BAX, GADD45A, RRM2B, etc. |
| **Pathway focus (Reactome)** | DSB processing, NHEJ, HR, checkpoint signalling | Cell‑cycle arrest, apoptosis, senescence, metabolic regulation |
| **Regulatory inputs** | Directly activated by MRN complex at DNA ends | Post‑translational modifications (phosphorylation by ATM/ATR, acetylation, methylation, ubiquitination) |
| **Unique core genes** | BER/MMR enzymes (OGG1, MSH6, etc.) – damage repair specialists | Chaperones (HSP90), tubulins – needed for p53 stability and mitotic execution |
| **Temporal role** | Immediate to early (seconds–minutes) | Delayed to sustained (minutes–hours) |

---

## 4. Evolutionary Rationale for Maintaining Two Distinct Pathway Modules

The separation of ATM‑ and TP53‑dominant pathways is not a redundancy but a **functional division of labour** that provides several evolutionary advantages:

1. **Speed vs. Durability**  
   ATM can phosphorylate pre‑existing proteins within seconds to initiate rapid repair and immediate cell‑cycle arrest (evidence: ATM GO term “mitotic G2 DNA damage checkpoint signaling” T1). TP53, in contrast, requires transcription and translation, producing durable effector proteins (e.g., p21) that sustain arrest or commit the cell to apoptosis (GO: “positive regulation of apoptotic process” and “replicative senescence” T1).

2. **Signal Integration and Noise Filtering**  
   ATM acts as a **high‑gain amplifier** of a physically detected lesion (double‑strand break). TP53 integrates multiple upstream signals (ATM, ATR, MAPKs, stress kinases) and translates them into a graded transcriptional output. This prevents a single erroneous break from triggering apoptosis while ensuring that severe or persistent damage leads to cell elimination.

3. **Cancer‑Suppressive Robustness**  
   The two‑tier architecture creates a bottleneck: inactivation of TP53 alone does not abolish all DDR (ATM can still promote repair), and loss of ATM can sometimes be compensated by ATR‑dependent TP53 activation (UniProt function of ATR: “activates checkpoint signalling upon DNA stress”). This redundancy reduces the likelihood that a single mutation disables the entire DDR network.

4. **Tissue‑Specific Adaptation**  
   The TP53‑only set (Q4) contains multiple tubulin isoforms and HSP90 chaperones, suggesting that TP53 co‑regulates genes involved in mitotic spindle integrity and protein homeostasis in a cell‑type‑specific manner. ATM, in contrast, is ubiquitously required for DNA repair in all dividing cells.

5. **Linkage to Other Cellular Processes**  
   TP53’s transcriptional targets extend beyond classical DDR—for example, metabolic genes (Q1: “TP53 Regulates Metabolic Genes”) and ferroptosis‑related pathways. This broad scope allows the DDR to couple with metabolism and redox balance, a versatility that a purely kinase‑based response (ATM) cannot achieve alone.

In summary, the two‑module design reflects an evolutionary strategy that combines **fast, local repair decisions** (ATM) with **slow, global cell‑fate decisions** (TP53). This layered control is a hallmark of the human DDR and explains why both proteins are among the most frequently mutated genes in gastric and other cancers.

---

# 中文版 (Chinese Version)

# ATM与TP53在DNA损伤反应中的不同通路、共享基因及进化意义

## 1. 为什么ATM和TP53属于不同通路？

知识图谱数据显示，ATM和TP53尽管在调控上紧密相连，但在功能层级上属于**不同的DDR模块**（基于Q1的Reactome通路注释）：

- **ATM** 主要出现在**DNA损伤感知与修复通路**以及**检查点信号通路**中。例如，ATM的T1叶子节点通路包括“DNA双链断裂末端的加工”、“非同源末端连接”、“ATM介导的修复蛋白磷酸化招募”和“双链断裂感知”（均为T1级）。这表明ATM作为**上游激酶**直接作用于损伤位点。

- **TP53** 则富集在**转录调控和效应通路**中。其T1叶子节点通路包括“TP53调控DNA修复基因的转录”、“TP53调控Caspase激活因子和Caspases的转录”、“TP53调控G1细胞周期阻滞基因的转录”，以及多个“TP53活性调控”节点（T1/T2）。这证实了TP53作为**转录因子**在接收活化信号（常由ATM传递）后协调下游基因表达程序的角色。

这一功能分离也得到Q2的GO注释的支持：
- ATM的GO术语包括“DNA损伤检查点信号传导”、“非同源末端连接的双链断裂修复”、“DNA损伤响应中的信号传导”和“组蛋白H2AXSer139激酶活性”（T1/T2级证据）。
- TP53的GO术语包括“p53类介导的DNA损伤响应信号传导”、“RNA聚合酶II驱动的转录正调控”和“DNA损伤诱导的固有凋亡信号通路”（T1/T2级）。

因此，ATM被归类为**传感器‑转导激酶**，而TP53被归类为**转录因子效应器**，这使得它们在Reactome层级中处于不同的通路类别。

## 2. 各自通路中的基因

数据将核心通路成员定义为参与≥2个DDR通路中的基因（Q3）。这些基因集高度重叠，但也存在独特子集。

### 2.1 ATM特有核心通路成员（ATM‑only）
仅存在于ATM核心集中的基因（Q4 `"ATM_only"`），排除共享者，包括：
- **MSH3** – 错配修复（MutSβ）
- **MSH6** – 错配修复（MutSα）
- **MUTYH** – 碱基切除修复（氧化损伤）
- **NEIL1** – 碱基切除修复（氧化嘧啶）
- **NEIL3** – 碱基切除修复（链间交联）
- **NTHL1** – 碱基切除修复（氧化嘧啶）
- **OGG1** – 碱基切除修复（8‑氧鸟嘌呤）

这七个基因均为**DNA糖基化酶或错配识别蛋白**，功能集中于BER和MMR，据给定数据，这些通路并非TP53直接转录调控的对象。

### 2.2 TP53特有核心通路成员（TP53‑only）
仅存在于TP53核心集中的基因（Q4 `"TP53_only"`）包括：
- **FKBPL** – p21蛋白稳定性调节因子
- **HSP90AA1**, **HSP90AB1** – 分子伴侣
- **微管蛋白基因**：*TUBA1A, TUBA1B, TUBA1C, TUBA3C, TUBA3D, TUBA3E, TUBA4A, TUBA4B, TUBA8, TUBAL3, TUBB1, TUBB2A, TUBB2B, TUBB3, TUBB4A, TUBB4B, TUBB6, TUBB8, TUBB8B*（19个基因）

这些主要是**分子伴侣**和**细胞骨架/有丝分裂纺锤体成分**。HSP90和微管蛋白并非经典DDR基因，但它们出现在TP53特有集中，表明它们参与了p53稳定性调控（如HSP90）或有丝分裂阻滞执行（如微管蛋白）等TP53特异性效应模块。

### 2.3 共享核心通路成员
绝大多数核心基因（约684个，占总数的绝大部分）是两方共享的（Q4摘要）。根据Q3列表，代表性共享基因包括：

| 类别 | 基因 |
|------|------|
| **传感器激酶** | ATM, ATR |
| **中介/接头蛋白** | MDC1, TP53BP1, TOPBP1, RHNO1 |
| **下游激酶** | CHEK1, CHEK2, PLK3 |
| **DSB修复** | MRE11, RAD50, NBN, BRCA1, BARD1, RBBP8, RAD51, RAD51B/C/D, XRCC2, XRCC3, BRCA2, PALB2, BRIP1, BLM, WRN, DNA2, EXO1, RMI1/2, TOP3A |
| **检查点夹钳** | RAD1, RAD9A/B, HUS1, RAD17, RFC2/3/4/5 |
| **复制因子** | RPA1/2/3, PCNA |
| **泛素化/染色质** | RNF8, RNF168, HERC2, UBC, UBB, UBA52, RPS27A, PIAS4, H2AX, 多个组蛋白基因 |
| **TP53调控因子** | MDM2, MDM4, CDKN2A/ARF |
| **CDK‑cyclin抑制因子** | CDKN1A (p21), CDKN1B (p27) |
| **凋亡/细胞死亡效应器** | BAX, BID, BBC3 (PUMA), PMAIP1 (NOXA), CASP2/6/10, APAF1, PIDD1, CRADD, BIRC5, TP53BP2, PPP1R13B |
| **细胞周期机器** | CDK1, CDK2, CCNA1, CCNA2, CCNB1, CCNE1, CCNE2, CDC25A/C, SFN, GTSE1 |
| **转录共调控因子** | EP300, CREBBP, PRMT5, HIPK1/2, HDAC1/2，以及多个TAF/TFIID亚基 |

这种广泛共享体现了ATM启动的信号传导与TP53执行的转录响应之间的深度整合。

## 3. ATM通路与TP53通路的关键差异

| 方面 | ATM核心通路 | TP53核心通路 |
|------|--------------|----------------|
| **主要功能** | DNA损伤传感器与信号放大器 | 决定细胞命运的转录效应器 |
| **分子活性** | Ser/Thr激酶；在S/T‑Q基序磷酸化靶蛋白（UniProt） | 序列特异性转录因子（GO：“DNA‑binding transcription factor activity”） |
| **直接底物/靶标** | H2AX, NBN, BRCA1, CHEK2, TP53, MDC1, RAD50（NCBI总结） | 结合DNA共有序列；激活CDKN1A, BAX, GADD45A, RRM2B等基因 |
| **通路焦点（Reactome）** | DSB加工，NHEJ，HR，检查点信号传导 | 细胞周期阻滞，凋亡，衰老，代谢调控 |
| **调控输入** | 由DNA末端MRN复合体直接激活 | 翻译后修饰（ATM/ATR磷酸化，乙酰化，甲基化，泛素化） |
| **特有核心基因** | BER/MMR酶（OGG1, MSH6等）– 损伤修复专一性 | 分子伴侣（HSP90），微管蛋白– 参与p53稳定性和有丝分裂执行 |
| **时序角色** | 即刻至早期（秒‑分钟） | 延迟至持续（分钟‑小时） |

## 4. 进化出两条不同通路的原因

ATM‑和TP53‑主导的通路分离并非冗余，而是一种**功能分工**，提供了多个进化优势：

1. **速度与持久性**  
   ATM可以在数秒内磷酸化已有蛋白质，启动快速修复和即时的细胞周期阻滞（证据：ATM的GO术语“有丝分裂G2 DNA损伤检查点信号传导” T1）。TP53则需要转录和翻译，产生持久的效应蛋白（如p21），维持阻滞或使细胞进入凋亡（GO：“凋亡过程的正调控”和“复制性衰老” T1）。

2. **信号整合与噪音过滤**  
   ATM作为物理损伤（双链断裂）的**高增益放大器**，TP53整合多个上游信号（ATM, ATR, MAPKs，应激激酶），将其转化为分级的转录输出。这可以防止单个错误断裂引发凋亡，同时确保严重或持续损伤导向细胞清除。

3. **抑癌鲁棒性**  
   两层架构创建了一个瓶颈：单独灭活TP53不会废除所有DDR（ATM仍可促进修复），而ATM缺失有时可由ATR依赖的TP53激活部分补偿（ATR的UniProt功能描述：“在DNA应激时激活检查点信号”）。这种冗余降低了单基因突变致使整个DDR网络失效的可能性。

4. **组织特异性适应**  
   TP53特有基因集（Q4）包含众多微管蛋白亚型和HSP90伴侣分子，提示TP53可能协同调控参与有丝分裂纺锤体完整性和蛋白质稳态的基因，且具有细胞类型特异性。而ATM在所有分裂细胞中普遍需要。

5. **与其他细胞过程的联系**  
   TP53的转录靶点超出经典DDR，例如代谢基因（Q1：“TP53 Regulates Metabolic Genes”）和铁死亡相关通路。这种广泛的范围使得DDR能够与代谢和氧化还原平衡相耦合，这是纯粹基于激酶的ATM响应无法单独实现的多样性和灵活性。

总之，双模块设计反映了一种进化策略：将**快速的局部修复决策**（ATM）与**缓慢的全局细胞命运决策**（TP53）相结合。这种分层控制是人类DDR的标志，也解释了为何这两个蛋白在胃癌及其他恶性肿瘤中均为最常突变的基因。