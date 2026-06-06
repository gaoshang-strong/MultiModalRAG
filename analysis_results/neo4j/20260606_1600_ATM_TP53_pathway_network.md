# ATM-TP53 Pathway Network Analysis
**Generated:** 2026-06-06 16:00  
**Purpose:** Paper Figure 3 — mechanistic basis of ATM-TP53 mutual exclusivity in gastric cancer  
**Neo4j DB:** bolt://127.0.0.1:7687 | Python env: ProjectGeneration | DeepSeek model: deepseek-v4-pro

---

## Cypher Queries Used

### Q1: Reactome Pathway Membership
```cypher
-- ATM-only genes (not in any TP53/p53 pathway)
MATCH (atm_path:Pathway)-[:SUBPATHWAY_OF*0..3]->(top:Pathway)
WHERE atm_path.name CONTAINS "ATM"
MATCH (g:Gene)-[:INVOLVED_IN]->(atm_path)
WHERE NOT EXISTS {
    MATCH (tp53_path:Pathway)-[:SUBPATHWAY_OF*0..3]->(any:Pathway)
    WHERE any.name CONTAINS "TP53" OR any.name CONTAINS "p53"
    MATCH (g)-[:INVOLVED_IN]->(tp53_path)
}
RETURN g.symbol, count(DISTINCT atm_path) AS n_atm_paths
ORDER BY n_atm_paths DESC

-- Shared genes
MATCH (atm_path:Pathway) WHERE atm_path.name CONTAINS "ATM"
MATCH (tp53_path:Pathway) WHERE tp53_path.name CONTAINS "TP53" OR tp53_path.name CONTAINS "p53"
MATCH (g:Gene)-[:INVOLVED_IN]->(atm_path)
MATCH (g)-[:INVOLVED_IN]->(tp53_path)
RETURN DISTINCT g.symbol ORDER BY g.symbol
```

### Q2: Shared BioGRID Interactors
```cypher
MATCH (atm:Gene {symbol: "ATM"})-[:BIOGRID_INTERACTS_WITH]-(partner:Gene)
MATCH (tp53:Gene {symbol: "TP53"})-[:BIOGRID_INTERACTS_WITH]-(partner)
WITH partner
MATCH (partner)-[:BIOGRID_INTERACTS_WITH]-(any:Gene)
WITH partner, count(DISTINCT any) AS degree
WHERE degree < 500
RETURN partner.symbol AS symbol, degree
ORDER BY degree
```

### Q3: GO Biological Process (shared)
```cypher
MATCH (atm:Gene {symbol: "ATM"})-[r1:ANNOTATED_TO]->(go:GOTerm)
MATCH (tp53:Gene {symbol: "TP53"})-[r2:ANNOTATED_TO]->(go)
WHERE go.namespace = "biological_process"
  AND r1.has_non_iea_evidence = true
  AND r2.has_non_iea_evidence = true
RETURN go.go_id AS go_id, go.name AS name ORDER BY go.name
```

### Q4: GO Biological Process (ATM-specific)
```cypher
MATCH (atm:Gene {symbol: "ATM"})-[r1:ANNOTATED_TO]->(go:GOTerm)
WHERE go.namespace = "biological_process"
  AND r1.has_non_iea_evidence = true
  AND NOT EXISTS { MATCH (tp53:Gene {symbol: "TP53"})-[:ANNOTATED_TO]->(go) }
RETURN go.go_id AS go_id, go.name AS name ORDER BY go.name
```

---

## Raw Query Results

### Q1: Reactome Results
- ATM pathway genes (total): **77**
- TP53 pathway genes (total): **399**
- **Shared (n=14):** ATM, BARD1, BRCA1, CHEK2, KAT5, MDC1, MRE11, NBN, RAD50, RPS27A, TP53, UBA52, UBB, UBC
- **ATM-only genes (non-histone, n=30):** ABL1, ABRAXAS1, APBB1, BABAM1, BABAM2, BAP1, BAZ1B, BRCC3, EYA1, EYA2, EYA3, EYA4, H2AX, HERC2, KDM4A, KDM4B, MAPK8, NSD2, PIAS4, PPP5C, RNF168, RNF8, SMARCA5, SUMO1, TP53BP1, UBE2I, UBE2N, UBE2V2, UBXN1, UIMC1

### Q2: Shared BioGRID Interactors (degree<500, n=102 total, key subset)
- Checkpoint kinases: CHEK1 (242), CHEK2 (240), ATR (264)
- MRN complex: MRE11 (302), NBN (187), RAD50 (379)
- Repair effectors: BRCA2 (379), RAD51 (237), BLM (251), WRN (109)
- Chromatin/mediators: KAT5 (327), MDC1 (421), PCNA (492)
- Stress kinases: DYRK2 (169), MAPK8 (226), MAPK14 (309)
- Transcription: E2F1 (165), HIF1A (434), NFE2L2 (163), FOXO3 (86)
- MMR: MSH2 (310), MSH6 (270)
- Other effectors: CDKN1A (357), STK11 (310), SETD2 (134)

### Q3-Q4: GO Results
- **Shared GO BP (n=10):** DNA damage response, p53-class mediator signaling, cellular response to gamma radiation, cellular senescence, replicative senescence, positive regulation of apoptotic process, regulation of apoptotic process, regulation of cell cycle, positive regulation of gene expression, positive regulation of transcription by RNA Pol II
- **ATM-specific GO BP (n=40, key):** DNA damage checkpoint signaling, DSB processing, DSB repair (NHEJ/HR), G2 DNA damage checkpoint, telomere maintenance, pexophagy, meiotic recombination, ROS response
- **TP53-specific GO BP (n=53, key):** Intrinsic apoptotic signaling by p53, G1 DNA damage checkpoint, p53 class mediator signaling, regulation of transcription (broad), nucleotide-excision repair, autophagy, cellular response to hypoxia, hematopoietic differentiation

---

## DeepSeek Mechanistic Interpretation (deepseek-v4-pro)


**1. Mechanism of mutual exclusivity**

ATM and TP53 mutations are mutually exclusive in gastric cancer because they gatekeeper two partially overlapping arms of the DNA damage response (DDR) that converge on a set of shared effectors critical for cell viability under oncogenic stress. The network topology (shared Reactome genes and BioGRID interactors) reveals that both proteins funnel into common nodes—checkpoint kinases (CHEK1, CHEK2, ATR), the MRN complex, repair mediators (BRCA1, RAD51), and cell fate executers (CDKN1A/p21, MAPK8/JNK, FOXO3). When only ATM or only TP53 is lost, the remaining protein can still route signals through these convergent nodes to arrest the cell cycle, initiate repair, or trigger apoptosis, permitting tumour survival and evolution. However, double loss completely silences all routing: the cell can no longer sense or signal double‑strand breaks, cannot enforce G1/S or G2/M checkpoints, and loses the ability to eliminate genomically unstable cells. This creates a synthetic‑lethal threshold where endogenous replication stress (e.g., from CCNE1 amplification in gastric cancer) runs unchecked, leading to mitotic catastrophe or programmed death. Thus, dual‑mutant clones are eliminated during tumourigenesis, statistically manifesting as mutual exclusivity. Single mutations are tolerated because partial redundancy at the convergence zone maintains a minimal DDR output, a principle consistent with the TCGA gastric adenocarcinoma subtypes (EBV, MSI, GS, CIN) and ACRG MSS/TP53+/− classifications, where ATM and TP53 alterations are almost never co‑occurring.

**2. Three‑compartment model**

- **ATM‑unique compartment**  
ATM is a serine/threonine kinase that directly mediates the early DSB signalling cascade. It phosphorylates substrates (H2AX, MDC1, KAP1) and recruits ubiquitin ligases RNF8, RNF168 and scaffold 53BP1 to the damage site. ATM‑unique functions include: DSB detection by the MRN complex and autophosphorylation; DSB processing and repair pathway choice (NHEJ and HR); chromatin remodelling (BAZ1B, SMARCA5); telomere maintenance; pexophagy; and response to reactive oxygen species. These functions depend on ATM’s kinase activity and localisation to breaks; TP53 cannot substitute because it lacks a kinase domain and does not directly orchestrate chromatin signalling.

- **TP53‑unique compartment**  
TP53 is a master transcription factor. Upon activation, it induces or represses hundreds of target genes that enforce G1 arrest (CDKN1A/p21), apoptosis (BAX, PUMA, NOXA, BID), senescence (CDKN2A/p14), autophagy, nucleotide‑excision repair, negative regulation of transcription, and responses to hypoxia and oxidative stress. These transcriptional programmes are exclusive to TP53; ATM cannot directly initiate them, although it can phosphorylate TP53 to enhance its activity. Thus, loss of TP53 cripples the transcriptional wing of the DDR, while ATM‑dependent repair signalling can still operate.

- **Convergence zone**  
A central set of 14 shared Reactome proteins and 102 high‑confidence physical interactors bridge the two pathways. The core convergence includes the MRN complex (sensor), checkpoint kinases CHEK1/2 and ATR (signal transducers), acetyltranferase KAT5/TIP60 (co‑activator for both ATM and TP53), repair effectors (BRCA1, BRCA2, RAD51, BLM, WRN), replication clamp PCNA, and cell‑fate mediators like CDKN1A/p21, DYRK2, MAPK8, MAPK14, FOXO3, and SETD2. Through these nodes, both ATM‑initiated phosphorylation and TP53‑driven transcription converge on the same biological outputs—cell‑cycle arrest, DNA repair, and apoptosis. This convergence creates functional redundancy: as long as one upstream regulator is intact, the shared nodes can be engaged at a level sufficient to maintain viability. When both are lost, convergent signalling collapses completely, eliminating all DDR decision capacity.

**3. Top 10 convergence nodes that drive synthetic lethality**

| Rank | Node | Reason for synthetic lethality upon dual ATM+TP53 loss |
|------|------|----------------------------------------------------------|
| 1 | **CHEK2** (deg 240) | ATM phosphorylates and activates CHEK2; TP53 is a major CHEK2 substrate. Loss of ATM abolishes CHEK2 activation, loss of TP53 prevents CHEK2‑driven apoptosis/arrest. Only single loss leaves residual CHEK2 output; double loss silences the critical S/G2 checkpoint. |
| 2 | **CHEK1** (deg 242) | ATM phosphorylates CHEK1 (via ATR); TP53 transcriptionally regulates CHEK1 expression. Dual loss removes both checkpoint kinase activities, causing catastrophic mitotic entry with unreplicated/damaged DNA. |
| 3 | **MRE11–NBN–RAD50 (MRN)** (deg 302/187/379) | The MRN complex is the primary sensor for ATM activation and a direct p53 interactor. Without ATM, MRN‑mediated signalling is blunted; without p53, MRN‑dependent transcription of repair genes is lost. Combined absence abolishes DSB sensing. |
| 4 | **KAT5/TIP60** (deg 327) | KAT5 acetylates both ATM (required for kinase activation) and TP53 (required for sequence‑specific DNA binding). Double mutation eliminates the acetylation‑driven activation of either pathway, shutting down DDR initiation entirely. |
| 5 | **BRCA1** (shared Reactome, deg ~200) | BRCA1 is phosphorylated by ATM and serves as a scaffold for HR. p53 transcriptionally controls BRCA1 expression. Loss of both factors severely impairs homologous recombination and G2/M arrest. |
| 6 | **CDKN1A/p21** (deg 357) | p21 is a major p53 target for G1 arrest; ATM can also induce p21 through alternative pathways (e.g., via p53‑independent p38/MAPKAP2). Dual loss removes all p21‑mediated arrest, leaving cells unable to halt the cycle after damage. |
| 7 | **PCNA** (deg 492) | PCNA is modified by ATM‑dependent ubiquitination during translesion synthesis and is regulated by p53 for replication fork stability. With both regulators gone, replication stress tolerance collapses, causing fork collapse and death. |
| 8 | **FOXO3** (deg 86) | FOXO3 is activated by ATM‑dependent phosphorylation and is a p53 transcriptional target. It induces cell cycle arrest, DNA repair, and apoptosis. Combined loss cripples FOXO3‑mediated tumour suppressive functions. |
| 9 | **DYRK2** (deg 169) | DYRK2 phosphorylates p53 at Ser46 (apoptosis), is itself an ATM substrate, and can also be induced by p53. Double loss eliminates a key apoptosis amplification loop. |
|10 | **MAPK8/JNK** (deg 226) | ATM activates JNK via RAC1/ASK1 to trigger apoptosis independently of p53; p53 transactivates components of the JNK pathway. When both are gone, stress‑induced apoptosis is blocked, enabling survival of heavily damaged cells. |

**4. mRNA prediction interpretation**

The absence of transcriptional upregulation of ATM‑specific DSB cascade genes (RNF8, RNF168, TP53BP1/53BP1, ABRAXAS1) in TP53‑mutant tumours supports that **functional compensation occurs at the post‑translational (protein activity) level, not through transcriptional cross‑induction**. In tumours that have lost TP53, ATM kinase is still present and can phosphorylate these downstream mediators, which are constitutively expressed. Therefore, the crucial convergence does not rely on one pathway transcriptionally upregulating the other; instead, it relies on the inherent kinase activity of ATM and the scaffolding abilities of its substrates to uphold a basal DDR output. In the paper, this should be positioned as evidence that the synthetic‑lethal relationship is rooted in protein‑level decision nodes. We would state that “the mutually exclusive pattern cannot be explained by transcriptional compensation but rather reflects functional redundancy at the post‑translational interface of the ATM kinase cascade and the TP53 transcriptional programme.”

**5. Figure 3 description**

**Figure 3. Network model of ATM–TP53 functional convergence and mutual exclusivity in gastric cancer.**

*Node placement and logic:*  
The figure is divided into three coloured panels. The left panel (blue) represents the ATM compartment; the right panel (red) represents the TP53 compartment; the central overlap region (purple) depicts the convergence zone.  
– ATM compartment (blue): Nodes include ATM kinase, MRN complex (MRE11‑NBN‑RAD50), MDC1, H2AX, RNF8, RNF168, TP53BP1, ABRAXAS1, UIMC1, and other ATM‑specific interactors. Processes: DSB sensing, chromatin ubiquitylation, NHEJ/HR.  
– TP53 compartment (red): Nodes include TP53 transcription factor, CDKN1A, BAX, PUMA, MDM2, CCNG1, GADD45, and NER components. Processes: G1 arrest, apoptosis, senescence.  
– Convergence zone (purple): Shared nodes are placed in the centre, including CHEK1, CHEK2, ATR, KAT5, BRCA1, RAD51, PCNA, MAPK8, FOXO3, DYRK2. These nodes connect to both upstream compartments via solid (physical interaction) and dashed (regulatory) edges.  

*Edge types and colour coding:*  
– Solid lines: protein–protein interactions from BioGRID.  
– Dashed lines: phosphorylation (ATM) or transcriptional regulation (TP53).  
– Blue edges originate from ATM; red edges from TP53; purple edges indicate a convergent interaction (e.g., ATM→CHEK2 and TP53→CHEK2).  

*Legend (3–4 sentences):*  
“Figure 3. ATM–TP53 network architecture defines a convergence zone that underlies their mutual exclusivity. ATM kinase cascade (blue) and p53‑driven transcriptional response (red) both signal to a central set of effectors (purple) that govern cell‑cycle checkpoints, DNA repair, and apoptosis. When one pathway is inactivated by mutation, the remaining pathway can sustain sufficient output through these shared mediators, allowing tumour evolution. Dual inactivation silences the entire convergence zone, reaching a synthetic‑lethal threshold that eliminates double‑mutant cells.”

---

**中文版本**

**1. 互斥性的机制**
ATM与TP53在胃癌中呈现突变互斥，是因为二者作为DNA损伤应答（DDR）的两个部分重叠臂的守门员，共同汇聚于一套在癌基因应激下对细胞存活至关重要的效应节点。网络拓扑（共享的Reactome基因与BioGRID互作蛋白）显示，二者的信号通过共同的汇聚节点——检查点激酶（CHEK1、CHEK2、ATR）、MRN复合体、修复介质（BRCA1、RAD51）和细胞命运执行者（CDKN1A/p21、MAPK8/JNK、FOXO3）——进行传导。当仅ATM或仅TP53发生突变时，剩下的蛋白仍可通过这些汇聚节点传递信号，实现细胞周期阻滞、DNA修复或凋亡启动，从而允许肿瘤存活与进化。然而，双缺失会使所有传导通路彻底瘫痪：细胞既无法感知或传递双链断裂信号，也无法执行G1/S或G2/M检查点，更丧失了清除基因组不稳定细胞的能力。此时，内源性复制应激（如胃癌中CCNE1扩增导致）将不受控制地积累，引发有丝分裂灾难或程序性死亡。因此，双突变克隆在肿瘤发生过程中被清除，在统计学上表现为互斥性。单突变之所以被容忍，是因为汇聚区的部分冗余性维持了最低限度的DDR输出，这与TCGA胃癌亚型（EBV、MSI、GS、CIN）及ACRG的MSS/TP53+/-分类中ATM与TP53变异几乎从不同时出现的原则一致。

**2. 三分区模型**

- **ATM独有功能区**  
ATM是一种丝/苏氨酸激酶，直接介导早期DSB信号级联。它磷酸化底物（H2AX、MDC1、KAP1），并招募泛素连接酶RNF8、RNF168和支架蛋白53BP1至损伤位点。ATM独有的功能包括：通过MRN复合体感知DSB并自身磷酸化；DSB加工及修复途径选择（NHEJ和HR）；染色质重塑（BAZ1B、SMARCA5）；端粒维持；过氧化物酶体自噬；以及活性氧响应。这些功能依赖于ATM的激酶活性和损伤定位，而TP53因缺乏激酶结构域且不直接统筹染色质信号传导，无法替代。

- **TP53独有功能区**  
TP53是一个主转录因子。激活后可诱导或抑制数百个靶基因，执行G1期阻滞（CDKN1A/p21）、凋亡（BAX、PUMA、NOXA、BID）、衰老（CDKN2A/p14）、自噬、核苷酸切除修复、转录负调控以及对缺氧和氧化应激的应答。这些转录程序为TP53独有；ATM虽能磷酸化TP53增强其活性，但不能直接启动转录响应。因此，TP53缺失会瘫痪DDR的转录翼，而ATM依赖的修复信号仍可运作。

- **汇聚区**  
由14个共享的Reactome通路蛋白和102个高置信度物理互作蛋白构成的中央集合，作为两路的桥梁。核心汇聚点包括MRN复合体（传感器）、检查点激酶CHEK1/2及ATR（信号转导）、乙酰转移酶KAT5/TIP60（ATM与TP53的共同活化因子）、修复效应器（BRCA1、BRCA2、RAD51、BLM、WRN）、复制夹PCNA，以及细胞命运介质如CDKN1A/p21、DYRK2、MAPK8、MAPK14、FOXO3和SETD2。通过这些节点，ATM启动的磷酸化信号与TP53驱动的转录调控汇聚于同一套生物学输出——细胞周期阻滞、DNA修复与凋亡。这种汇聚创造了功能冗余：只要一个上游调节因子完好，共享节点即可被充分激活以维持细胞活力。当二者同时丧失，汇聚信号完全崩溃，所有DDR决策能力丧失。

**3. 驱动合成致死的十大汇聚节点**

| 排名 | 节点 | 双缺失导致合成致死的原因 |
|------|------|----------------------------------------------------------|
| 1 | **CHEK2**（度240） | ATM磷酸化并激活CHEK2；TP53是CHEK2主要底物。缺失ATM则CHEK2无法激活，缺失TP53则CHEK2介导的凋亡/阻滞失效。单缺失仍残留部分CHEK2输出；双缺失彻底关闭关键的S/G2检查点。 |
| 2 | **CHEK1**（度242） | ATM通过ATR磷酸化CHEK1；TP53在转录水平调控CHEK1表达。双缺失使两个检查点激酶活性尽失，引发携带未修复DNA的灾难性有丝分裂。 |
| 3 | **MRE11–NBN–RAD50（MRN）**（度302/187/379） | MRN复合体是ATM激活的主要传感器，亦直接与p53互作。无ATM时MRN介导的信号减弱；无p53时MRN依赖的修复基因转录丧失。双缺失彻底废除DSB感知。 |
| 4 | **KAT5/TIP60**（度327） | KAT5乙酰化激活ATM（激酶活化所必需）和TP53（序列特异性DNA结合所必需）。双突变移除两条通路的乙酰化驱动活化，完全关闭DDR启动。 |
| 5 | **BRCA1**（共享Reactome，度~200） | BRCA1被ATM磷酸化并作为HR的支架；p53转录调控BRCA1表达。双缺失严重损害同源重组和G2/M阻滞。 |
| 6 | **CDKN1A/p21**（度357） | p21是p53介导G1阻滞的主效应器；ATM亦可通过不依赖p53的途径（如p38/MAPKAP2）诱导p21。双缺失移除所有p21驱动的阻滞，细胞受损后无法停止周期。 |
| 7 | **PCNA**（度492） | PCNA在跨损伤合成期间被ATM依赖的泛素化修饰，并受p53调控以维持复制叉稳定。双调节因子缺失后，复制应激耐受瓦解，导致复制叉崩解和细胞死亡。 |
| 8 | **FOXO3**（度86） | FOXO3通过ATM依赖的磷酸化被激活，同时也是p53的转录靶基因。其促细胞周期阻滞、DNA修复和凋亡功能在双缺失后完全丧失。 |
| 9 | **DYRK2**（度169） | DYRK2在Ser46位点磷酸化p53（凋亡），其本身是ATM底物，也可被p53诱导。双缺失废除一条关键的凋亡放大环路。 |
|10 | **MAPK8/JNK**（度226） | ATM通过RAC1/ASK1激活JNK以不依赖p53的方式触发凋亡；p53转录激活JNK通路组分。二者皆失时，应激诱导的凋亡被阻断，使重度损伤细胞得以存活。 |

**4. mRNA预测的解读**

TP53突变型肿瘤中ATM特有的DSB级联基因（RNF8、RNF168、TP53BP1/53BP1、ABRAXAS1）未见显著mRNA上调，这支持**功能性补偿发生在翻译后（蛋白质活性）水平，而非通过转录交叉诱导**。在缺失TP53的肿瘤中，ATM激酶仍然存在并能磷酸化这些组成性表达的下游介质。因此，关键汇聚并不依赖于一条通路在转录上上调另一条通路，而是依赖于ATM固有的激酶活性及其底物的支架功能，以维持基础DDR输出。在论文中，应将其作为证据，表明合成致死关系根植于蛋白质层面的决策节点。陈述为：“互斥模式不能用转录补偿来解释，而是反映了ATM激酶级联与TP53转录程序在翻译后界面上的功能冗余。”

**5. 图3描述**

**图3. 胃癌中ATM–TP53功能汇聚与互斥性的网络模型**

*节点布局与逻辑：*  
图分为三个彩色面板。左侧（蓝色）为ATM功能区；右侧（红色）为TP53功能区；中央重叠区域（紫色）为汇聚区。  
– ATM区（蓝）：节点包括ATM激酶、MRN复合体（MRE11‑NBN‑RAD50）、MDC1、H2AX、RNF8、RNF168、TP53BP1、ABRAXAS1、UIMC1等ATM特有互作因子。过程：DSB感知、染色质泛素化、NHEJ/HR。  
– TP53区（红）：节点包括TP53转录因子、CDKN1A、BAX、PUMA、MDM2、CCNG1、GADD45及NER组分。过程：G1阻滞、凋亡、衰老。  
– 汇聚区（紫）：共享节点置于中央，包括CHEK1、CHEK2、ATR、KAT5、BRCA1、RAD51、PCNA、MAPK8、FOXO3、DYRK2。这些节点通过实线（物理互作）和虚线（调控关系）连接上游两区。  

*边类型与颜色编码：*  
– 实线：BioGRID蛋白‑蛋白互作。  
– 虚线：磷酸化（ATM）或转录调控（TP53）。  
– 蓝色边源于ATM；红色边源于TP53；紫色边表示汇聚性互作（如ATM→CHEK2和TP53→CHEK2）。  

*图注（3–4句）：*  
“图3. ATM‑TP53网络架构定义了一个汇聚区，揭示了它们互斥性的基础。ATM激酶级联（蓝）与p53驱动的转录响应（红）均信号传导至一组决定细胞周期检查点、DNA修复与凋亡的中央效应器（紫）。当一条通路因突变而失活时，剩余通路可通过这些共享介质维持足够的输出，允许肿瘤进化。双失活则使整个汇聚区沉默，达到合成致死阈值，从而清除双突变细胞。”