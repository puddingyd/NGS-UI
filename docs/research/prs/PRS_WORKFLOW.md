# PRS 流程說明：從基因資料到疾病風險

更新：2026-11-14｜文件版本 v10.36（原 v10.32 的 review／response 見[文件索引](README.md)）

這份文件說明平台如何利用全基因體定序資料，逐步得到「多基因風險分數」、「百分位」與「未來發病風險」，以及開發時為什麼選擇目前的方法。

各章先說概念；供生資審查的實作方法、參數與程式位置，放在可展開的「生資方法與程式對照」中。文內程式路徑均以NGS-UI repo根目錄為起點，完整本機／伺服器路徑及入口索引列在文末。本次補充描述現有實作，沒有更動計分或風險算法。

PRS（Polygenic Risk Score，多基因風險分數）的概念，是把許多常見基因差異對疾病的影響合併考慮。平台採用已發表的研究模型計分，再結合參考族群與台灣疾病統計，提供六種疾病的風險估計。**目前已確認程式能重現指定的計分與研究結果，但還沒有用台灣族群的實際追蹤資料，驗證預測是否準確，結果僅供研究。**

2026-11-14 的後續查核與工具／呈現調整，另見[執行紀錄](PRS_VALIDATION_2026-11-14.md)。本版將已有結果的來源查核納入第七章，並在第九章說明院內健康人可以做的分布比較。全部癌症的分數方向與 HR 對應仍未結案，平台及報告會明示這項限制；本次沒有變更權重、排名或風險參數。

## 一、先分清楚三種結果

從基因資料到發病風險，中間需要回答三個不同的問題：

| 結果 | 代表的意思 | 使用時怎麼理解 |
|---|---|---|
| 原始 PRS score（分數） | 依疾病模型，把個案各個基因位置的影響加總 | 它是計算結果，需要放進參考族群中比較，才能知道高低；不同疾病的分數不能直接比較 |
| PRS percentile（百分位） | 這個分數在參考族群中的排名 | 第80百分位表示分數高於或等於約80%的參考者；不代表有80%的發病機率 |
| Absolute risk（絕對風險） | 從目前年齡起，在指定期間內首次罹患該疾病的估計機率 | 要連同年齡與期間閱讀，例如「未來10年」或「目前至80歲」 |

因此，得到百分位後還不能直接當成發病風險。同樣位於第80百分位的人，會因評估年齡、性別及該疾病在人群中的常見程度不同，而有不同的風險估計。接下來的每一步，就是逐漸補上這些資訊。

假設模型估計「未來10年風險5%」，意思是每100位具有相同評估條件的人，約有5位會在這段期間首次發病。這是解釋用的假設數字；對個案而言，它表達發病的可能性，不是預告一定會不會發病。

## 二、整個流程分成兩段

```text
三級分析
原始 gVCF → 讀取可靠基因資料 → PRS 分數 → 百分位 → 完成祖源與品質檢查

PRS 頁面
百分位＋年齡、性別與各疾病病史
  → 結合 INTERVENE 相對風險及台灣疾病統計
  → 未來10年、未來20年、目前至80歲的風險 → 圖表與健檢報告
```

第一段只需要定序資料，適合在三級分析的後處理階段先完成。第二段需要知道「從幾歲開始評估」及「目前是否已經患病」，所以放在個案的 PRS 頁面，由使用者補充資料後按「風險計算」。

六個疾病會共用讀取基因資料、計分與祖源檢查的工作，避免把同一份全基因體資料重複處理六次。不過，每個疾病仍有自己的計分規則、參考分數分布與風險參數。

### 目前使用的疾病模型

「模型」在這裡指一套已發表的計分規則：它指定要看哪些基因位置，以及每個位置對分數的影響方向與大小。平台使用 PGS Catalog 資料庫收錄的 INTERVENE_MegaPRS 系列，沒有自行重新訓練這些基因權重。

| 疾病 | 使用的 PGS 模型 |
|---|---|
| 冠狀動脈心臟病（CHD） | [PGS004879](https://www.pgscatalog.org/score/PGS004879/) |
| 第二型糖尿病（T2D） | [PGS004870](https://www.pgscatalog.org/score/PGS004870/) |
| 全部癌症（All cancers） | [PGS004875](https://www.pgscatalog.org/score/PGS004875/) |
| 乳癌 | [PGS004869](https://www.pgscatalog.org/score/PGS004869/) |
| 心房顫動（AF） | [PGS004878](https://www.pgscatalog.org/score/PGS004878/) |
| 癲癇 | [PGS004881](https://www.pgscatalog.org/score/PGS004881/) |

<details>
<summary>生資方法與程式對照：各疾病如何建立自己的計分位點與參考分布</summary>

這是一次性的模型準備，之後個案分析沿用準備好的版本。**新增五疾病各自對完整reference比對；沒有把它們限制在CHD位點內。** 共用處理的是六模型位點的聯集，並非六病共同交集。

1. `scripts/setup_prs_multi.py`是準備入口，下載並核對官方GRCh38 harmonized權重的SHA256，來源記在`research/prs/multidisease/sources.json`。`prepare.py`的`score_rows()`目前只接受常染色體、ACGT雙等位SNV及加性權重，不接受indel、交互作用或隱性計分規則。
2. `research/prs/multidisease/prepare.py`的`prepare()`逐一讀取完整HGDP＋1kGP v1的GRCh38 PVAR，以座標及allele pair匹配每個疾病模型。`candidates()`優先直接匹配；反向互補候選須經dbSNP157確認身分。相同優先級有多個候選時排除，不按列順序任選。CHD另核對凍結版master的ID、REF、ALT及effect allele，維持610,305位點。
3. 五個新模型分別要求匹配數／原始權重列數≥99.5%，各病實際固定數量記入`manifest.json`，不直接把官方列數當成可用數量。SQLite的`locus`表存聯集位置，`model`表保留疾病所屬、effect allele與權重，因此一次讀基因型仍能分別計分。
4. 使用相同的3,330名無明顯親緣參考者，以PLINK 2 `--freq`估頻率、六欄`--score`建立參考分數；每病各自取其中735名EAS的分布。CHD全3,330人須與舊Model A分數的最大差≤`1e-5`，且CHD正式分布保留原值；其餘五病使用新算的各病分布。完成這些檢查才寫出模型manifest。

**審查重點：** allele harmonization與不唯一匹配的排除是否充分、被排除位點是否集中於重要權重，以及使用同一批參考者是否適合各疾病。CHD回歸通過只表示CHD校準保留，不能當作五個新模型的疾病預測驗證。

</details>

## 三、從原始基因資料算出分數

### 1. 先判斷每個位置是否讀得可靠

gVCF 是定序分析產生的基因資料檔。除了記錄偵測到的變異，它也保留許多未偵測到變異區域的資訊，讓平台有機會區分「這裡確實沒有該變異」和「這裡沒有足夠資料可以判斷」。這個區分很重要：資料不足的位置不能直接當成沒有風險變異。

平台會依模型指定的位置讀取個案資料，檢查定序品質，並確認基因位置及變異方向能與模型對上。可信的結果才用來計分；低品質、未涵蓋或無法明確對應的位置，會保留為缺失。資料轉換前後也會再核對一次，避免轉換格式時改變基因型。

<details>
<summary>生資方法與程式對照：gVCF重建、品質門檻與格式往返核對</summary>

`research/prs/multidisease/reconstruction.py`的`scan_chromosome()`用bcftools按染色體串流讀取gVCF，並以samtools `faidx`查GRCh38 FASTA。`run.py`的`reconstruct()`處理22條染色體，每條只掃一次六病聯集。

- `reconstruction.py`的`rules()`先核對凍結版vendor檔案的checksum，再載入其中Python重建函式；規則原檔是`research/prs/scoring_v211/vendor/PGS004879_FINAL_gVCF_resume.sh`中的`source_record()`、`observe_at()`與`choose_observation()`。六病流程重用這些規則，沒有把整支舊shell流程執行六次。
- 目前要求單一樣本、可判讀二倍體GT、`FILTER=PASS`、DP≥10及GQ≥20。reference block須為0/0，深度優先取`MIN_DP`，缺少時才用DP；沒有DP／GQ、低品質或無覆蓋者記missing。`FILTER=.`也不當成PASS。
- 多等位紀錄只投影目標SNV；GT含非目標ALT時不把它折成REF。會處理同座標多筆及重疊紀錄；個案FASTA的REF若與已準備reference不符，不另外猜測翻股。A/T、C/G回文SNP仍採`assume-direct`假設，沒有藉此解決所有股向歧義。
- `research/prs/multidisease/run.py`的`run()`用bcftools串接各染色體VCF，經PLINK 2匯入PGEN再匯出VCF；`roundtrip()`逐位點核對位置、REF／ALT、缺失及ALT dosage。這驗證格式轉換保留資料，不是對原始定序基因型的獨立真值驗證。

**審查重點：** 不同gVCF caller的block、FILTER及品質欄位是否符合這些假設，GRCh38／FASTA是否一致，以及回文SNP與多等位排除對各病分數的影響。

</details>

### 2. 按照模型，把各位置的影響加總

在這裡分析的常染色體上，每個位置通常有兩份DNA，一份來自父親、一份來自母親。同一位置的DNA可能有不同版本，稱為等位基因。平台先計算個案帶有幾份模型指定的版本，再乘上模型給這個位置的權重，最後加總所有位置的貢獻。有些權重增加分數，有些降低分數，不能把所有變異都當成危險因素。

計分時仍會遇到少量無法可靠判讀的位置。為了讓每個個案都能在同一組位置上計算，平台用參考資料中該位置的平均值補入。這樣做可以維持分數的可比性，但補入的是人群平均，並沒有找回這個人的真實基因型。因此，缺失太多時，即使算得出分數，也不會允許繼續估算風險。

目前參考資料來自 HGDP 與 1000 Genomes 兩個人類基因組研究計畫。為避免親屬造成重複影響，先篩選彼此沒有明顯親緣關係的參考者，固定使用3,330人估計補值所需的平均值；下一步比較百分位，則使用其中735名東亞參考者。前者用來處理缺失，後者用來比較排名。

平台表格的「原始 PRS score」指的就是**固定計分位置、包含上述平均補值的總分**。它可能與完整工具紀錄裡其他計分階段的數字不同，因為那些階段採用的位置集合或缺失處理方式不同；閱讀個案結果時，以正式表格為準。

<details>
<summary>生資方法與程式對照：加權計分、平均補值與獨立對帳</summary>

計分由`research/prs/multidisease/run.py`的`run()`及`summarize()`執行。對每個疾病，`beta`是原模型權重，`d_effect`是可靠位點的effect-allele dosage（0／1／2），`p_effect`是固定3,330人的effect-allele頻率：

```text
score = Σ_observed(beta × d_effect) + Σ_missing(beta × 2 × p_effect)
```

- PLINK 2以`--read-freq reference.afreq`及`--score-col-nums 3-8`一次計六欄。為共用ALT dosage矩陣，若某病effect allele為REF，準備階段將該項權重改為負號，另加`2 × beta`常數offset，恢復原effect-allele定義。
- Python另按各病原effect allele逐項計算，將observed與imputed貢獻分開加總；與PLINK校正offset後的總分差須≤`1e-5`。CHD補值頻率沿用凍結master，其餘疾病由新參考模型準備。
- 結果存`summary.json`、`scores.tsv`及`missing_loci.tsv.gz`，保留各病可靠／缺失位點、補值貢獻、較高權重缺失及分組邊界提醒。模型檔由`research/prs/multidisease/common.py`的`load_models()`核對manifest與SHA256，禁止改檔後默默沿用舊校準。

**審查重點：** 補值使用3,330人的跨族群平均，而排名使用735名EAS；應評估這項選擇、較高權重缺失及個案缺失模式對排名的影響。平均補值不是定相或單倍型imputation，也不是個案真實GT。

</details>

## 四、把分數放進參考族群，得到百分位

單看一個分數，例如0.36，並不知道它在人群中算高還是低。平台因此把個案的分數，放進735名東亞參考者的同疾病分數中排序。

舉例來說，假設有588名參考者的分數低於或等於某個案，588除以735就是80%，因此這個案位於第80百分位。這個計算只告訴我們相對排名，還沒有用到疾病發生率，所以無法從這一步得知發病機率。

選擇東亞參考者，是為了讓比較族群較接近目前預期的使用對象。平台也會檢查個案的基因資料是否適合與這批參考者比較，這就是祖源分析的用途。若無法明確判定適合，或該疾病的可靠計分位置不足，就不會開放風險計算。

因此，**畫面可以先出現分數和百分位，風險則要等祖源與品質檢查通過**。這是流程分段完成的正常狀況。還要留意，735人是有限的東亞參考資料，並不是台灣人口常模；畫面顯示的小數位也不代表排名具有同等程度的精密度。

<details>
<summary>生資方法與程式對照：ECDF、共用祖源分析與風險放行條件</summary>

**百分位：** `research/prs/multidisease/run.py`的`summarize()`讀各病`reference_scores.json`，以`bisect_right()`計算`100 × count(reference score ≤ target score) / 735`，包含同分者。沒有用常態分布Z值換算排名；`eas_z`只是另外保存的描述值。風險分組使用未四捨五入的百分位，邊界歸右側組。

**祖源：** 同檔`ancestry()`從聯集PGEN擷取`ancestry.ids`；這份清單只由固定CHD位點產生，先以`--geno 0`去掉該個案缺失位置，再交給`pgscatalog/pgsc_calc v2.3.0`、Nextflow及Apptainer。這只共用祖源判斷，不會把五個疾病的分數或參考分布換成CHD。

- 參考panel為HGDP＋1kGP v1；現有範例設定使用FRAPOSA OADP投影、`SuperPop`標籤、`pca_maf_target=0`及`pca_geno_miss_target=0.1`，實際值由該job的config決定。PCA最終使用經QC／LD pruning的子集，不是把全部610,305點直接當獨立標記。
- `research/prs/multidisease/ancestry_input.py`的`stage_score()`只補完整PGS metadata、逐列保留變異與權重，再由`preflight()`呼叫pgscatalog-format檢查；預設格式工具映像為pgscatalog-utils 2.0.4。
- `research/prs/absolute_risk/multi.py`的`ancestry_ready()`要求樣本相符、`MostSimilarPop=EAS`、非低信心、可接受QC狀態及正的PCA標記數；`score_rows()`另要求每病observed／fixed loci≥99.5%。這是個案門檻，與準備模型時的99.5% reference匹配率不同。
- 距風險分組邊界≤2百分點、缺失位點的`abs(weight)≥0.001`會記錄提醒；提醒本身不等於自動判定失敗。這些研究用門檻沒有經台灣疾病結果驗證。

**審查重點：** 以CHD位置選出的祖源panel是否有選位偏差、能否充分辨識混合或非EAS個案，以及低信心旗標是否足夠。這不是獨立全基因組祖源分析，仍是可再改進的部分。

</details>

## 五、加入疾病統計，估算未來發病風險

### 1. 需要兩種互補的資料

知道某人位於第80百分位，仍然不知道這個疾病有多常見。平台需要同時取得「分數與發病的關係」及「台灣人群的疾病基準」，才能把排名轉成機率。

- **INTERVENE 提供分數分組與發病的關係。** 本版採用原研究愛沙尼亞生物資料庫（EstBB）的結果，估計不同分數組別相對於中間組的發病速率。研究中的 HR（hazard ratio）就是這種相對發病速率的倍數。例如HR為2，表示在研究設定下，發病速率約為參考組的兩倍，但不表示10年發病機率必然乘以二。
- **GBD 提供台灣疾病基準。** GBD是全球疾病負擔研究。本版沿用作者公開的2019年資料，依疾病、年齡及性別，取得台灣的新發病例、既有病例及死亡統計。新發病例用來建立發病基準，既有病例用來調整尚未患病者的估計，死亡資料則用來考慮可能先因其他原因死亡的情況。

這兩種資料各自補上一部分資訊：台灣基準告訴我們疾病的常見程度，INTERVENE則告訴我們如何依分數分組調整這個基準。**這仍是把歐洲研究結果用到台灣的推估，並不等於已在台灣驗證。**

<details>
<summary>生資方法與程式對照：HR、GBD資料整理與疾病對應</summary>

風險參數由`research/prs/absolute_risk/build_multi.py`的`build()`／`bins()`離線建立，沒有在每次風險計算時重新下載或重新擬合。`registry.py`記疾病、PGS、GBD cause、適用性別與HR分層；原始來源的SHA256固定在builder的`PINS`，來源URL與完整清單見`research/prs/releases/tw_multi_v1/README.md`。

| 疾病 | GBD2019 cause | EstBB HR選擇 |
|---|---|---|
| CHD | Ischemic heart disease | 男性age＋sex；女性sex |
| T2D | Diabetes mellitus type 2 | age＋sex |
| 全部癌症 | Total cancers | age，不另按性別選HR |
| 乳癌 | Breast cancer | full sample HR，僅計算女性風險 |
| AF | Atrial fibrillation and flutter | age＋sex |
| 癲癇 | Idiopathic epilepsy | age，不另按性別選HR |

- **HR處理：** 使用EstBB原始Age＋Sex、Female、Age及Full Sample四種表。對age／age＋sex模型，每個非參考分組以四個年齡層的`MedianAAO`為x、`log(HR)`為y、`1/SE²`為權重做線性回歸，再取各五歲區間中點的HR；超出觀察年齡時固定在端點，不無限外推。此函式為`research/prs/absolute_risk/build_taiwan.py`的`male_hr()`；名稱來自舊CHD版本，現由六病builder重用。sex／full模型直接用對應HR，40–60組固定為1。
- **GBD處理：** 讀作者的`Sex_Stratified_{Incidence,Prevalence,Mortality}.csv`；乳癌另讀`BreastCancerProstateCancer_*.csv`。只選台灣、2019年、`metric=Rate`及指定cause／All causes，依年齡、性別鍵值連接。每十萬人的率除以100,000，發生／死亡為每人每年，盛行則為比例；不使用All ages、年齡標準化率或直接以病例數代替率。
- **資料例外：** 全部癌症使用`CANCER_SOURCE_GROUPS`反向對應七個來源分組標籤，不倒數或重新擬合HR；AF四格20–29歲疾病死亡率僅依已稽核的作者補零例外接受。原始六份GBD CSV相關台灣500列保存在`taiwan_raw_gbd_audit.json`，builder逐格核對正式包使用的基準率，其他缺列拒絕。

**審查重點：** HR的SE與年齡回歸定義、端點固定的合理性、疾病endpoint是否與GBD cause相容，以及全部癌症標籤對應是否有足夠證據。重現發表值可支持實作，但無法單獨證明這些外推假設成立。

</details>

### 2. 從目前年齡開始，逐段計算

平台先依百分位，把個案分進 INTERVENE 使用的七組：低於20、20–40、40–60、60–80、80–90、90–95及95以上。40–60是比較用的中間組。各疾病採用原研究選定的年齡、性別分層方式，並非所有疾病都套相同的相對風險。

接著，平台從評估年齡開始，將台灣人口發生率依這些相對風險分配到各組，同時維持各組加權平均與人口基準一致，再逐段累計到指定年齡。計算過程也會考慮：有些人可能先因其他原因死亡，而不再有機會發生這個疾病。這稱為「競爭風險」，也是不能單純把百分位、人口發生率與年數相乘的原因。

目前使用五歲年齡區間的統計資料，假設同一區間內的發病及死亡速率固定，並假設未來仍維持2019年的年齡別情況。這是模型的簡化，沒有預測醫療進步或未來發病率的變動。同一分組、年齡與性別會得到相同估計；分組邊界附近的小幅排名變化，可能影響風險估計，應留意品質提醒。

七個分組也會帶來階梯式的差異：跨組時風險數字可能跳動，來源估計的相鄰組也不一定依序增加。因此不能把每一小段百分位增加，都解讀成風險必然上升。平台會依實際評估期間提示非單調情況，接近組界的提醒則不代表已估出個人風險的誤差區間。

<details>
<summary>生資方法與程式對照：競爭風險公式與數值方法</summary>

`research/prs/absolute_risk/engine.py`的`group()`決定七分組，`integrate()`累計條件式風險。令`I`為人口發生率、`P`為盛行比例、`HR_g`為目標組相對風險、`q_g`為分組比例；所有量均取該年齡區間及性別的資料：

```text
q = (0.20, 0.20, 0.20, 0.20, 0.10, 0.05, 0.05)
H_bar   = Σ(q_g × HR_g)
lambda_g = I × HR_g / H_bar / (1 − P)
mu       = all_cause_mortality − disease_mortality

deltaF = S × lambda_g/(lambda_g+mu) × [1 − exp(−(lambda_g+mu) × dt)]
S_new  = S × exp(−(lambda_g+mu) × dt)
```

`S`是尚存活且未發病的機率，從評估年齡設為1；各段`deltaF`加總成累積風險，`dt`是實際跨到的區間長度。總hazard為零時增量為零。程式用`math.expm1()`計算小增量，採`competing_exact_v1`在每段固定hazard下積分，可切分非五歲邊界的起訖年齡。

原文重現另用`intervene_discrete_v1`：每個完整五年區間的增量為`S × [1 − exp(−5 × lambda_g)]`，S仍以疾病與其他死因hazard更新。這兩條數值路徑沒有混用，也不能把原文重現的誤差當成UI算法與原文的差距。

**審查重點：** 各PRS組共用盛行率與其他原因死亡率、固定分組比例、五歲區間內速率固定，以及將人口率轉成未患病者hazard的假設。百分位不依年齡／性別重新排名，年齡／性別調整發生在風險參數。歐洲HR能否套用到固定子集、平均補值後的EAS分組，仍是未完成台灣追蹤驗證的關鍵。

</details>

### 3. 三種風險的期間與適用條件

若在50歲評估，未來10年指50至60歲，未來20年指50至70歲，目前至80歲則指50至80歲。三者都從同一個評估年齡開始；20年風險會重新累計，不能用10年風險乘以二。

| 顯示的風險 | 可提供的評估年齡 |
|---|---|
| 未來10年 | 20至70歲 |
| 未來20年 | 20至60歲 |
| 目前至80歲 | 20至未滿80歲 |

本版計算到80歲為止，所以超出終點的完整10年或20年會顯示不適用。所有估計都以「目前仍存活、尚未罹患該病」為前提，用來估算首次發病；已有疾病者不計算復發風險。乳癌風險模型目前只適用女性。

在平台操作時，使用者核對年齡、性別及各疾病病史後，按「風險計算」即可保存評估。病史預設「否」只是方便輸入，仍須確認。修改臨床資訊後要重新計算；有有效評估的個案，才可把PRS納入健檢報告。

<details>
<summary>生資方法與程式對照：執行入口、來源綁定與風險輸出</summary>

平台worker是`backend/app/workers/prs_run.py`的`Worker.run_multi()`，依序呼叫`research/prs/multidisease_cli.py --stage score`與`--stage ancestry`；CLI轉入`research/prs/multidisease/run.py`。計分完成先發布summary，祖源完成再發布QC，不以祖源工具的輔助分數覆蓋正式分數。

風險API由`backend/app/routers/prs.py`接收請求，經`backend/app/services/prs.py`分派至`prs_multi.py`的`assess()`。`research/prs/absolute_risk/multi.py`的`person_for()`驗臨床輸入，`evaluate()`逐病檢查病史、性別、計分及祖源條件，再呼叫`engine.integrate()`計算10年、20年、至80歲與40–60參考組曲線。API不採用使用者自行傳入的score／percentile；超出完整期間的欄位回null，病史／性別／品質不符的疾病另附原因。

`prs_multi.py`固定核對`research/prs/releases/tw_multi_v1/bundle.json`的release hash，並把臨床輸入、score、bundle、來源hash及結果保存為不可變評估。`report_assessment()`依當前資料重現結果並比較快照；健檢匯出由`backend/app/services/docx_export.py`使用同份有效評估，不另算一套風險。個案輸入或來源變更時，舊快照保留但不能靜默當成新結果。

</details>

## 六、開發時遇到的問題，如何影響目前的選擇

### 固定計分位置：從 Method 1 到 Method 2

開發初期，工具會依每個個案實際可用的位置，讓該個案與參考者使用同一批位置計分。由於個案之間可用的位置不同，同一名參考者在不同個案的分析裡，也會得到不同的分數。這表示每次排名的比較標準都可能改變，使不同個案的結果不容易放在一起比較。

**Method 1**先取三個開發樣本與參考資料共同可用的607,445個位置。固定這個集合後，同一名參考者在不同個案分析中算出的分數就一致了，支持「位置集合不同是差異來源」的判斷。不過，這個集合是依這三案的資料選出來的，適合作為問題驗證，仍不足以作為長期處理新個案的標準。

因此後來採用**Method 2**：在分析新個案之前，先依官方模型與參考資料的相容性，固定CHD的610,305個位置。新個案沿用這個集合，少數缺失依共同規則補值，參考分布也維持固定。現行六疾病流程延續這個原則，並將能共用的處理整合在一起。

### 比較缺失資料的影響：Model A 與 Model B

Method 2決定「固定位置再計分」的原則，A／B則是在這個原則下，選擇「固定哪些位置」的兩個版本。Model A保留610,305個相容位置；Model B進一步排除參考資料較常缺失的位置，只保留603,438個。這項比較是想了解「參考資料的缺失，會不會明顯改變排名」。

| 歷史驗證樣本 | Model A百分位 | Model B百分位 | 分組是否改變 |
|---|---:|---:|---|
| 1 | 35.1 | 35.4 | 否，皆為20–40組 |
| 2 | 60.5 | 61.5 | 否，皆為60–80組 |
| 3 | 79.5 | 79.5 | 否，皆為60–80組 |

三案的排名差異都約在1個百分點以內，沒有改變分組，但第三案接近80的邊界，仍需留意。最後採用A作主分析，以保留較多相容位置及既有校準；B保留作敏感度比較，也就是檢查換一種合理處理方式後，結果是否穩定。這些結果支持這三案的技術一致性，沒有證明A在台灣的疾病預測能力優於B。

上述名稱都是**計分方法的開發比較**，不是讓使用者在不同疾病風險算法中挑選。詳細分數、位置規則與原始紀錄列在文末的延伸資料。

## 七、目前驗證到哪裡

驗證需要分成三個層次，因為「程式算得一致」與「能準確預測台灣人的疾病」需要不同的證據。

**第一層是計分是否正確且穩定。** 開發時用兩套程式獨立加總，確認結果在允許誤差內一致。歷史紀錄中，三個CHD樣本從原始gVCF重新分析，A／B共6組結果全部重現；兩個樣本同時執行時，4組結果也一致。六疾病共用流程另以合成資料測試。這些證據用來檢查計分與執行流程，不能當作疾病預測能力的證明。

**第二層是能否重現原研究。** 使用原作者的相對風險與愛沙尼亞疾病基準，核對論文發表的1,232筆風險結果；其中乳癌112筆，其餘五病各224筆，全部誤差都小於0.005個百分點，符合原表兩位小數的精度。台灣基準資料也逐格核對來源。這表示資料讀取與研究計算可以重現，仍不是台灣實際發病驗證。

原研究以五年為單位近似累計；平台為了支援任意評估年齡，採用在每個區間內同時計算發病與其他原因死亡的方式。兩種算法分開保留：前者用來重現原文，後者用於畫面估計。它們不會每個數字都完全相同，既有CHD比較測試的最大差約1.44個百分點，測試範圍見延伸紀錄。這個數字描述的是兩種算法的差距；對台灣實際發病的預測準確性，仍要靠下一層驗證。

**第三層是台灣人的實際疾病結果，目前尚未完成。** 要確認這一層，必須取得台灣個案的診斷日期、追蹤期間與死亡資料，檢查預測的發病人數是否接近實際人數。前兩層通過，不能取代這一層。

執行可靠性也仍有待補驗的項目。已測過完成工作的續跑與部分失敗情境，但不能宣稱所有分析中途的中斷都已驗證可恢復；六疾病完整院內WGS流程與大量個案負載，也仍需持續驗收。

### 外部 review 後，哪些疑問已經有查核結果？

2026-11-14 我們針對 review 的部分疑問，重新讀取固定版本的來源檔與作者程式。以下是已取得的結果；它們的用途是縮小問題範圍，不是宣告台灣疾病風險已驗證。

| 查核問題 | 已取得的結果 | 結果能說明什麼？ |
|---|---|---|
| 原始權重中，有沒有容易混淆方向的 A/T、C/G 回文位點？ | CHD 610,677、T2D 930,506、全部癌症 818,872、乳癌 872,060、AF 785,779、癲癇 605,432 列，六個來源的回文位點計數皆為 **0**。 | 對這六份原始權重，不需要為了來源回文位點另算一個排除版本；其他等位基因對應及個案讀取仍需核對。這些列數不是個案實際可靠位點數。 |
| 這次補查的 CHD 來源，是否就是當初固定的那份？ | 原始完整權重解壓後的 SHA256 與原 model lock **完全相符**。 | 確認來源檔版本一致；院內 610,305 位點 master 尚未重新讀取，不能把來源一致說成院內模型已重新驗收。 |
| 全部癌症是否在 GRCh37→GRCh38 harmonization 時被整體翻轉？ | PGS004875 原始與 harmonized 檔的 818,872 列，rsID、effect allele、other allele、weight 逐列比較，**差異為 0 列**。 | 這兩個檔案之間沒有整體變號；仍不能確認作者產生 EstBB HR 表時使用了哪個分數方向。 |
| 作者公開程式有沒有說全部癌症要取負號？ | 固定版本的計分、染色體加總與年齡分組模板，未找到全部癌症專用整體取負號；分組依分數升序建立。 | 公開模板不足以補齊作者當次實際分析紀錄。反向對應 HR 標籤仍是能重現 224 筆曲線結果的推論，沒有改動平台 raw score，也沒有把 HR 取倒數。 |
| AF 年輕年齡的零值，應怎麼解讀？ | 來源在 20–29 歲的發生與盛行資料列為零，部分疾病死亡率另由作者補零；呈現邏輯已按完整評估期間檢查。 | 全期間落在來源零值區間時，顯示「來源資料在此年齡設為零，無法據此確認零風險」；例如25→35歲跨入非零區間，仍照實計算完整期間。零值不是低於偵測極限的證據。 |

逐筆來源證據在 [source_audit.json](validation_20261114/source_audit.json)，作者程式與 HR 對照在[執行紀錄第二節](PRS_VALIDATION_2026-11-14.md#二全部癌症方向查核證據與決策)。來源回文為零與全部癌症方向，是兩個不同問題，不能因前者通過就把後者結案。

獨立計分與補值比較工具已提供，但**院內保留模型、標準品基準比較、真實缺失遮罩與 EAS／台灣分布比較尚無本輪實測結果**。已核對標準品名單為 HG001 對應 VAL-45～49、HG002 對應 VAL-50～54，共兩個基準個體、十次檢測；這是資料對照，不是十人驗證通過。健康人盤點工具也只建立檔案清單，不會因此產生 PRS 或品質結果。[本輪設計與狀態](PRS_VALIDATION_STUDY_2026-11-14.md)

<details>
<summary>生資方法與程式對照：驗證資料與可重現範圍</summary>

| 想核對的項目 | 程式／紀錄 | 證據範圍 |
|---|---|---|
| CHD固定A／B與三案／平行執行 | `research/prs/scoring_v211/METHOD2_INTEGRATION_SOP.md`、`VALIDATION_REPORT_v211_2026-11-02.md`（同目錄） | 歷史交付結果；沒有因本次補寫文件重跑院內WGS |
| 六病gVCF重建、計分與往返 | `tests/test_prs_multi_tools.py` | 小型合成22染色體資料，含真正bcftools／samtools／PLINK的測試；缺工具時該測試會skip，不能把skip算通過 |
| 風險解析解、期間與輸入條件 | `research/prs/tests/test_absolute_risk.py`、`tests/test_prs_multi.py` | 數值與條件邏輯驗證，不是疾病預測效能 |
| 原文1,232筆點估計 | `research/prs/absolute_risk/build_multi.py`、`research/prs/releases/tw_multi_v1/reproduction.json` | 使用原HR與Estonia基準；公開表格只作expected值，不作係數擬合 |

CHD歷史A／B計分由`research/prs/scoring_v211/method2_score.py`讀取已固定的模型。官方原始610,677位點中，A保留610,305個reference相容位置；B再要求EAS及全體參考者的位點缺失率均≤1%，保留603,438個。這些原始篩選證據見歷史SOP與model lock；大型master／score／reference資料在設定的`method2_model_dir`，未隨repo附上，審查時須另取得。六病`prepare.py`只沿用CHD A並做回歸，沒有替五疾病產生B模型。

builder要求每病完整覆蓋兩性／適用性別、七組及16個年齡終點，最大誤差≤0.0050001百分點才產生正式參數。重現作者出生起曲線時，第一格沿作者將1–4歲資料對應五年區間的約定；正式台灣包僅從20歲起，不因此提供出生風險。原作者R bootstrap未在此流程執行，個人不確定性區間仍為null。

Review時請同時查看原文離散重現與UI連續積分的測試，不把「點估計重現」「流程能跑完」「台灣疾病預測準確」視為同一個驗證結論。

</details>

## 八、閱讀結果時，最需要注意的限制

**東亞排名加上台灣統計，仍然包含跨族群推估。** 主要模型與相對風險來自歐洲研究，735名東亞參考者則用來比較分數。兩者與台灣真實發病關係是否相符，仍需前述追蹤驗證。此外，平台目前沒有納入血壓、血脂、血糖、吸菸、家族史、罕見致病變異或治療效果，也尚未提供個人風險估計的誤差範圍，因此這不是完整的個人臨床風險評估。

**疾病名稱相近，資料中的定義仍可能不同。** CHD原研究包含不穩定心絞痛、心肌梗塞、其他急慢性缺血性心臟病及心跳停止；GBD基準則按急性心肌梗塞、慢性缺血性心臟病等估計，所以不能把輸出只理解為心肌梗塞或心血管死亡風險。其他疾病也有對應上的限制：全部癌症原研究的診斷碼沒有涵蓋所有血液惡性腫瘤；AF基準包含心房撲動；癲癇則使用GBD「特發性癲癇」基準。

另外有兩項資料問題需要保留說明。全部癌症的原始相對風險表，其高低分組標籤與論文曲線呈現的方向不一致；開發時調整標籤對應後，重現了224筆結果，但這仍是依公開資料推得的重建方式，尚需作者確認。AF在20–29歲的原始發生與盛行資料為零，部分疾病死亡資料由作者補零，因此這段年齡的模型零值不能解讀為不會發病。

本次已確認全部癌症原始／harmonized 權重沒有整體翻轉符號，但作者公開模板仍不足以連結實際 EstBB 輸入與該 HR 表，不能據此宣稱方向已驗證。AF 則逐一看評估期間：若全段位於來源零值區間，顯示「來源零值」並說明無法確認零風險；若例如25→35歲跨入非零區間，保留整段的計算數值並附前段限制。這不是將零值說成低於偵測極限。

目前六個原始權重檔均已核對沒有 A/T 或 C/G 回文位點；CHD 來源也符合原 lock 的解壓後 hash。這只排除來源的特定疑慮，院內保留模型、公開基準樣本及真實缺失型態仍須另用 `research/prs/validation.py` 核對。[方法與指令](../../ops/PRS_VALIDATION.md)

最後，六個疾病的風險各自估算，不能相加成「總風險」。尤其全部癌症與乳癌涵蓋的事件會重疊。比較曲線時也要讀取百分比刻度，因為各圖的縱軸範圍可能不同。

## 九、後續如何更貼近台灣使用者

後續有兩個不同的改進方向：一個是讓排名更能代表台灣，另一個是確認風險預測是否符合台灣的實際發病情況。

針對排名，約1,200名成人健檢／社區個案已經可以開始做**沒有疾病追蹤也能進行的技術與分布比較**。先盤點實際檔案、去除重複個體、確認祖源與計分品質，再把六病分數放進現行 EAS735 參考，觀察排名是否有系統性偏移。檔名26T0開頭目前只是使用者指定的健康人候選規則，不等於已完成臨床病史確認；有效人數也要等盤點與QC後才知道。

**同日收到的院內盤點清單已有實際數量：1,516個26T0候選檢體ID、1,516個gVCF路徑，分布於27個run目錄（目錄日期2026-06-11至2026-11-09）。** 每個ID在回傳清單中均只有一份、回報為非空且名稱相符的檔案，沒有需要二選一的來源；這尚不是去親緣／去重後人數或品質通過數。另外10份標準品均找到，HG001與HG002各5份，都在2026-04-28同一個run，現階段可研究同run重複性，不能稱為跨run重複性驗證。使用者另回傳終端摘要，確認本次搜尋完成、讀取錯誤與未追蹤目錄連結皆為0，名稱衝突與跨ID共用實體也皆為0；header及PRS結果仍未知。只保存[不含個案明細的盤點摘要](validation_20261114/cohort_inventory_summary.json)，原始清單留院內。

Review 的 M4 提議檢查百分位均勻性，M5 則提議用真實缺失模式量測補值影響。這兩項現在就能規劃進行，但要先分清楚「分數分布」與「排名分布」：**原始分數不要求均勻；如果台灣個案與 EAS 的分數分布相近，轉成 EAS 百分位後，各排名區間的人數才會大致接近。** 例如20–40組約占20%、最高5%約占5%，不是要求每個組別人數一樣。

| 先做的比較 | 看什麼結果 | 如何解讀 |
|---|---|---|
| 六病分數與 EAS735 並列 | 平均、標準差、中位數、尾部、分布／Q–Q圖 | 看整體位置、寬度或尾端是否不同；差異可能來自族群、收案或技術，不能直接判哪個一定錯。 |
| 現行 EAS 百分位的分布 | 平均排名、各十分位與七個風險組的比例 | 當成 EAS 參考是否貼近這批使用者的線索；分組比例的理想參照為20／20／20／20／10／5／5。 |
| 按技術批次與品質分層 | 各批次分數／缺失率及高權重缺失，與祖源PC的關係 | 先排除某一批次整體偏高、偏低或缺失較多的情況。年齡／性別資料未齊時，這些分層先留待補資料。 |
| 候選台灣排名與現行排名並列 | 同一人的百分位位移、跨組與尾端穩定度 | 了解換參考會影響多少結果；此階段不替換平台參考或覆寫舊報告。 |

百分位來自有限735人的ECDF，包含同分與離散階梯。統計比較應保留這些特性及參考抽樣的不確定性，可從兩組原始分數的分布差與適當重抽樣著手；不能只把百分位當連續均勻分布，憑單一KS檢定的P值決定採用。[KS兩樣本檢定的假設](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ks_2samp.html)

**用同一批人建立台灣ECDF，再觀察這批人的台灣百分位近似均勻，是排名定義造成的結果，不能當成外部驗證。** 候選常模需用保留／後續新樣本或按家庭與批次的重抽樣檢查穩定度。健康人資料可改善「跟誰比較」，卻不能單憑分布判定全部癌症的疾病方向，也不能估計台灣HR或證明發病機率準確。

本輪先盤點26T0候選與兩種標準品，並進行技術／分布及補值的資料準備。癌症病例對照因年齡、性別及癌別資料尚未齊全，依使用者決定暫緩；全部癌症方向繼續保留待確認。實際執行入口與資料需求見[健康人盤點及比較操作](../../ops/PRS_COHORT_INVENTORY.md)。未來若另取得合適的台灣參考，例如TWB，亦可按同一套方法與現行參考並列評估。

針對風險，則需要台灣的疾病追蹤驗證：除了確認高分組是否較容易發病，也要檢查模型預測10%風險的一群人，實際發病比例是否接近10%。若有偏差，就需要評估相對風險或台灣基準是否須調整。這項驗證可以先針對現行版本進行，不必等到換成台灣的百分位參考資料。

更後續可評估加入臨床因素、更新人口統計年度與呈現估計的不確定性，並完成更大規模的流程測試。每一種更新都需要保留版本與重新驗證，才能知道新舊結果的差異來自哪裡。目前TWB等台灣個人資料尚未納入這套模型。

## 延伸資料：工具、來源與詳細驗證

基因資料主要由bcftools與samtools讀取、檢查，PLINK 2負責計分，Python另做獨立加總確認。祖源分析使用pgsc_calc v2.3.0與FRAPOSA；Nextflow安排分析步驟，Apptainer提供工具執行環境。風險估算由平台的Python程式執行。這些工具讓流程可重複執行，前面介紹的研究假設與限制仍需另外驗證。

### 程式完整路徑與入口索引

本文的`research/prs/...`、`scripts/...`及`backend/...`都從**repo根目錄**起算。依既有checkout，本機根目錄是`/Users/changym/Desktop/NGS-UI/`，院內伺服器是`/home/n102968/NGS_UI/NGS-UI/`；若部署位置不同，請換成實際repo根目錄。例如個案運算入口的完整路徑為：

```text
本機：/Users/changym/Desktop/NGS-UI/research/prs/multidisease_cli.py
院內：/home/n102968/NGS_UI/NGS-UI/research/prs/multidisease_cli.py
```

| 要審查的步驟 | Repo內程式路徑 | 主要入口／函式 |
|---|---|---|
| 一次性模型準備 | `scripts/setup_prs_multi.py` → `research/prs/multidisease/prepare.py` | `main()` → `prepare()`、`candidates()` |
| 健康人與標準品gVCF盤點 | `scripts/inventory_prs_gvcfs.py` | `main()`、`scan()`；只讀檔案metadata，不讀gVCF內容或計分 |
| 來源／獨立計分／補值敏感度 | `research/prs/validation.py` | `source_audit()`、`model_audit()`、`benchmark()`、`aggregate()`；真實資料結果仍待執行 |
| 新個案計分與祖源 | `research/prs/multidisease_cli.py` → `research/prs/multidisease/run.py` | `main()` → `run(stage='score'/'ancestry')` |
| gVCF判讀規則 | `research/prs/multidisease/reconstruction.py` | `rules()`、`scan_chromosome()`、`scan_records()` |
| 格式往返／分數／百分位／祖源 | `research/prs/multidisease/run.py` | `roundtrip()`、`summarize()`、`ancestry()` |
| 祖源輸入格式預檢 | `research/prs/multidisease/ancestry_input.py` | `stage_score()`、`preflight()` |
| 建立風險包／原文重現 | `research/prs/absolute_risk/build_multi.py`、`research/prs/absolute_risk/build_taiwan.py` | `build()`、`bins()`、`male_hr()` |
| 疾病／PGS／GBD／HR對應 | `research/prs/absolute_risk/registry.py` | `DISEASES` |
| 個案風險條件與數值積分 | `research/prs/absolute_risk/multi.py`、`research/prs/absolute_risk/engine.py` | `person_for()`、`evaluate()`、`group()`、`integrate()` |
| 平台啟動與發布 | `backend/app/services/prs_jobs.py`、`backend/app/workers/prs_run.py` | `start()`、`Worker.run_multi()`、`Worker.publish()` |
| 風險保存與匯出核對 | `backend/app/services/prs_multi.py`、`backend/app/services/docx_export.py` | `assess()`、`report_assessment()`、`build_health_docx()` |

未啟用六病模型的舊CHD個案，仍可能經`research/prs/scoring_v211/pipeline_v21.py`與`method2_score.py`執行。它們也是追溯CHD歷史A／B的入口；不能把舊CHD腳本誤當成新增五疾病的個案入口。上述科學程式可對照[本次核對的repo版本](https://github.com/puddingyd/NGS-UI/tree/114281b/research/prs)閱讀。

### 設定與分析產物在哪裡

PRS設定由`NGS_UI_PRS_CONFIG`指定，未設時讀`NGS_UI_HOME/data/prs_config.json`；`multi_model_dir`決定是否啟用六病。工具／FASTA／dbSNP／reference panel等實際路徑均由此設定取得。現有PLINK容器檔名標示`2.00a5.10`；完整環境重現還須一併查看該job的設定、容器身分與工具紀錄。

| 資料 | 預設位置或設定 |
|---|---|
| 準備完成的六病reference模型 | `NGS_UI_HOME/biotools/prs/multi_v1/`，由`multi_model_dir`指定；含manifest、SQLite、頻率與各病參考分布 |
| 模型準備暫存 | `NGS_UI_HOME/nf_work/prs/_reference_prepare_v1/` |
| 個案重建、audit、Nextflow及暫存 | `NGS_UI_HOME/nf_work/prs/{個案ID}/`；沿既有work root設定及SSD掛載，容器tmp在該案`.runtime/tmp/` |
| 正式分數／QC／風險評估 | `NGS_UI_TERTIARY_ROOT/{個案ID}/08_postprocessing/`；分數連結為`{個案ID}.prs.summary.json`，發布版本與快照在`{個案ID}.prs/` |
| 固定台灣風險包與驗算紀錄 | repo的`research/prs/releases/tw_multi_v1/` |

院內典型的`NGS_UI_HOME`是`/home/n102968/NGS_UI`，三級主root預設是`/home/datalake_Intermediate/pipeline/tertiary_output`；環境變數可覆寫，因此查實際job時以保存的config及來源紀錄為準。舊PRS專案只供讀取既有資源，不接收新個案輸出。`research/prs/scoring_v211/config.example.json`保留歷史配置，不能把其中舊`project`路徑直接當成平台目前的個案輸出位置；平台路徑配置實作在`backend/app/services/prs_paths.py`。

### 研究來源與歷史紀錄

- [INTERVENE原論文與附件](https://doi.org/10.1038/s41467-024-48938-2)：模型架構、疾病定義與原文風險結果。
- [作者公開的GBD資料與風險程式](https://github.com/intervene-EU-H2020/flagship/tree/179986db2f9af3b6eb32e4e16dc9201a23d2a29c/AbsoluteRiskEstimation)：本平台使用固定版本，以維持結果可追溯。
- [Method 1／2與A／B詳細方法](https://github.com/puddingyd/NGS-UI/blob/7ef35a9/research/prs/scoring_v211/METHOD2_INTEGRATION_SOP.md)、[CHD歷史驗證紀錄](https://github.com/puddingyd/NGS-UI/blob/7ef35a9/research/prs/scoring_v211/VALIDATION_REPORT_v211_2026-11-02.md)、[兩種風險算法的比較紀錄](https://github.com/puddingyd/NGS-UI/blob/7ef35a9/research/prs/VALIDATION.md)：供需要查看細節者查閱，不代表已完成台灣疾病驗證。
- [六疾病參數與資料來源](https://github.com/puddingyd/NGS-UI/blob/7ef35a9/research/prs/releases/tw_multi_v1/README.md)、[1,232筆逐筆重現紀錄](https://github.com/puddingyd/NGS-UI/blob/7ef35a9/research/prs/releases/tw_multi_v1/reproduction.json)：保存來源版本與驗算結果。
