# PRS 健康人與標準品：院內驗證操作

v10.37，2026-11-14。接續[已完成的盤點](PRS_COHORT_INVENTORY.md)：1,516個26T0健康人候選檢體ID、27個run，以及HG001／HG002各5份標準品。這份工具把既有六疾病計分接成可接續的院內批次，先回答「重複檢測是否一致、缺失補值影響多大、這批人的分數與EAS參考差多少」。它不計算新的疾病HR或更換平台常模；癌症病例對照仍暫緩。

## 一、先理解三種不同的比較

| 比較 | 怎麼做 | 可以回答什麼 |
|---|---|---|
| 標準品重複性 | 同一HG001或HG002的5份gVCF各自計分；比較分數、百分位、風險組與共同讀到位置的基因型 | 同一基準個體重跑是否穩定。十份皆在4/28同run，不能稱跨run重複性，也不是與GIAB答案比較。 |
| 補值頻率敏感度 | 保留個案實際讀到的基因型，只把缺失位置的平均補值從現行3,330人改成EAS735人平均 | 換平均值會移動多少分數／百分位、是否跨組；沒有已知答案，不能據此判哪種補值較準。 |
| 參考者遮罩試驗 | 從EAS參考基因型匯出單一參考者，先算可確認的原始分數；遮住真實個案讀不到的位置，分別補兩種平均 | 在有答案的位置上衡量補值誤差。參考者自己缺失的位置仍標為未知；不完整真值不報全模型排名。 |

健康人分布比較另外使用**現行EAS735排名**，看六病分數的位置、寬度、十分位、七組比例及run差異。分數不要求均勻；EAS排名接近均勻也不能證明風險預測準確。人數仍是檢體ID數，未做個人身分去重或親緣分析。

## 二、建立這一輪的固定輸入

在院內repo更新後，用原PRS執行帳號、既有`.venv`及工具環境執行，不需重啟服務：

```bash
cd ~/NGS_UI/NGS-UI
.venv/bin/python scripts/run_prs_cohort_validation.py prepare \
  --inventory-dir /home/n102968/NGS_UI/nf_work/prs/_validation_20261114/inventory_20261114T070046Z_6hjyoj6n \
  --study /home/n102968/NGS_UI/nf_work/prs/_validation_20261114/cohort_v1
```

`prepare`需要原目錄內的`inventory.json`及`samples.tsv`，不能只給下載回本機的TSV。它會交叉核對清單、六病固定來源權重與模型、gVCF單樣本名稱、22條常染色體的GRCh38長度、GT／DP／GQ與參考區塊宣告、相鄰索引及caller header；每50份顯示進度。這只是前置檢查，不讀全基因型、不證明樣本身分或coverage，也不替原始gVCF建立索引。單份header失敗會記錄於`plan.json`；後續若選到失敗樣本則停止，避免默默排除。

設定沿用`NGS_UI_PRS_CONFIG`，未設時讀`NGS_UI_HOME/data/prs_config.json`；若終端與服務的環境不同，可在上面命令加`--config /實際路徑/prs_config.json`。須已啟用`multi_model_dir`，並能執行現有的bcftools、samtools、Apptainer及Nextflow。保留DP≥10、GQ≥20及現行股向政策；不重新準備或啟用模型。

原始權重預設讀既有PRS work下`_reference_prepare_v1/scorefiles/`，CHD讀設定中的凍結來源；若來源權重另存，請用`prepare --score-dir /實際來源目錄`指定。缺少來源會停止，不下載別的版本替代。EAS頻率與參考試驗還需要模型目錄保留的`reference.pgen/.pvar/.psam`及兩份keep清單。

成功後先看`cohort_v1/preflight_summary.json`；逐樣本錯誤看`plan.json`。只有前置核對通過才進入下一步。

## 三、先跑十份標準品

```bash
.venv/bin/python scripts/run_prs_cohort_validation.py run \
  --study /home/n102968/NGS_UI/nf_work/prs/_validation_20261114/cohort_v1 \
  --group standards --stage score

.venv/bin/python scripts/run_prs_cohort_validation.py report \
  --study /home/n102968/NGS_UI/nf_work/prs/_validation_20261114/cohort_v1 \
  --genotypes
```

每份直接呼叫平台現行`multidisease_cli.py`，共用六病位點重建、PLINK 2與Python對帳、固定EAS735百分位，沒有另寫一套計分。`report`只讀校驗和及執行後核對均通過的結果；看到中途產生的`summary.json`不會當成完成。

報告放在新建的`report_.../`。`summary.json`的`standard_score_repeatability`逐病列重複份數、run數、分數／百分位範圍及跨組情況；`standard_genotype_repeatability`比較兩兩共同觀測位置的基因型，另外列單側缺失與雙側缺失。觀測一致率不把缺失算成一致；只有一份完成時不報零差異冒充重複性通過。

**目前這份批次工具尚未整合GIAB高信賴區間與真值的外部比較。** 標準品結果只能先解讀為重複性；GIAB可靠範圍、位點覆蓋與外部誤差需另接資料。原[獨立VCF驗證工具](PRS_VALIDATION.md)也不能把未列出的VCF位置當0/0。

## 四、每個run先取一份健康人候選

```bash
.venv/bin/python scripts/run_prs_cohort_validation.py run \
  --study /home/n102968/NGS_UI/nf_work/prs/_validation_20261114/cohort_v1 \
  --group healthy --per-run 1 --stage score

.venv/bin/python scripts/run_prs_cohort_validation.py report \
  --study /home/n102968/NGS_UI/nf_work/prs/_validation_20261114/cohort_v1
```

本次清單預計選27份，每個run以固定ID雜湊排序選取，不按分數挑人。先看run分數、缺失比例及高權重缺失有沒有異常；每run一份不能估計完整批次分布。預設只計分，祖源狀態仍是未確認。

若要確認這27份的祖源，另執行相同選樣的`run --group healthy --per-run 1 --stage ancestry`，其他參數與上面相同，再產生一次報告。此階段需每份各跑既有Nextflow祖源流程，會重複處理共用參考，耗時與磁碟需求較高；目前沒有新增整批祖源共用快取。也可以一開始用`--stage both`順序完成計分與祖源。標準品不要求判成EAS。

`summary.json`會分開列全部已計分者的**初步描述**，與已確認EAS且逐病觀測位點達99.5%者的比較。沒跑祖源時後者為0，不當成沒有找到健康人。各run的初步分布仍可用來找技術問題，但不能作為已篩選的台灣常模。

## 五、計算兩種平均補值的差異

```bash
.venv/bin/python scripts/run_prs_cohort_validation.py frequency \
  --study /home/n102968/NGS_UI/nf_work/prs/_validation_20261114/cohort_v1
```

用保留的3,330人參考PGEN及EAS735 keep重新估EAS等位基因頻率，核對成員／股向／觀測數，再對每個已完成個案的實際缺失清單計算差異：`Σ 缺失位置權重 × 2 ×（EAS頻率 − 現行頻率）`，依effect allele處理REF與ALT。比較時固定現行EAS ECDF，才不會把「補值變了」與「排名參考也變了」混在一起。

每次產生新的`frequency_.../`；`summary.json`列逐病最大百分位位移與跨組數，`sensitivity.tsv`含個案明細、兩種分數、百分位位移與缺失數。標準品及健康人候選分組保留。這是敏感度試驗，尚不能宣稱EAS補值比較正確，也不會套用到UI。

接著可用同一個`frequency_...`做有參考答案的遮罩試驗。以下`本次frequency目錄`需換成上一步終端印出的完整路徑：

```bash
.venv/bin/python scripts/run_prs_cohort_validation.py reference \
  --study /home/n102968/NGS_UI/nf_work/prs/_validation_20261114/cohort_v1 \
  --frequency-dir /本次frequency目錄 \
  --output /home/n102968/NGS_UI/nf_work/prs/_validation_20261114/cohort_v1/reference_pilot_v1 \
  --reference-limit 10 --mask-limit 10
```

先固定選10名EAS參考者、最多10份已完成樣本的缺失模式；有多run時優先分散取樣，未完成的健康人不會提供缺失模式。總共最多100組比較，並逐組保存結果；這是小規模試跑，不代表已評估全部735人或全部個案。對完整真值疾病，獨立未遮罩分數須與凍結參考分數差≤1e−5；不符即停。原始參考基因型缺失時保留`PARTIAL_TRUTH_COVERAGE`，不強填成答案、不報全模型百分位。結果列完整／部分比較數與兩種補值誤差，不能把部分比較當完整確效。

## 六、試跑確認後擴到1,516份

```bash
.venv/bin/python scripts/run_prs_cohort_validation.py run \
  --study /home/n102968/NGS_UI/nf_work/prs/_validation_20261114/cohort_v1 \
  --group healthy --all-healthy --stage score
```

沿用相同模型及已完成檢體的checkpoint，完成後重新產生`report`與`frequency`。不需年齡、性別或病史即可完成這些技術比較；目前不啟動風險計算。全體祖源可再用`--all-healthy --stage ancestry`，但先根據27份試跑的耗時／磁碟與品質判斷，不自動接著啟動。

本輪沒有實作身分／親緣去重、候選台灣ECDF、重抽樣信賴區間、Q–Q圖或臨床效能驗證。報告先提供可追溯的分布數字；不得把1,516個檢體ID稱為1,516名無親緣、品質合格的健康人。

## 七、執行位置、接續與交付

- 所有新增產物在既有`NGS_UI_PRS_WORK_ROOT/_validation_20261114/cohort_v1/`；上面命令沿用已確認的`/home/n102968/NGS_UI/nf_work/prs/`邏輯路徑及既有掛載。自訂root時同步換`--study`及`--output`，不要另開pipeline頂層目錄。
- 每案工作、容器`/tmp`／`/var/tmp`、Nextflow及工具暫存在`cases/{檢體ID}/`下，沿用平台容器隔離規則；只讀原始gVCF、模型與來源。這批研究結果不登錄新個案、不發布到08、不覆寫UI結果。
- 與平台PRS共用執行名額，一次一份；平台有工作時等待，每案完成後釋出名額。不需停止平台。整批耗時由試跑的`state.json.elapsed_seconds`觀察，不能拿單一PLINK計分秒數估全流程。
- `Ctrl-C`或SIGTERM會終止本工具啟動的程序組並保存checkpoint。第一個錯誤即停；看該案`state.json`與`tools.log`，處理後用相同命令接續。來源／設定／模型／科學程式變更時拒絕混跑，請用新的study名稱；文件更新不影響接續。
- 若`prepare`尚未寫出`plan.json`就失敗或被中斷，保留該目錄供查錯，改用新的study名稱重新prepare，不宣稱可從半份前置清單接續。參考遮罩的選樣或參數改變時也需新的`--output`。
- `plan.json`、`scores.tsv`、`sensitivity.tsv`、缺失清單與原始基因資料留院內。先交回`preflight_summary.json`，以及各`report_.../summary.json`、`frequency_.../summary.json`、`reference_pilot_v1/summary.json`供解讀；摘要不列健康人個別ID。

程式入口是`scripts/run_prs_cohort_validation.py`，實作在`research/prs/cohort_validation.py`；方法與歷史證據看[PRS workflow](../research/prs/PRS_WORKFLOW.md)。本版已做本機合成及既有工具整合測試，尚未在院內對這批WGS執行，沒有新的真實補值誤差或臨床通過結論。

## 八、本機驗證紀錄

本輪相關測試130項通過、7項未執行。包含來源／header與索引變動拒絕、研究輸出與SSD根目錄連結、重複性及缺失分母、兩種頻率的REF／ALT方向、參考者未遮罩分數／遮罩比較／checkpoint遭改拒絕，以及既有PRS工作、QC、模型準備與真實bcftools／samtools／PLINK的小型六病整合測試。新批次的頻率與參考匯出測試使用合成資料及模擬工具輸出，不能代替院內全模型容器驗收。未執行者為6項需指定實際Nextflow的選配檢查及1項pgscatalog格式套件檢查。

另核對平台說明下載：未登入401、登入200且內容與新版workflow逐byte相同、不讀個案、檔案缺失503。文件連結與差異格式檢查通過；沒有連入院內主機或重啟服務。
