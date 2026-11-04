# CHD PRS 研究風險 UI — v9.85

## 本版範圍

平台採兩階段；第一個疾病為 CHD，未建立 DM／cancer 的空白模型。

1. **三級分析**勾選「PRS」：00–08 發布完成後，自動排入原始 gVCF → 固定 v2.1.1 Model A/B score／percentile → 祖源 QC 的背景任務。不需要病人年齡或病史，不計算 absolute risk。PRS 啟動或執行失敗不回滾已完成三級結果。
2. **Secondary findings → PRS**：上方單行輸入「評估時年齡、性別、已有 CHD？」後按「風險計算」。日期自動記為當日；每次保存獨立評估快照。改臨床輸入只重算 risk，不重跑 gVCF。

主表每列一個疾病，固定欄位為「疾病／原始 PRS score／PRS percentile／未來10年風險／目前至80歲累積風險」。未計算的兩個風險欄顯示「請輸入臨床資訊後計算」；曲線放表格下方。未安裝 risk bundle 或祖源尚未通過時，合格 score／percentile 仍可獨立顯示，風險按鈕停用並說明原因。

表格的原始 score 指 **Model A total_sum（含固定模型缺失補值）**，與該列 EAS735 ECDF percentile 使用同一定義；不是先前610,677 loci的 observed-only SUM，也不是 pgsc_calc ancestry 的 normalized score。A/B完整結果仍保留於計分 output，B 暫不作UI風險主分析。

### 為何 Nextflow 還在跑，表格已經有數字？（v9.85）

平台先執行 raw gVCF 重建與固定610,305位點計分，將完成的分數／固定735名東亞參考者百分位發布至08；接著才執行另一個祖源／品質檢查Nextflow。第二階段不取代表格分數，完成後只是補上檢查證據。第一階段完成就能查看score與percentile；祖源未確認符合東亞比較或QC未符合要求，仍不能計算風險。重跑期間已有舊結果時，明示「上一次已發布」，不把它當成新run的結果。

以26T00969的既有驗證紀錄及使用者Log為例：

| 數值 | 定義與用途 |
|---|---|
| 0.474838562769280／60.5442176871 | 固定位點加總（包含參考平均基因型補值）及固定EAS ECDF；表格顯示0.474839／第60.54百分位，風險分組使用這個未四捨五入的百分位 |
| raw_sum 0.471545 | raw重建中實際可用位點的加總，不是固定計分補值後總分 |
| normalized_sum 0.475192／percentile 60.6802721088 | 祖源pgsc_calc依該案匹配後的位點集合，另算target與reference分數／百分位；欄名normalized_sum不代表表格總分的換算，不用於UI風險 |

v9.85的PRS主Log保留Nextflow進度、完成時間、警告與錯誤，僅從同job的有效發布summary顯示「表格採用的結果」，並保留補值、分組邊界與QC提醒。原工具的JSON及COMPLETED SUM/PCT不再混入摘要；原始stdout/stderr逐byte保留，可由「下載完整工具Log」取得。新API `GET /api/samples/{sid}/prs/jobs/{job}/raw-log`沿用登入保護並核對case/job，只提供固定log.txt快照。舊Log也可讀取摘要，不重寫原檔。舊job若已被另一run取代，不能把新score貼到舊job下。

### 畫面上的方法與限制（v9.85）

`prs_presentation.explanation()`提供純顯示文字，獨立於已鎖定bundle及歷史評估；UI不再展示內部Model A／model_id，也不從舊snapshot載入難懂的limitations。說明區位於曲線下、分析工具區上方，即使尚未輸入臨床資訊也可閱讀；Log下的重複註腳已移除。科學算法、builder、bundle checksum、score／risk快照完全不變，部署重載即可。

方法分三步：PGS004879（INTERVENE_MegaPRS_CHD）固定610,305位點加權及缺值補值 → HGDP／1000 Genomes固定735名東亞參考排名 → INTERVENE EstBB七組相對風險搭配台灣GBD2019年齡／性別率及競爭死亡，估算首次CHD。未納入血壓、血脂、吸菸、糖尿病與治療效果，也未估計不確定性範圍。

疾病定義具體依據：原研究Supplementary Data 8（`MOESM11`，第9列）CHD ICD10為`I200|I21|I22|I23|I24|I25|I46`，含不穩定心絞痛、心肌梗塞及心跳停止等；ICD9亦含`413`。Data 9（`MOESM12`）另列`I20.1/I20.8/I20.9/I46`等不完全重疊的碼，UKB的重疊指標82%不能當成台灣準確率。[GBD2019的IHD定義](https://pubmed.ncbi.nlm.nih.gov/34922374/)則包含急性MI與慢性IHD（穩定心絞痛、MI後無症狀IHD）。因此不能將這個首次CHD研究估計宣稱為單獨MI、心跳停止或心血管死亡的預測。畫面提供[PGS004879](https://www.pgscatalog.org/score/PGS004879/)、[原研究與附件](https://pmc.ncbi.nlm.nih.gov/articles/PMC11169548/)及GBD定義連結。

風險適用20至79歲、確認尚未有 CHD、具備合格 EAS 祖源證據與相符 model lock 的個案。>70歲不將截短期間冒稱10年風險，仍可計算至80歲。以「評估時存活且無 CHD」為條件，不是出生累積率相減。TWB／BBJ未納入；目前不收血壓、血脂、吸菸或家族史，因為本模型不使用這些欄位。結果不加入診斷／健檢 DOCX 或列印報告。

## 執行環境與首次啟用

程式本身已接上 queue／worker；**Git pull 不會安裝大型計分 reference，也不會自動複製原伺服器的 config**。使用原 v2.1.1 已驗證的設定檔，讓服務讀取：

- `NGS_UI_PRS_CONFIG`：既有完整 config.json 的絕對路徑；未設定時讀 `$NGS_UI_HOME/data/prs_config.json`。欄位見 `research/prs/scoring_v211/config.example.json`，勿直接沿用不存在的示範路徑。
- 原 config 的 `project`／`raw_work` 保留作來源紀錄，**新平台工作不使用它們作輸出位置**。平台產生獨立 config snapshot，覆寫這兩欄；不要修改原專案 config 或重新 freeze model lock。
- `NGS_UI_PRS_WORK_ROOT`：大型 raw／Nextflow 工作根目錄，預設為既有 `NGS_UI_TERTIARY_NF_WORK_ROOT` 下的 `prs/`，一般部署為 `/home/n102968/NGS_UI/nf_work/prs/`。三級 work root 未覆寫時是 `$NGS_UI_HOME/nf_work/`；PRS會跟隨其設定，不在 `/home/datalake_Intermediate/pipeline/` 下新增資料夾。服務帳號需有建立個案目錄的權限。不能指到原 PRS 專案內或包含參考資源的目錄。若先前已明寫 `NGS_UI_PRS_WORK_ROOT`，需移除或改成新位置才會採用此預設。
- scorefile、model lock、Model A/B reference、FASTA、dbSNP、containers、兩份 Nextflow config 與 ancestry reference 維持既有有效路徑，僅讀取；不搬動 model lock 綁定的 model/source 絕對路徑。
- 平台 worker 的 PATH 需可執行 `nextflow, apptainer, bcftools, samtools`。`raw_script` bytes 必須等於 repo 中已驗證的 vendor 腳本；既有十檔核心 checksum 仍由 pipeline_v21 驗證。
- `NGS_UI_PRS_BUNDLE`：第二階段的風險參數包，安裝方式見下節；第一階段不依賴它。

在服務相同使用者／環境中檢查（只讀，不執行個案分析）：

```bash
PYTHONPATH=backend python3 -c 'from app.services.prs_jobs import settings; c,h=settings(); print("PRS settings OK",h)'
```

三級視窗「PRS」預設未勾選；勾選時先檢查基本執行環境，缺設定立即回報，不會讓操作者等完整三級後才知道缺設定。每個樣本在三級發布後核對 gVCF；WES、找不到來源或多個候選均記錄到三級 Log／`prs_jobs` 子任務摘要，其他三級結果保留。載入個案後可在 PRS 區補跑或確認 gVCF。

### 原始 gVCF 定位

綁定 `pipeline_source.json` 中的 source sample、原 VCF anchor、來源類型與新保存的 `seq_type`，不從 UI 顯示姓名或 LIS 字串猜檔案。僅WGS；舊 marker 缺類型時可讀已登錄 metadata 的 test_type，仍不明則拒絕。

DRAGEN 預設精確搜尋同一 anchor 目錄或同 run 的 `gvcf.gz/` 中 `{source}.hard-filtered.gvcf.gz`／`{source}.gvcf.gz`，以及同 run 的 `other/{source}/germline_seq/{source}.hard-filtered.gvcf.gz`（v9.78新增院內DRAGEN實際佈局）。不同實體檔案不自動選最新；多候選由操作者在 UI 下拉確認。單候選自動使用。

In-house 沒有猜測的預設 caller／目錄；院內管理者確認原始gVCF後，在同一 config 加入 `gvcf_templates`。這個欄位會取代預設 mapping，需保留所需的 DRAGEN templates；若已自訂，新增預設路徑不會自動合併進自訂值，須自行加入所需模板。只支援精確模板，不接受 glob 或 browser 提供任意檔案路徑：

```json
{
  "gvcf_templates": {
    "dragen": [
      "{anchor_dir}/{source}.hard-filtered.gvcf.gz",
      "{run_dir}/gvcf.gz/{source}.hard-filtered.gvcf.gz",
      "{run_dir}/other/{source}/germline_seq/{source}.hard-filtered.gvcf.gz"
    ],
    "inhouse": []
  }
}
```

`anchor_dir` 是原 SNV VCF 的目錄，`run_dir` 是它的上一層，`source` 是原 sequencing ID。In-house 只有核對 caller、組裝及實際命名後才能填入模板。每個選定輸入都要求單樣本 header 精確符合 source ID、GRCh38 chr1長度248956422、gVCF reference-block allele宣告，以及 tbi/csi index。正式執行另由原 pipeline／bcftools 驗證資料。

### 排程、重試與保存

- 工作狀態／Log：`$NGS_UI_HOME/data/jobs/prs/{jobUUID}/`；保存平台衍生 config snapshot、source／gVCF／resource fingerprints、路徑版本、操作者與 parent tertiary job。`NGS_UI_DATA_ROOT` 自訂時跟隨該位置。
- raw／Method2：`{NGS_UI_PRS_WORK_ROOT}/{個案ID}/integrated_runs/{runID}/{source}/`；個案ID與UI一致，例如 `26T00970-dragen`。這裡的 raw 包含重建 VCF、roundtrip、稽核、計分 results，**不只是 Nextflow 暫存**；真正 task 暫存為其下 `raw/nextflow_work/`。
- ancestry：相同 raw 完成後，以原 pipeline.py 的 reuse 模式執行於 `ui_prs_{UUID}_ancestry/`。**不使用它的 normalized percentile 覆蓋 Method2 percentile**。原始 Method2 summary 不改寫；平台 envelope 另記 QC hash 並附上同源 ancestry evidence。
- 正式結果：原子複製完整 Method2 六檔至 `08_postprocessing/{個案ID}.prs/runs/{jobUUID}/score/`，祖源 QC 四檔至同 run 的 `ancestry/`；各自保存 SHA256 manifest。UI入口 `{個案ID}.prs.summary.json` 記錄 artifact paths、hash、工作路徑版本與衍生 config hash。缺檔或驗證失敗不替換上一份有效結果。
- Nextflow repository／plugin／container cache、工具臨時檔放在個案工作目錄的 `.runtime/`；runner 的工作目錄在 `launch/{runID}/`。平台產生的 Nextflow wrapper 先讀原 config，再覆寫 cache／report 等執行路徑（包含 singularity profile），原設定檔不改寫。模型、篩選、計分及祖源參數不變。
- 每次先原子發布已完成 score／percentile，再跑祖源。祖源失敗時顯示「計分完成，祖源檢查未完成」；score 保留，risk 不啟用。非 EAS／低信心亦不啟用風險。
- 全站同時只跑一個 PRS 計分工作，其他排隊。與三級／Somatic共用 UI及source sample locks，防止同個案來源替換。三級完成不代表其獨立 PRS 子工作已完成；在 PRS 區看進度／Log。
- 「重新分析／重試」在輸入、資源與平台路徑 fingerprints 完全相同時沿用原 run ID，交給 v2.1.1 guarded resume；不刪除 raw audit。原 pipeline 不能恢復的 early partial run 可用「重新開始」建立新 run ID，原檔保留。jobUUID 每次不同，retry 可沿用 runID。
- 取消由 worker 收到 cancel flag 後終止自己啟動的工具 process group；原有 score／快照保留。完整刪除個案 pipeline 時，active PRS 阻擋刪除；完成工作狀態、08結果及有平台 ownership marker 的該個案 `nf_work/prs/{個案ID}/` 一起移除。取消登錄不刪這些分析產物。舊位置資料不自動清除。

預設佈局（以 `26T00970-dragen` 為例）：

```text
/home/datalake_Intermediate/pipeline/tertiary_output/26T00970-dragen/08_postprocessing/
  26T00970-dragen.prs.summary.json            UI讀取入口
  26T00970-dragen.prs/
    runs/{jobUUID}/score/                    A/B score、percentile、缺失清單、狀態及報告
    runs/{jobUUID}/ancestry/                 QC JSON/TSV、報告與高權重缺失清單
    assessments/{assessmentUUID}.json       臨床輸入、風險與比較曲線快照

/home/n102968/NGS_UI/nf_work/                 沿用三級Nextflow暫存根目錄
  contexts/shared-dragen/{launch,work}/      既有DRAGEN三級暫存
  contexts/shared-nckuh/{launch,work}/       既有NCKUH三級暫存
  prs/26T00970-dragen/
    .ngs-ui-prs.json                         平台個案 ownership marker
    .runtime/                               config wrappers、Nextflow及container caches、tmp
    launch/{runID}/{score,ancestry}/         runner工作目錄
    integrated_runs/{runID}/26T00970/
      raw/                                  重建與audit/results；nextflow_work/為task暫存
      method2/                              上游完整原始計分結果
    integrated_runs/{runID}_ancestry/26T00970/
      ancestry/{inputs,results,work}/
      qc/
```

正式 artifacts 以一般檔案複製，不留指向 work 的 symlink；風險計算只需08已發布結果與固定風險 bundle。工作目錄保留供稽核與 guarded resume。三級既有「Nextflow暫存清理」會**明確跳過 `nf_work/prs/`**，並保護自訂PRS根目錄所在單位；PRS即使完成也不把 raw audit／接續分析資料當一般task暫存刪除。PRS worker與直接runner沿用三級清理保護lease（scope=prs）；可辨識的PRS執行不會因缺少DRAGEN／NCKUH的lineage pointer而讓全部閒置三級暫存被保留，仍維持未納管程序一律保守保留。本版未新增 PRS 暫存清理按鈕。

### 已啟動的舊路徑工作

部署不會搬移或改寫正在執行的工作，舊 worker 已讀取舊 config，仍會沿舊位置完成。若要立即停止該個案繼續寫舊位置，先在 PRS 頁面取消並等狀態變成「已取消」，再更新／重啟平台、選「重新開始」。**不要只在它執行途中改 config 或搬目錄**。新工作不接續舊 project 或v9.79預設 `pipeline/prs_work/` 的 run；現存 score／風險快照仍可讀，舊位置歷史資料保留，不自動刪除或遷移。完整刪除個案僅管理目前設定的PRS工作根目錄，舊位置須另行確認處理。

啟動新工作後，Log 應先列平台的「個案工作目錄」「正式結果」；raw 指令的 `--work` 與 `--out` 均應位於 `nf_work/prs/{個案ID}/`（未自訂根目錄時）。指令中的 scorefile／model／FASTA 等唯讀資源仍可指向原專案，這不是個案輸出。

Nextflow 路徑設定依官方 [config include／優先序](https://docs.seqera.io/nextflow/config)、[runtime環境變數](https://docs.seqera.io/nextflow/reference/env-vars)及 [Singularity cacheDir](https://docs.seqera.io/nextflow/reference/config/singularity)；仍需由操作者於院內執行實際 WGS 驗收。

### 容器暫存與根磁碟空間（v9.81）

祖源 `INTERSECT_VARIANTS` 使用 Python `tempfile` 建立排序中間檔。v9.80只設定主機端TMPDIR時，Nextflow的 `env - ... singularity exec --no-home` 不保證容器能看到 `.runtime/tmp`；指定目錄不可用時，Python可能改用容器 `/tmp`，而它通常連到主機根磁碟。曾觀察到 `/home` 尚餘4.4TB、inode充足，但根磁碟97%使用率，祖源出現 `OSError: [Errno 28] No space left on device`。這不是 `APPTAINERENV… preferred` 提示本身造成失敗。

v9.81固定主機暫存位置為：

```text
/home/n102968/NGS_UI/nf_work/prs/{個案ID}/.runtime/tmp/
```

如自訂PRS root則跟隨該案project。Nextflow raw／ancestry兩份衍生config會保留原 `singularity.runOptions`，追加以下掛載與 `--env`（含選用的singularity profile，避免重複追加）：

- 主機個案 `.runtime/tmp` → 容器同一絕對路徑。
- 主機個案 `.runtime/tmp` → 容器 `/tmp` 及 `/var/tmp`，讓未使用TMPDIR的工具也留在個案work。
- 容器 `TMPDIR`／`TMP`／`TEMP` 明確設為該個案暫存絕對路徑；不依賴Nextflow是否轉送主機環境。

原raw腳本直接呼叫Apptainer時，由worker的 `APPTAINER_BIND*`／`SINGULARITY_BIND*` 及相容的 `*ENV_TMP*` 設定套用相同掛載，保留管理者原有bind，不修改凍結腳本或原參考config。Apptainer支援的掛載與環境參數見 [bind paths](https://apptainer.org/docs/user/latest/bind_paths_and_mounts.html)、[runtime environment](https://apptainer.org/docs/user/latest/environment_and_metadata.html)。

計分前新增「容器暫存掛載與寫入檢查」：使用已設定PLINK image，在主機暫存目錄建立隨機marker，進容器核對同一絕對路徑、`/tmp`、`/var/tmp` 都可讀到該marker，並逐一實際建立／刪除probe檔。失敗立即停止該PRS工作，不啟動計分，也不覆蓋既有結果。Log應先看到個案暫存路徑與 `CONTAINER_TMP_PASS`；Nextflow的 `.command.run` 應包含上述三個 `--bind`。真正院內掛載由此執行時檢查確認；本機沒有執行Linux容器WGS。

本次衍生config的hash有變更；舊工作不原地改寫、不強制放寬guarded resume。部署後在PRS區「重新開始」，建立使用新設定的PRS run；可能重新執行raw，但不必重跑三級。舊run、已發布score與風險快照保留，失敗不覆蓋。程式不會自動清理主機 `/tmp` 或搬移舊暫存，避免影響其他工作。

### 三級自動 PRS 的環境交接與 Log（v9.84）

曾觀察到三級自動PRS的容器預檢失敗，但個案PRS頁面手動重試成功：三級worker設定的一次性 `contexts/shared-dragen/.runtime/jobs/{三級job}/tmp` 掛載，被傳入detached PRS。PRS等待個案鎖釋放後才執行，此時三級finally已清除該tmp，Apptainer因此以 `mount source ... doesn't exist` 結束。這與PRS自己的暫存目錄及模型資料無關。

平台現於建立PRS worker前就套用獨立個案環境，每個PRS stage再次核對同樣的掛載：保留reference與工具bind，移除所有屬於三級一次性jobs樹的source（無論它當時是否還存在），替換舊 `/tmp`、`/var/tmp`，依容器destination去重。使用BINDPATH相容變數，清空會重複添加的BIND變數；容器預檢已透過argv明寫的三個tmp掛載，不再由env加第二次。三級本身也使用相同的bind整理規則，避免site script與多次環境套用累積相同掛載。

不保留／重建已清理的父tmp，不放寬輸入或model lock檢查。這次只改環境交接，未改PRS科學config及既有run fingerprint；更新後可在個案PRS頁面重試失敗工作，原guarded resume規則照常生效，不需重跑三級。手動已成功的PRS結果可以沿用。

後處理主Log使用下列格式，stopgaps、GPN、CNV整合註解及PRS排程／執行共用步驟摘要：

```text
[2026-11-04 10:00:00] [postprocessing] sample=S1 step=genebe status=start
[2026-11-04 10:01:00] [postprocessing] sample=S1 step=genebe status=done elapsed=60.0s
```

AnnotSV配置、欄位清單、GPN signature差異與完整第三方輸出保留在 `$NGS_UI_HOME/data/jobs/tertiary/{三級job}/postprocessing.log`，主Log顯示狀態、耗時、必要數量與警告。三級進度視窗及分析清單的「查看Log」均可下載這份工具Log；新增登入保護的 `GET /api/dragen/jobs/{job}/postprocessing-log` 只提供該job固定檔案的當下快照。舊工作的原始log不重新改寫。

三級done代表00–08已發布，PRS仍為獨立背景工作。前端持續輪詢其狀態，API／三級Log會依parent job與sample識別核對子工作，顯示實際step、失敗／警告及FATAL等原因；PRS失敗不將三級改成failed，也不回滾已發布結果。手動補跑是另一個工作，不把舊子工作紀錄改成成功。PRS詳細工具Log仍由個案PRS頁面查看。

## 固定參數與科學驗算

- HR：INTERVENE 原研究 Estonian Biobank（EstBB）歐洲祖源，非跨 cohort meta-analysis。男性使用六個非參考 PRS 組各四個年齡組的 HR/SE/MedianAAO，按作者 `1/SE²` 加權線性回歸擬合 log-HR；外側年齡固定於最近組的**擬合值**。女性使用六組 sex-specific HR；參考組 HR=1。
- 基準：作者固定 commit 的台灣 GBD2019 IHD incidence/prevalence/all-cause mortality/IHD mortality。
- 原方法驗算：以 Estonia baseline 與作者五年離散算法，逐一比對 Supplementary Data 4 的 **224筆 CHD 結果，全部通過**。最大絕對誤差 `0.004995100317904511` 百分點，小於表格保留兩位小數所容許的 `0.0050001` 百分點。沒有重現 bootstrap CI，也沒有聲稱原 R 執行環境完全重建。
- UI 算法：`competing_exact_v1`，同一組 HR／基準率，採每個五年區間內固定疾病／死亡 hazard 的連續競爭風險積分，支援50.5歲等非五年邊界。**它不是與論文五年離散近似逐值相同的數值方法**。台灣兩性七組、20–75歲每5歲起點的10年／至80歲測試網格，兩方法最大差異為 `1.4436119070748987` 百分點；不同方法必須保留名稱，不得混稱。
- 成人參數包只含 [20,80)，未將來源的1–4歲改寫成0–4歲。驗算出生曲線時，才明記沿用作者第一格乘5年的約定。
- Model A、EAS735 參考分布及原始 model lock 固定；分位組跨族群與 A 子集的 outcome 相容性仍未驗證。原研究歐洲 HR 套用台灣率是研究外推，並非完成台灣臨床校正。

模型 ID：`intervene_estbb_tw2019_A_exact_v1`。`clinical_validation=NOT_VALIDATED`；`uncertainty_interval=null`。不生成吸菸／血壓／血脂或治療情境，也不套臨床高低風險門檻。

## 安裝資料（由操作者部署，不由 Codex SSH）

程式與文件隨 Git 更新；參數資料不放 Git。此次本機已產生：

```
research/prs/runtime/tw_chd_v1/bundle.json
research/prs/runtime/tw_chd_v1/reproduction.json
research/prs/runtime/tw_chd_v1.zip
```

將資料包解壓到伺服器 `$NGS_UI_HOME/biotools/prs/tw_chd_v1/`，或以 `NGS_UI_PRS_BUNDLE` 指定 `bundle.json` 的完整路徑。正式 production checkout 的 NGS_UI_HOME 通常是 `~/NGS_UI`。包內 `reproduction.json` 保留逐筆 expected/calculated/error 與來源 hash。

本版 bundle SHA256：

```
322d4e68d49a2047d7326172a609b3898b8754f658f800f50806f4e2276a0a27
```

後端固定檢查這個 digest 與模型版本；synthetic、被改寫或僅自行標記 PASS 的參數不能啟用。原檔 bytes 必須保留，不可重新排版 JSON。後續更新模型需另做驗算、版本與 release digest 更新。

讀取資料包檢查（在 repo 根目錄）：

```bash
PYTHONPATH=backend python3 -c 'from app.services.prs import load_bundle; b,h=load_bundle(); print(b["model_id"],h)'
```

## 相容入口：匯入既有計分結果（新流程不需要手動匯入）

```bash
PYTHONPATH=backend python3 scripts/import_prs_summary.py \
  --sample-id '實際已登錄個案ID' \
  --summary '/伺服器上既有PRS結果/summary.json'
```

`summary.sample` 必須符合 UI 個案 ID，或該個案 `pipeline_source.json` 內的 `source_sample_id`；不接受任意指定 alias。僅接受 `PGS004879_METHOD2_V0_20261031`、`METHOD2_TECHNICAL_PASS`。匯入保留原始 summary、原始 SHA256、來源識別與匯入時間，連結到當時 pipeline source marker 的 hash。

「成功匯入」不等於「適合風險計算」。UI 還會檢查 Model A、locus counts、percentile/group 一致性、model lock 與合格 EAS 祖源證據。原 pipeline 預設 raw-only 的 `NOT_ASSESSED` 無法估算；應取得既有流程的真實 ancestry QC，不得手改為 EAS 或偽造 PASS。

使用者自行更新／重啟平台後，載入個案，在「Secondary findings → PRS」填入年齡、性別與 CHD 病史後按「風險計算」。未安裝參數／未匯入 score 時，UI 明示缺件，不顯示合成數值。

## API 與保存

- `GET /api/samples/{sample_id}/prs`：登入保護；回傳 available、issues、score 摘要與有效的上次評估。
- `POST /api/samples/{sample_id}/prs/jobs`：啟動／重試；只接受伺服器列出的 candidate_id 與布林 new_run。
- `POST /api/samples/{sample_id}/prs/jobs/{job_id}/cancel`、`GET .../log`：個案綁定，禁止跨個案讀取／取消。
- `POST /api/samples/{sample_id}/prs/assessments`：只接受 `age, sex, assessment_date, prior_chd`。percentile、HR、參數路徑皆由伺服器讀取，不能由 request 指定。
- 若輸入缺失／不合格，不會寫評估快照。future date、NaN、bool age、未知性別、既往 CHD 及超出20–79歲皆拒絕。

依既有 state resolver 選擇 v3 的 `08_postprocessing/` 或 legacy state directory；新 PRS 檔名一律加個案前綴：

```
{LIS_ID}.prs.summary.json
{LIS_ID}.prs/runs/{jobUUID}/score/{summary.json,scores.tsv,report.md,status.json,inputs.lock.json,missing_loci.tsv.gz,manifest.json}
{LIS_ID}.prs/runs/{jobUUID}/ancestry/{qc.json,qc.tsv,report.md,missing_high_weight.tsv,manifest.json}
{LIS_ID}.prs/assessments/{UUID}.json
```

每次評估使用新 UUID 原子保存：完整 bundle／summary／person、輸入 hash、建立時間、操作者、兩個風險與兩條曲線。未寫入 SNV review／ACMG／report state。summary 或模型、個案 metadata 改變時，不再顯示舊快照為目前結果；pipeline source marker 或 worker 記錄的 gVCF／index signature 改變則要求重新分析。舊快照留作歷史稽核。前端改輸入即隱藏舊數值，切換個案會中止／丟棄舊 request，防止交叉顯示。

## 來源與重建

原 HR 來自作者 GitHub `MetaAnalysis/FullSample_MetaAnalysisandHetero_GoogleDrive.R` 引用的[原研究資料根目錄](https://drive.google.com/drive/folders/1bwedgU4lb4Y4i1pLsHXzjVaLBusjfcOQ)，其中 EstBB Full_cohort → Update_v3_211122：

- [男性年齡／性別 HR](https://drive.google.com/file/d/18ylwVsjMsaaQ-FASwHat_XuhTyItmI3i/view)：`HR_AgeandSexStratified_EstBB_211122.csv`。
- [女性 HR](https://drive.google.com/file/d/14IQVWEdxTvp2OloAji2637VDmI8VaCRC/view)：`HR_FemaleSample_EstBB.csv`。

同一資料根目錄的 MGB 候選檔與發表曲線未吻合，因此未採用；FinnGen 候選的年齡性別檔缺 Sex 欄，不猜測男女順序。EstBB 是本版唯一採用的 HR。

本機原始資料在 `research/prs/runtime/original_hr_20261102/`、GBD 與論文表在 `runtime/sources_20261102/`。builder 先檢查所有 input SHA256，再依精確 schema 解析；錯檔、缺組、驗算失敗皆不建立輸出包：

```bash
python3 -m research.prs.absolute_risk.build_taiwan \
  --sources research/prs/runtime/sources_20261102 \
  --hr-dir research/prs/runtime/original_hr_20261102 \
  --output-dir /tmp/tw_chd_v1_new
```

輸出目錄必須尚不存在。資料使用與再散布授權未因程式 MIT 授權而自動成立；原 HR／GBD 檔未提交 Git。

## 本機驗證

- 原文 Estonia CHD：224/224符合發表表格精度。
- v9.72 Python：PRS API／queue／worker／引擎及相關三級回歸，共66項通過；外部WGS工具以隔離假工具驗證發布流程，非真實WGS驗收。原版計分27項離線測試亦通過。
- 既有 frontend Node 回歸：93項通過。
- 獨立瀏覽器：桌面單行輸入／疾病表格、390px手機表格捲動、待計算文案、兩曲線、API數值、修改輸入清除舊值、>70歲、祖源不合格仍顯示score、開始／Log／取消、切換個案晚到回應，全部通過。預覽使用明示的合成個案，不是實際病人結果。
- 未執行院內 fresh WGS、未測試院內部署、未進行台灣 clinical validation。
- v9.79 路徑隔離／完整發布／重試／刪除及相關三級、risk 回歸：81項通過。測試以合成輸入與隔離假工具執行，不代表院內 Nextflow／WGS 已驗收；原 PRS 專案逐檔 bytes／mtime 不變、舊路徑不resume、祖源失敗保留score、缺完整artifact不覆蓋既有結果均已驗證。
- v9.80 將PRS工作納入既有三級 `nf_work/prs/`；相關回歸與Nextflow清理測試共114項通過，涵蓋預設／自訂根目錄、執行中與已完成PRS均保留、其他閒置三級暫存可清、runner實際繼承lease。未在院內執行WGS。
- v9.81 容器暫存／queue／worker／risk／三級與清理回歸共127項通過，原版計分27項離線測試通過，核心checksum未變。以真正Nextflow 24.10.4驗證空白／既有config及profile合併，並以模擬容器engine執行Nextflow產生的 `env -` launcher，確認Python暫存寫到個案目錄；這驗證指令產生與環境傳遞，不是實際Apptainer掛載或院內WGS驗收。
- v9.84 相關198項Python與94項前端測試通過，含父tmp清除後自動PRS完成、與手動環境一致、reference保留、目的地去重、AnnotSV raw保留與摘要、失敗根因、下載權限及PRS持續輪詢。實際Nextflow24.10.4的設定／resume／launcher測試仍以模擬engine執行，院內Linux容器與WGS待部署後驗收。
- v9.85 相關78項Python測試通過（5項未啟用的Nextflow外部工具測試略過），原版27項離線測試及94項前端單元測試通過；Playwright以合成資料驗證桌機／390px手機、計分已發布但QC未完成仍停用風險、重跑保留上次結果提示、方法與限制在評估前可見、無Model A用語、下載連結及既有曲線／輸入失效／切案競態。科學核心、bundle及快照未更動，未跑院內WGS。

可重複執行的測試入口：

```bash
PYTHONPATH=backend python -m pytest tests/test_prs.py tests/test_prs_jobs.py research/prs/tests/test_absolute_risk.py tests/test_dragen_same_directory_batch.py tests/test_dragen_rerun.py tests/test_dragen_delete_output.py -q
PYTHONPATH=backend python -m pytest tests/test_prs_container_tmp.py -q
# 設定 NGS_UI_TEST_NEXTFLOW=/已安裝的nextflow 可另執行真實Nextflow設定／launcher測試；未設定時這5項略過。
TMPDIR=/private/tmp python3 -m unittest discover -s research/prs/scoring_v211/tests -q
node --test tests/frontend_*.test.cjs
# 需可載入 Playwright；CHROME_EXECUTABLE_PATH 可指定已安裝 Chrome。
# 測試自建本機靜態預覽，所有API為合成資料，不連正式平台。
node tests/browser_prs_ui.cjs
```
